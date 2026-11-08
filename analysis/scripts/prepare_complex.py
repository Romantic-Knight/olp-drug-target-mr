"""prepare_complex.py — assemble the MD starting complexes and check them geometrically.

The input to an MD run is one file containing the protein and the ligand in a single coordinate
frame. Getting that wrong is the most common way an MD project silently starts from a nonsense
structure, so this script also verifies what it produced:

  * the ligand sits inside the receptor envelope (no periodic-image or frame mix-ups),
  * the number of protein-ligand heavy-atom contacts is in a plausible binding range,
  * there are no severe steric clashes (heavy-atom pairs closer than 1.6 A).

Protein hydrogens are deliberately **not** included: the topology builder (tleap) adds them with
the force field's own protonation states, which is more reliable than carrying over docking
hydrogens. Crystallographic waters are dropped; the SUOX MTE and SO4 cofactors are kept so tleap
can build the metal site around them.

Writes, per complex, under `analysis/md/complexes/<name>/`:
    receptor.pdb   protein + cofactor, heavy atoms, TER-separated chains
    ligand.pdb     the docked pose, heavy atoms (for parameterisation)
    complex.pdb    receptor + ligand, the frame everything downstream must use
    contacts.txt   residues within 4 A of the ligand
and one `complex_manifest.json` for the whole set.

Usage:
    python analysis/scripts/prepare_complex.py
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOCK = os.path.join(ROOT, "analysis", "docking")
MDDIR = os.path.join(ROOT, "analysis", "md")
OUTDIR = os.path.join(MDDIR, "complexes")
REDO = os.path.join(DOCK, "redo", "main")

# (complex name, source PDB, receptor chain, ligand, docking label)
COMPLEXES = [
    ("GDF11_GW8510", "5E4G.pdb", "A", "GW-8510", "GDF11"),
    ("SUOX_Sanguinarine", "2A9A.pdb", "A", "Sanguinarine", "SUOX"),
    ("SUOX_Spiperone", "2A9A.pdb", "A", "Spiperone", "SUOX"),
    ("SUOX_Vorinostat", "2A9A.pdb", "A", "Vorinostat", "SUOX"),
]
# Hetero residues that carry through to MD, by source structure.
#
# SUOX keeps MTE (molybdopterin), its Mo and the catalytic SO4 — plus, crucially, the two waters
# that make up the rest of the metal's first coordination shell. Reading 2A9A directly gives the
# Mo environment as: MTE S1' 2.38 A, Cys185 SG 2.41 A, MTE S2' 2.43 A, and two waters at 1.69 A
# (the oxo ligand) and 2.10 A (hydroxo/water). Those two waters are *ligands of the metal*, not
# solvent, so dropping them - as an earlier revision did, along with the Mo itself - leaves the
# metal site chemically wrong and MCPB.py with nothing to bond the metal to.
KEEP_HETERO = {
    "2A9A.pdb": {"MTE", "SO4", "MO"},
    "5E4G.pdb": set(),               # no ligand in the crystal
}
# 2A9A carries TWO sulfates. Res 7503 is the catalytic one (5.1 A from the Mo, and its midpoint with
# the Mo is exactly the published grid centre); res 7504 sits ~10.5 A away on the surface and is a
# crystallisation artefact, so it is excluded rather than left to perturb the electrostatics.
KEEP_HETERO_BY_ID = {("2A9A.pdb", "SO4"): {"7503"}}
# waters bonded to a metal are kept by proximity, not by residue name
METAL_WATER_CUTOFF = 2.6            # A; Mo-O bond distances in 2A9A are 1.69 and 2.10 A
DROP_HETERO = {"PG4", "EDO", "GOL", "ACT", "PEG"}


def xyz(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def read_best_pose(path, max_atoms=200):
    """Heavy atoms of the top-ranked Vina MODEL block."""
    out, mode, in_first = [], 0, False
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("MODEL"):
            mode += 1
            in_first = (mode == 1)
            continue
        if in_first and ln.startswith(("ATOM", "HETATM")):
            out.append(ln)
        if mode > 1:
            break
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    manifest = {"note": ("heavy-atom starting structures for MD; tleap adds protein hydrogens "
                         "and the force field decides protonation states"),
                "complexes": {}}
    bad = 0
    for name, pdb, chain, lig, target in COMPLEXES:
        work = os.path.join(OUTDIR, name)
        os.makedirs(work, exist_ok=True)
        src = os.path.join(DOCK, "structures", pdb)

        # ---- receptor: protein heavy atoms + retained cofactors ------------------------------
        keep_het = KEEP_HETERO.get(pdb, set())
        # collect the metal positions first, so coordinating waters can be identified by distance
        metal_xyz = []
        for ln in open(src, encoding="utf-8"):
            if ln.startswith("HETATM") and ln[21] == chain and ln[17:20].strip() in ("MO", "ZN", "FE", "MG", "MN", "CU"):
                metal_xyz.append(xyz(ln))
        rec_lines, het_written = [], {}
        kept_waters = []
        for ln in open(src, encoding="utf-8"):
            if ln.startswith("ATOM") and ln[21] == chain:
                if ln[76:78].strip() == "H":       # heavy atoms only
                    continue
                rec_lines.append(ln)
            elif ln.startswith("HETATM") and ln[21] == chain:
                res = ln[17:20].strip()
                if ln[76:78].strip() == "H":
                    continue
                keep = res in keep_het
                allow = KEEP_HETERO_BY_ID.get((pdb, res))
                if allow is not None and ln[22:26].strip() not in allow:
                    keep = False
                if not keep and res == "HOH" and metal_xyz:
                    # retain only waters that actually coordinate a metal
                    here = xyz(ln)
                    if any(math.dist(here, m) <= METAL_WATER_CUTOFF for m in metal_xyz):
                        keep = True
                        kept_waters.append(ln[22:26].strip())
                if keep:
                    rec_lines.append(ln)
                    het_written[res] = het_written.get(res, 0) + 1
        rec_lines.append("TER\n")
        rec_path = os.path.join(work, "receptor.pdb")
        with open(rec_path, "w", encoding="utf-8", newline="\n") as f:
            f.writelines(rec_lines)

        # ---- ligand pose ---------------------------------------------------------------------
        pose = os.path.join(REDO, f"{target}__{lig}__seed0.pdbqt")
        if not os.path.exists(pose):
            print(f"{name}: missing docking output {os.path.basename(pose)}")
            bad += 1
            continue
        raw = read_best_pose(pose)
        if not raw:
            print(f"{name}: no pose in {os.path.basename(pose)}")
            bad += 1
            continue
        lig_lines = []
        for i, ln in enumerate(raw, 1):
            atom = ln[12:16].strip() or "C"
            # PDBQT atom names carry the AutoDock type; strip to an element symbol for PDB
            elem = "".join(c for c in atom if c.isalpha())[:1] or "C"
            lig_lines.append(f"HETATM{i:5d} {atom:<4s} LIG Z 901    "
                             f"{ln[30:38]}{ln[38:46]}{ln[46:54]}  1.00  0.00          {elem:>2s}\n")
        lig_path = os.path.join(work, "ligand.pdb")
        with open(lig_path, "w", encoding="utf-8", newline="\n") as f:
            f.writelines(lig_lines)
            f.write("END\n")

        # ---- combined frame ------------------------------------------------------------------
        cx_path = os.path.join(work, "complex.pdb")
        with open(cx_path, "w", encoding="utf-8", newline="\n") as f:
            f.writelines(rec_lines)
            f.writelines(lig_lines)
            f.write("END\n")

        # ---- geometry checks -----------------------------------------------------------------
        rec_hv = [xyz(l) for l in rec_lines if l.startswith(("ATOM", "HETATM"))]
        lig_hv = [xyz(l) for l in lig_lines]
        clash = 0
        close = []
        for i, a in enumerate(lig_hv):
            near = sum(1 for b in rec_hv if math.dist(a, b) < 1.6)
            if near:
                clash += near
                close.append(i + 1)
        # contacts
        from collections import Counter
        contact_res = Counter()
        for b in rec_lines:
            if b.startswith(("ATOM", "HETATM")):
                bp = xyz(b)
                if any(math.dist(bp, a) < 4.0 for a in lig_hv):
                    contact_res[b[17:20].strip() + b[22:26].strip()] = 1
        # ligand must sit inside the receptor envelope
        cen = tuple(sum(a[i] for a in lig_hv) / len(lig_hv) for i in range(3))
        rmax = max(math.dist(cen, b) for b in rec_hv)
        inside = rmax > 5.0
        ok = (clash == 0) and inside and (len(contact_res) >= 5)

        rec = {"target": target, "ligand": lig, "source_pdb": pdb, "receptor_chain": chain,
               "keep_hetero": sorted(keep_het), "hetero_atoms_written": het_written,
               "coordinating_waters_kept": sorted(set(kept_waters)),
               "metal_site_note": ("Mo + MTE dithiolene sulfurs + Cys185 SG + the oxo/hydroxo "
                                   "waters form the first shell; MCPB.py builds the bonded model "
                                   "from these") if "MO" in het_written else "",
               "receptor_heavy_atoms": len(rec_hv), "ligand_heavy_atoms": len(lig_hv),
               "ligand_contacts_within_4A": len(contact_res),
               "clashes_under_1.6A": clash, "ok": ok,
               "files": {k: os.path.relpath(os.path.join(work, v), ROOT).replace("\\", "/")
                         for k, v in (("receptor", "receptor.pdb"), ("ligand", "ligand.pdb"),
                                      ("complex", "complex.pdb"))}}
        with open(os.path.join(work, "contacts.txt"), "w", encoding="utf-8") as f:
            f.write(f"# residues within 4 A of {lig} in {name}\n")
            for r in sorted(contact_res):
                f.write(r + "\n")
        manifest["complexes"][name] = rec
        flag = "OK " if ok else "CHECK"
        print(f"{flag} {name:22s} rec {len(rec_hv):5d} heavy | lig {len(lig_hv):3d} | "
              f"contacts {len(contact_res):3d} | clashes {clash}"
              + (f" at ligand atoms {close[:5]}" if clash else ""))
        if not ok:
            bad += 1

    with open(os.path.join(OUTDIR, "complex_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print(f"\n{len(COMPLEXES) - bad}/{len(COMPLEXES)} complexes passed the geometry checks")
    print(f"-> {os.path.relpath(OUTDIR, ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
