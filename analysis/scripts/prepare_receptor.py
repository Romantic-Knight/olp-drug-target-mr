"""prepare_receptor.py — build AutoDock Vina receptors for the OLP docking redo.

Both receptors are built from the raw crystal files with the same recipe, so the two targets are
comparable and the provenance is recorded in one place.

GDF11 (PDB 5E4G)
    Human GDF11 mature domain, residues 299-407, one chain in the asymmetric unit. Four intrachain
    disulfides (Cys6-Cys16, Cys15-Cys74, Cys43-Cys106, Cys47-Cys108) are kept; SSBOND 5
    (Cys73-Cys73, symmetry 7555) is the INTERCHAIN disulfide that builds the biological dimer, so
    Cys73 is left free on the monomer. Two PG4 (polyethylene glycol) copies are crystallisation
    artefacts and are dropped. The previous docking receptor also had NO hydrogens (901 atoms =
    the hydrogen-free chain); the fix here is the same as for SUOX.

SUOX (PDB 2A9A)
    Chicken sulfite oxidase, chain A. Its catalytic centre is the molybdopterin cofactor MTE plus
    a molybdenum ion. AutoDock Vina cannot parse molybdenum at all - appending the crystal Mo gives
    "Atom type Mo is not a valid AutoDock type", and the uppercase spelling is rejected identically,
    while Fe/Zn/Mg parse fine. MTE itself contains only C/N/O/P/S, all standard AutoDock types, so
    the single Mo atom is dropped and MTE plus the substrate analogue SO4 are KEPT: this restores
    the real shape of the pocket instead of emptying the cavity, which is what the published
    receptor did. The Mo is restored at the MD stage (MCPB.py), where it becomes a test of whether
    a pose docked into the Mo-free pocket survives once the metal is present.

Charges: Vina's default scoring function reads no partial charges, so a rigid receptor needs only
coordinates and AutoDock types. Explicit charges would be required only for the AD4 scoring
function (which Vina supports as receptor-less flexible docking only).

Hydrogens: X-ray structures carry none, but Vina's hydrogen-bond term needs donor/acceptor typing on
BOTH sides, so both receptors are protonated with OpenBabel at the given pH.

Output records `receptors/receptor_manifest.json`: source file, pH, atoms, hydrogens, kept/excluded
hetero groups, disulfides and atom-type inventory - the audit trail the previous docking lacked.

Usage:
    python analysis/scripts/prepare_receptor.py            # both targets
    python analysis/scripts/prepare_receptor.py gdf11 suox
"""
import argparse
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = os.path.join(ROOT, "analysis", "docking")
OBABEL = os.path.join(BASE, ".pixi", "envs", "default", "Library", "bin", "obabel.exe")
STRUCT = os.path.join(BASE, "structures")
RECDIR = os.path.join(BASE, "receptors")

AD_TYPES = {"A", "C", "N", "NA", "OA", "SA", "S", "P", "HD", "H", "F", "Cl", "Br", "I",
            "Fe", "Zn", "Mg", "Mn", "Ca"}

TARGETS = {
    "gdf11": {
        "source": "5E4G.pdb",
        "out": "5E4G_chainA_H.pdbqt",
        "pdb_out": "5E4G_chainA_H.pdb",
        "chain": "A",
        "keep_hetero": set(),                    # no functional cofactor
        "drop_hetero": {"PG4", "HOH"},           # PEG = crystallisation artefact
        "expect_excluded": {},                   # nothing essential is removed
    },
    "suox": {
        "source": "2A9A.pdb",
        "out": "2A9A_chainA_MTE_SO4.pdbqt",
        "pdb_out": "2A9A_chainA_MTE_SO4.pdb",
        "chain": "A",
        "keep_hetero": {"MTE", "SO4"},           # molybdopterin + substrate analogue
        "drop_hetero": {"MO", "HOH"},            # Mo is not a Vina atom type (see module docstring)
        "expect_excluded": {"MO": 1},
    },
}


def read_disulfides(pdb):
    out = []
    for ln in open(pdb, encoding="utf-8"):
        if ln.startswith("SSBOND"):
            parts = ln.split()
            out.append({
                "cys1": f"{parts[3]}{parts[4]}", "cys2": f"{parts[6]}{parts[7]}",
                "sym1": parts[8], "sym2": parts[9],
                "interchain": parts[8] != parts[9],   # symmetry mismatch = interchain bond
            })
    return out


def extract(cfg, pdb, out_pdb):
    kept, het, dropped, excluded = [], {}, {}, {}
    for ln in open(pdb, encoding="utf-8"):
        if ln.startswith("ATOM"):
            if ln[21] == cfg["chain"]:
                kept.append(ln)
        elif ln.startswith("HETATM") and ln[21] == cfg["chain"]:
            res = ln[17:20].strip()
            if res in cfg["keep_hetero"]:
                kept.append(ln)
                het[res] = het.get(res, 0) + 1
            elif res in cfg["drop_hetero"]:
                excluded[res] = excluded.get(res, 0) + 1
            else:
                dropped[res] = dropped.get(res, 0) + 1
    with open(out_pdb, "w", encoding="utf-8") as f:
        f.writelines(kept)
        f.write("END\n")
    return kept, het, dropped, excluded


def protonate_and_rigidify(pdb, raw, out, ph):
    r = subprocess.run([OBABEL, "-ipdb", pdb, "-opdbqt", "-O", raw, "-h", "-p", str(ph)],
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(raw):
        raise SystemExit(f"obabel failed ({r.returncode}): {(r.stderr or '').strip()[:200]}")
    types, unknown, n = {}, set(), 0
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(raw, encoding="utf-8") as fi, open(out, "w", encoding="utf-8") as fo:
        for ln in fi:
            # OpenBabel writes ligand-style BRANCH/ROOT/TORSDOF tags for HETATM groups and Vina's
            # rigid-receptor parser rejects the file as a whole ("Unknown or inappropriate tag").
            # A rigid receptor needs coordinates and AutoDock types only, so this reduction is
            # exact, not a workaround.
            if not ln.startswith(("ATOM", "HETATM")) or len(ln.rstrip("\n")) < 79:
                continue
            at = ln[77:].strip()
            types[at] = types.get(at, 0) + 1
            if at not in AD_TYPES:
                unknown.add(at)
            fo.write(ln if ln.endswith("\n") else ln + "\n")
            n += 1
        fo.write("ENDMDL\n")
    os.remove(raw)
    return n, types, unknown


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("targets", nargs="*", default=None,
                    help="subset of " + ", ".join(TARGETS) + " (default: all)")
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()
    names = a.targets or list(TARGETS)
    for n in names:
        if n not in TARGETS:
            raise SystemExit(f"unknown target {n!r}; choose from {list(TARGETS)}")

    manifest_path = os.path.join(RECDIR, "receptor_manifest.json")
    manifest = {}
    if os.path.exists(manifest_path):
        manifest = json.load(open(manifest_path, encoding="utf-8"))

    for name in names:
        cfg = TARGETS[name]
        src = os.path.join(STRUCT, cfg["source"])
        pdb_out = os.path.join(STRUCT, cfg["pdb_out"])
        rec_out = os.path.join(RECDIR, cfg["out"])

        print(f"=== {name.upper()}  ({cfg['source']})")
        sss = read_disulfides(src)
        intra = [s for s in sss if not s["interchain"]]
        inter = [s for s in sss if s["interchain"]]
        print(f"  disulfides: {len(intra)} intrachain, {len(inter)} interchain"
              + (f"  <- {inter[0]['cys1']}-{inter[0]['cys2']} (dimer interface, left free)"
                 if inter else ""))

        kept, het, dropped, excluded = extract(cfg, src, pdb_out)
        print(f"  chain {cfg['chain']}: {len(kept)} atom lines | kept hetero {het or '{}'} | "
              f"excluded {excluded or '{}'} | other hetero dropped {sorted(dropped)}")
        for res, want in cfg["expect_excluded"].items():
            got = excluded.get(res, 0)
            if got != want:
                raise SystemExit(f"expected {want} {res} to be excluded, found {got}")

        raw = os.path.join(STRUCT, "_obabel_raw.pdbqt")
        n, types, unknown = protonate_and_rigidify(pdb_out, raw, rec_out, a.ph)
        n_h = sum(v for k, v in types.items() if k in ("HD", "H"))
        if n_h == 0:
            raise SystemExit("no hydrogens were added - protonation failed")
        if unknown:
            raise SystemExit(f"non-standard AutoDock types present: {sorted(unknown)}")
        print(f"  pH {a.ph}: {n} atoms, {n_h} hydrogens, types {dict(sorted(types.items()))}")
        print(f"  -> {os.path.relpath(rec_out, ROOT)}")

        manifest[name] = {
            "source_pdb": cfg["source"], "chain": cfg["chain"], "ph": a.ph,
            "receptor_pdbqt": os.path.relpath(rec_out, ROOT).replace("\\", "/"),
            "prepared_pdb": os.path.relpath(pdb_out, ROOT).replace("\\", "/"),
            "atoms": n, "hydrogens": n_h, "atom_types": dict(sorted(types.items())),
            "hetero_kept": het, "hetero_excluded": excluded,
            "hetero_dropped": sorted(dropped),
            "disulfides_intrachain": intra, "disulfides_interchain": inter,
            "excluded_note": ("Mo omitted: AutoDock Vina cannot parse molybdenum; restored at the "
                              "MD stage via MCPB.py") if excluded.get("MO") else "",
        }
        print()

    os.makedirs(RECDIR, exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"manifest -> {os.path.relpath(manifest_path, ROOT)}  ({len(manifest)} target(s))")


if __name__ == "__main__":
    main()
