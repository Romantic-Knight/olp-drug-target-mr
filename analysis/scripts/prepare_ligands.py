"""prepare_ligands.py — regenerate every docked ligand reproducibly at pH 7.4.

Why the protonation step is explicit rather than delegated
----------------------------------------------------------
`dimorphite_dl` was tried first and is reliable for single ionisable groups: at pH 7.4 it returns
phenol neutral, acetic/benzoic acid -1, methylamine +1, pyridine neutral (9/9 single-group probes
correct in a wider test). But it treats each acidic group independently, ignoring the electrostatic
penalty of neighbouring charges, so it **over-ionises polyphenols**: catechol comes back as -2 and
Quercetin as C15H5O7^-5, when the real sequential pKa values (catechol 9.4 / 12.8; flavonols ~7-11)
leave both close to neutral at pH 7.4. Three ligands here are polyphenols, so that error would have
propagated straight into the docking.

Policy used instead (uniform, so it is not a per-ligand judgement call):

    carboxylic acid            -> deprotonated (-1)      pKa ~4-5, fully ionised at 7.4
    aliphatic 1/2/3-amine      -> protonated   (+1)      pKa ~9-11, fully ionised at 7.4
    phenol, polyphenol         -> neutral                first pKa ~7-10, i.e. largely protonated;
                                                         keeping the neutral form also retains every
                                                         OH as both donor and acceptor
    sulfonamide, hydroxamic
    acid, thiol, amide         -> neutral                pKa > 9
    permanent charges          -> kept as drawn          e.g. the quaternary N of sanguinarine
    counter-ions               -> removed                e.g. the Cl- of hemin

Every ligand records which rules fired, so the state is auditable. The three polyphenols
(Quercetin, Apigenin) and Niclosamide are flagged `protonation_approximate`: their true species at
pH 7.4 is partly anionic, and the neutral form is the conservative choice for docking (it never
removes a hydrogen-bond donor).

Geometry and PDBQT are unchanged from the earlier revision and are the parts that work:
RDKit ETKDG + MMFF94 for coordinates (OpenBabel's `--gen3d` cannot run here - its fragment library
is missing from this environment), then OpenBabel SDF -> PDBQT for AutoDock typing and the torsion
tree.

Usage:
    python analysis/scripts/prepare_ligands.py [--ph 7.4]
"""
import argparse
import json
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = os.path.join(ROOT, "analysis", "docking")
OBABEL = os.path.join(BASE, ".pixi", "envs", "default", "Library", "bin", "obabel.exe")
LIGDIR = os.path.join(BASE, "ligands")
OUTDIR = os.path.join(LIGDIR, "prepared")
WORK = os.path.join(LIGDIR, "_work")
SOURCES = {"GDF11": "GDF11.smi", "SUOX": "SUOX.smi"}
SEED = 42

# SMARTS -> (action, description). Applied in order; each match is modified once.
RULES = [
    ("[CX3](=O)[OX2H1]", "deprotonate", "carboxylic acid -> carboxylate"),
    # H0 included so tertiary aliphatic amines (e.g. the piperidine N of spiperone) are caught;
    # `!$(N-[a])` keeps aniline-type amines neutral and `[CX4]` requires an sp3 neighbour.
    ("[NX3;H0,H1,H2;!$(NC=O);!$(NS=O);!$(N=O);!$([N+]);!$([nX3]);!$([NX3]=*);!$(N-[a])][CX4]",
     "protonate", "aliphatic amine -> ammonium"),
]
# Ligands whose true pH 7.4 species is partly anionic at a phenol the uniform policy leaves neutral.
APPROXIMATE = {"Quercetin", "Apigenin", "Niclosamide"}


def read_smi(path):
    out = []
    for ln in open(path, encoding="utf-8"):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t") if "\t" in ln else ln.split(None, 1)
        if len(parts) == 2:
            out.append({"smiles": parts[0].strip(), "name": parts[1].strip()})
    return out


def apply_rules(mol, Chem):
    """Apply the protonation policy in place. Returns the list of descriptions that fired."""
    from rdkit.Chem import Descriptors
    fired = []
    for smarts, action, desc in RULES:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            raise RuntimeError(f"bad SMARTS: {smarts}")
        n = 0
        while True:
            match = mol.GetSubstructMatch(patt)
            if not match:
                break
            rw = Chem.RWMol(mol)
            if action == "deprotonate":
                # the O with the hydrogen is the last matched atom for this SMARTS layout
                idx = [i for i in match if mol.GetAtomWithIdx(i).GetSymbol() == "O"
                       and mol.GetAtomWithIdx(i).GetTotalNumHs() > 0][-1]
                atom = rw.GetAtomWithIdx(idx)
                atom.SetFormalCharge(-1)
                atom.SetNumExplicitHs(0)
            else:
                idx = [i for i in match if mol.GetAtomWithIdx(i).GetSymbol() == "N"][0]
                atom = rw.GetAtomWithIdx(idx)
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
            mol = rw.GetMol()
            Chem.SanitizeMol(mol)
            n += 1
            if n > 12:                      # safety valve
                break
        if n:
            fired.append(f"{desc} x{n}")
    return mol, fired


def embed_metalloporphyrin(mol, Chem, AllChem):
    """Embed a metalloporphyrin whose metal was written as a separate SMILES component.

    RDKit cannot assign geometry to an unbonded metal: hemin's SMILES lists [Fe+3] as its own
    component, so ETKDG leaves it at the origin and the ligand is unusable. Here the organic
    macrocycle is embedded and optimised with UFF on its own, then the metal is placed at the
    centroid of the four pyrrole nitrogens - the position it occupies in the real cofactor.

    Returns (mol_with_metal_placed, note).
    """
    from rdkit.Geometry import Point3D
    metal_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Fe"]
    if len(metal_atoms) != 1:
        raise RuntimeError(f"expected exactly one Fe, found {len(metal_atoms)}")
    fe = metal_atoms[0]
    organic = Chem.RWMol(mol)
    organic.RemoveAtom(fe)
    organic = organic.GetMol()
    Chem.SanitizeMol(organic)
    if AllChem.EmbedMolecule(organic, randomSeed=SEED) != 0:
        raise RuntimeError("ETKDG embedding of the porphyrin failed")
    try:
        AllChem.UFFOptimizeMolecule(organic)
    except Exception:
        pass                                    # a metal-free macrocycle still embeds acceptably

    n_idx = [a.GetIdx() for a in organic.GetAtoms()
             if a.GetSymbol() == "N" and a.GetTotalNumHs() == 0]
    if len(n_idx) < 4:
        raise RuntimeError(f"found only {len(n_idx)} coordinating nitrogens")
    conf = organic.GetConformer()
    pts = [conf.GetAtomPosition(i) for i in n_idx[:4]]
    c = Point3D(sum(p.x for p in pts) / 4.0, sum(p.y for p in pts) / 4.0,
                sum(p.z for p in pts) / 4.0)

    n_org = organic.GetNumAtoms()
    combined = Chem.RWMol(organic)
    combined.AddAtom(mol.GetAtomWithIdx(fe))
    combined = combined.GetMol()
    Chem.SanitizeMol(combined)
    conf = Chem.Conformer(n_org + 1)
    for i in range(n_org):
        conf.SetAtomPosition(i, organic.GetConformer().GetAtomPosition(i))
    conf.SetAtomPosition(n_org, c)
    combined.RemoveAllConformers()
    combined.AddConformer(conf, assignId=True)
    note = ("porphyrin embedded separately; Fe placed at the centroid of the four pyrrole N")
    return combined, note


METALS = ("Fe", "Zn", "Mg", "Mn", "Cu", "Mo", "Co", "Ni", "Ca")


def strip_counterions(mol, Chem):
    """Desalt: keep the largest organic fragment plus any metal-bearing fragment.

    Two distinct cases matter here.

      * **Single-atom counter-ions** - the Cl- of hemin, alkali halides as listed in the NCI
        database. Dropped.
      * **Multi-fragment database entries** - 1:1 salts or co-crystals, e.g. four NCI decoys are
        recorded as "benzoic acid . hydroxyquinoline". These make OpenBabel emit two ROOT blocks
        and AutoDock Vina then refuses the ligand outright ("Unknown or inappropriate tag found in
        flex residue"), so the smaller fragment is dropped too.

    A fragment containing a metal is **always** kept: hemin's Fe is written as its own SMILES
    component, and it has to survive so `embed_metalloporphyrin` can place it.

    Returns (mol, dropped_smiles).
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) <= 1:
        return mol, []

    scored = []
    for f in frags:
        syms = {a.GetSymbol() for a in f.GetAtoms()}
        scored.append({"mol": f, "metal": bool(set(syms) & set(METALS)),
                       "heavy": f.GetNumHeavyAtoms(), "smiles": Chem.MolToSmiles(f)})

    # the largest fragment that is more than a bare ion, so a salt like "Na.Cl" does not end up
    # represented by a single atom
    organic = [t for t in scored if t["heavy"] > 1]
    biggest = max(organic or scored, key=lambda t: t["heavy"])

    keep, dropped = [], []
    for t in scored:
        if t is biggest or t["metal"]:
            keep.append(t["mol"])
        else:
            dropped.append(t["smiles"])
    if len(keep) == len(frags) or not keep:
        return mol, []

    out = keep[0]
    for f in keep[1:]:
        out = Chem.CombineMols(out, f)
    out = Chem.RWMol(out).GetMol()
    Chem.SanitizeMol(out)
    return out, dropped


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    ap = argparse.ArgumentParser()
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)

    manifest = {"ph": a.ph, "seed": SEED,
                "hydrogen_policy": "polar hydrogens only (PDBQT standard)",
                "protonation_policy": [d for _, _, d in RULES]
                                      + ["phenol/sulfonamide/thiol/amide left neutral",
                                         "single-atom counter-ions removed"],
                "ligands": {}}
    ok = 0
    for target, smi_file in SOURCES.items():
        ligs = read_smi(os.path.join(LIGDIR, smi_file))
        print(f"=== {target}: {len(ligs)} ligands from {smi_file}")
        for lg in ligs:
            name = lg["name"].replace(" ", "_")
            rec = {"target": target, "name": name, "smiles": lg["smiles"]}
            try:
                mol = Chem.MolFromSmiles(lg["smiles"])
                if mol is None:
                    raise RuntimeError("RDKit could not parse the SMILES")
                mol, dropped = strip_counterions(mol, Chem)
                if dropped:
                    rec["counterions_removed"] = dropped
                mol, fired = apply_rules(mol, Chem)
                rec["protonation_rules_fired"] = fired
                # the authoritative record of the species: parameterisation must build from this,
                # not from the pose file, whose hydrogens are incomplete and charge-blind
                rec["smiles_protonated"] = Chem.MolToSmiles(mol)
                rec["net_charge"] = Chem.GetFormalCharge(mol)
                if name in APPROXIMATE:
                    rec["protonation_approximate"] = (
                        "true pH 7.4 species is partly anionic at a phenol left neutral here")

                mol = Chem.AddHs(mol)
                metal = None
                if any(a.GetSymbol() == "Fe" for a in mol.GetAtoms()):
                    mol, metal = embed_metalloporphyrin(mol, Chem, AllChem)
                else:
                    if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
                        raise RuntimeError("ETKDG embedding failed")
                    AllChem.MMFFOptimizeMolecule(mol)
                if metal:
                    rec["geometry_note"] = metal
                rec.update({"formula": rdMolDescriptors.CalcMolFormula(mol),
                            "heavy_atoms": mol.GetNumHeavyAtoms(),
                            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                            "mw": round(Descriptors.MolWt(mol), 2)})

                sdf = os.path.join(WORK, name + ".sdf")
                w = Chem.SDWriter(sdf)
                w.write(mol)
                w.close()

                out = os.path.join(OUTDIR, name + ".pdbqt")
                r = subprocess.run([OBABEL, "-isdf", sdf, "-opdbqt", "-O", out],
                                   capture_output=True, text=True)
                if r.returncode != 0 or not os.path.exists(out):
                    raise RuntimeError(f"pdbqt conversion failed: {(r.stderr or '').strip()[:120]}")

                metalfree = None
                if any(x.GetSymbol() == "Fe" for x in mol.GetAtoms()):
                    # AutoDock Vina has no ligand-side Fe atom type: it rejects such a ligand with
                    # "Unknown or inappropriate tag found in flex residue". Emit the metal-free
                    # porphyrin as an explicit, documented surrogate so the ligand list stays
                    # complete instead of silently dropping one entry.
                    rw = Chem.RWMol(mol)
                    for idx in sorted([x.GetIdx() for x in mol.GetAtoms()
                                       if x.GetSymbol() == "Fe"], reverse=True):
                        rw.RemoveAtom(idx)
                    mf = rw.GetMol()
                    Chem.SanitizeMol(mf)
                    mf = Chem.AddHs(mf, addCoords=True)
                    mf_sdf = os.path.join(WORK, name + "_metalfree.sdf")
                    w2 = Chem.SDWriter(mf_sdf)
                    w2.write(mf)
                    w2.close()
                    mf_out = os.path.join(OUTDIR, name + "_metalfree.pdbqt")
                    r2 = subprocess.run([OBABEL, "-isdf", mf_sdf, "-opdbqt", "-O", mf_out],
                                        capture_output=True, text=True)
                    if r2.returncode == 0 and os.path.exists(mf_out):
                        metalfree = os.path.relpath(mf_out, ROOT).replace("\\", "/")
                        rec["metal_free_surrogate"] = metalfree
                        rec["metal_free_note"] = ("Vina cannot dock a ligand containing Fe; the "
                                                  "metal-free porphyrin is docked instead and "
                                                  "reported as a surrogate")

                n_atoms = n_h = n_zero = torsdof = 0
                for ln in open(out, encoding="utf-8"):
                    if ln.startswith(("ATOM", "HETATM")):
                        n_atoms += 1
                        if ln[76:79].strip() in ("H", "HD"):
                            n_h += 1
                        if ln[30:54].strip() == "0.000   0.000   0.000":
                            n_zero += 1
                    elif ln.startswith("TORSDOF"):
                        torsdof = int(ln.split()[1])
                if n_zero:
                    raise RuntimeError(f"{n_zero} atoms at the origin - coordinates lost")
                rec.update({"pdbqt_atoms": n_atoms, "polar_hydrogens": n_h,
                            "torsdof": torsdof, "ok": True,
                            "pdbqt": os.path.relpath(out, ROOT).replace("\\", "/")})
                approx = " (approx)" if name in APPROXIMATE else ""
                print(f"  {name:28s} {rec['formula']:15s} q={rec['net_charge']:+d} "
                      f"heavy={rec['heavy_atoms']:3d} rot={rec['rotatable_bonds']:2d} "
                      f"polarH={n_h:2d} torsdof={torsdof:2d}{approx}")
                ok += 1
            except Exception as e:
                rec.update({"ok": False, "error": f"{type(e).__name__}: {e}"})
                print(f"  {name:28s} FAILED: {rec['error'][:88]}")
            # keyed per (target, ligand): two ligands are docked into both receptors, and a
            # name-only key silently overwrote the first target's record.
            manifest["ligands"][f"{target}:{name}"] = rec

    path = os.path.join(LIGDIR, "ligand_manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")
    bad = [k for k, v in manifest["ligands"].items() if not v.get("ok")]
    n_uniq = len({v["name"] for v in manifest["ligands"].values()})
    print(f"\n{ok}/{len(manifest['ligands'])} target-ligand pairs prepared "
          f"({n_uniq} unique ligands); failures: {bad or 'none'}")
    print(f"manifest -> {os.path.relpath(path, ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
