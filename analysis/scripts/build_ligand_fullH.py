#!/usr/bin/env python3
"""build_ligand_fullH.py — correct, charge-aware hydrogen completion for parameterisation.

Why this exists
---------------
Two successive bugs came from deriving ligand connectivity out of the docking pose file:

  1. `ligand.pdb` came from the AutoDock PDBQT, which follows the united-atom convention and writes
     only **polar** hydrogens. antechamber then inferred the wrong valence, produced an **odd
     electron count**, and sqm refused outright. Fixed by adding hydrogens.
  2. Adding them with `obabel -h` is charge-blind: for a permanent cation such as sanguinarine
     (quaternary N+) it adds a hydrogen to a nitrogen that has no hydrogen, turning the intended
     +1 species into something else and producing an odd electron count again (173).

So neither the pose file nor a generic hydrogen-adder can be trusted for the protonation state.
The authoritative source is `smiles_protonated` in the docking ligand manifest, which records the
explicit, auditable protonation policy applied to each molecule (carboxylic acids -> -1, aliphatic
amines -> +1, phenols/sulfonamides/thiols neutral, permanent charges kept).

What this does
--------------
1. Reads the protonated SMILES from `analysis/docking/ligands/ligand_manifest.json`.
2. Builds that molecule in RDKit, so bonds **and formal charges** are exactly as documented.
3. Overwrites the heavy-atom coordinates with the **docked pose** coordinates, then calls
   `AddHs(addCoords=True)` so the hydrogens are placed on the real pose geometry.
   The heavy atoms keep the pose bit-for-bit, which is the point: MD must start from the docked
   binding mode, not from a fresh conformer.
4. Verifies the heavy-atom coordinates were not disturbed, and reports the electron count that will
   be handed to sqm, so an odd count can never slip through again.

Writes `ligand_fullH.pdb` and `ligand_fullH.sdf` in each complex directory. Parameterise
from the **SDF**: a PDB carries no bond orders, and antechamber's bond perception then fails on the
sulfonamide sulfur of GW-8510 ("the best APS is not zero ... bonds involving this atom are frozen").

Usage:
    python analysis/scripts/build_ligand_fullH.py [complex ...]
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MD = os.path.join(ROOT, "analysis", "md")
MANIFEST = os.path.join(ROOT, "analysis", "docking", "ligands", "ligand_manifest.json")


def read_pose_heavy(path):
    """[(element, x, y, z)] for heavy atoms, in file order."""
    out = []
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith(("ATOM", "HETATM")):
            name = ln[12:16].strip()
            elem = (ln[76:78].strip() or name[:1]).upper()
            elem = "".join(c for c in elem if c.isalpha())[:1] or "C"
            if elem == "H":
                continue
            out.append((elem, float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return out


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D
    RDLogger.DisableLog("rdApp.*")

    complexes = sys.argv[1:] or ["GDF11_GW8510", "SUOX_Sanguinarine",
                                 "SUOX_Spiperone", "SUOX_Vorinostat"]
    man = json.load(open(MANIFEST, encoding="utf-8"))
    by_cx = json.load(open(os.path.join(MD, "complexes", "complex_manifest.json"),
                           encoding="utf-8"))["complexes"]

    bad = 0
    for cx in complexes:
        d = os.path.join(MD, "complexes", cx)
        if cx not in by_cx:
            print(f"{cx}: not in complex_manifest.json"); bad += 1; continue
        lig = by_cx[cx]["ligand"]
        rec = None
        for k, v in man["ligands"].items():
            if v["name"] == lig and v["target"] == by_cx[cx]["target"]:
                rec = v
                break
        if rec is None:
            print(f"{cx}: no manifest record for {lig}"); bad += 1; continue

        smi = rec.get("smiles_protonated") or rec["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"{cx}: RDKit could not parse {smi!r}"); bad += 1; continue

        pose = None
        for cand in ("ligand.pdbqt", "ligand.pdb"):
            p = os.path.join(d, cand)
            if os.path.exists(p):
                pose = read_pose_heavy(p)
                break
        if not pose:
            print(f"{cx}: no pose file"); bad += 1; continue

        n_heavy = mol.GetNumHeavyAtoms()
        if len(pose) != n_heavy:
            print(f"{cx}: pose has {len(pose)} heavy atoms but {smi} has {n_heavy} - "
                  f"refusing to guess the correspondence")
            bad += 1
            continue

        # put the pose geometry on the RDKit molecule (heavy atoms stay put), then add H onto it
        conf = Chem.Conformer(n_heavy)
        for i, (_e, x, y, z) in enumerate(pose):
            conf.SetAtomPosition(i, Point3D(x, y, z))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        mol_h = Chem.AddHs(mol, addCoords=True)

        # electron bookkeeping, so an odd count cannot slip through
        z = sum(a.GetAtomicNum() for a in mol_h.GetAtoms())
        q = Chem.GetFormalCharge(mol_h)
        n_elec = z - q
        n_h = sum(1 for a in mol_h.GetAtoms() if a.GetSymbol() == "H")

        out = os.path.join(d, "ligand_fullH.pdb")
        Chem.MolToPDBFile(mol_h, out)
        # SDF carries explicit bond orders, so antechamber never has to infer them
        sdf = os.path.join(d, "ligand_fullH.sdf")
        w = Chem.SDWriter(sdf)
        w.write(mol_h)
        w.close()

        # verify the heavy atoms did not move
        after = read_pose_heavy(out)
        moved = sum(1 for (_e, *c1), (_e2, *c2) in zip(pose, after)
                    if sum(abs(a - b) for a, b in zip(c1, c2)) > 0.01)
        parity = "even-OK" if n_elec % 2 == 0 else "ODD-BAD"
        flag = "OK" if (moved == 0 and n_elec % 2 == 0) else "CHECK"
        if flag != "OK":
            bad += 1
        print(f"{flag:5s} {cx:20s} {lig:14s} q={q:+d} heavy={n_heavy} H={n_h:2d} "
              f"electrons={n_elec} ({parity})  moved={moved}")

    print(f"\n{len(complexes) - bad}/{len(complexes)} ligands built with a valid closed-shell parity")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
