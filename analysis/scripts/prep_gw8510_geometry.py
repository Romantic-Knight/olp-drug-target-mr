#!/usr/bin/env python3
"""prep_gw8510_geometry.py — relax the ligand geometry before the AM1 step.

Why
---
GW-8510's AM1 SCF fails to converge in sqm's hard limit of 1000 iterations, stalling at a density
change of ~1.5e-3 regardless of the requested tolerance. One common cause is a strained input
geometry: the SDF came from RDKit's distance-geometry embedding, which is bonded correctly but not
energetically relaxed, and a poor starting geometry makes SCF convergence much harder.

This writes a MMFF94-optimised conformer, keeping every atom and the bond orders intact, so the QM
step starts from a reasonable structure. Hydrogens are explicit and are re-optimised with the heavy
atoms.

The output replaces ligand_qm.sdf; nothing that the docking used is touched.
"""
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

path = sys.argv[1] if len(sys.argv) > 1 else "ligand_fullH.sdf"
out = sys.argv[2] if len(sys.argv) > 2 else "ligand_qm.sdf"

mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=True)
if mol is None:
    raise SystemExit(f"could not read {path}")

n_atoms = mol.GetNumAtoms()
n_bonds = mol.GetNumBonds()
charges = sum(a.GetFormalCharge() for a in mol.GetAtoms())
print(f"read {path}: {n_atoms} atoms, {n_bonds} bonds, formal charge {charges:+d}")

if AllChem.MMFFHasAllMoleculeParams(mol):
    ff = "MMFF94"
    props = AllChem.MMFFGetMoleculeProperties(mol)
    field = AllChem.MMFFGetMoleculeForceField(mol, props)
else:
    ff = "UFF"
    field = AllChem.UFFGetMoleculeForceField(mol)

e0 = field.CalcEnergy()
field.Minimize(maxIts=2000)
e1 = field.CalcEnergy()
print(f"  {ff}: E {e0:.2f} -> {e1:.2f} kcal/mol")

# confirm the chemistry survived the relaxation
mol2 = Chem.Mol(mol)
Chem.SanitizeMol(mol2)
if mol2.GetNumAtoms() != n_atoms or mol2.GetNumBonds() != n_bonds:
    raise SystemExit("atom/bond count changed during optimisation - aborting")
print(f"  verified: still {mol2.GetNumAtoms()} atoms, {mol2.GetNumBonds()} bonds")

Chem.MolToMolFile(mol, out)
print(f"wrote {out} ({os.path.getsize(out)} bytes)")
