"""prepare_decoys.py — prepare the decoy set to PDBQT with the SAME pipeline as the ligands.

The decoys must be prepared exactly like the candidates, otherwise a score difference could come
from the preparation rather than from binding. The protonation rules, geometry generation and
PDBQT conversion are therefore imported from `prepare_ligands.py` rather than reimplemented.

Writes analysis/docking/decoys/prepared/*.pdbqt and decoys/decoy_prepared_manifest.json.

Usage:
    python analysis/scripts/prepare_decoys.py
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import prepare_ligands as pl          # noqa: E402  (reuse the validated rules)


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    src = os.path.join(pl.LIGDIR, "..", "decoys", "decoys.smi")
    src = os.path.abspath(src)
    outdir = os.path.join(os.path.dirname(src), "prepared")
    os.makedirs(outdir, exist_ok=True)
    work = os.path.join(os.path.dirname(src), "_work")
    os.makedirs(work, exist_ok=True)

    if not os.path.exists(src):
        print(f"missing {src} - run build_decoys.py first")
        return 1

    entries = pl.read_smi(src)
    print(f"{len(entries)} decoys to prepare (same rules as the candidates)")

    manifest, ok, zero, noH = {}, 0, 0, 0
    for e in entries:
        name, smi = e["name"], e["smiles"]
        rec = {"smiles": smi}
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise RuntimeError("unparseable")
            mol, _dropped = pl.strip_counterions(mol, Chem)
            mol, fired = pl.apply_rules(mol, Chem)
            rec["net_charge"] = Chem.GetFormalCharge(mol)
            if fired:
                rec["protonation_rules_fired"] = fired
            mol = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol, randomSeed=pl.SEED) != 0:
                raise RuntimeError("ETKDG embedding failed")
            AllChem.MMFFOptimizeMolecule(mol)
            rec["mw"] = round(Descriptors.MolWt(mol), 2)
            rec["heavy_atoms"] = mol.GetNumHeavyAtoms()
            rec["rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)

            sdf = os.path.join(work, name + ".sdf")
            w = Chem.SDWriter(sdf)
            w.write(mol)
            w.close()

            out = os.path.join(outdir, name + ".pdbqt")
            r = subprocess.run([pl.OBABEL, "-isdf", sdf, "-opdbqt", "-O", out],
                               capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(out):
                raise RuntimeError(f"pdbqt: {(r.stderr or '').strip()[:80]}")

            n = nz = nh = tors = 0
            for ln in open(out, encoding="utf-8"):
                if ln.startswith(("ATOM", "HETATM")):
                    n += 1
                    if ln[76:79].strip() in ("H", "HD"):
                        nh += 1
                    if ln[30:54].strip() == "0.000   0.000   0.000":
                        nz += 1
                elif ln.startswith("TORSDOF"):
                    tors = int(ln.split()[1])
            zero += nz
            if nh == 0:
                noH += 1
            rec.update({"pdbqt_atoms": n, "polar_hydrogens": nh, "torsdof": tors,
                        "ok": nz == 0})
            if nz:
                rec["error"] = f"{nz} atoms at the origin"
            else:
                ok += 1
        except Exception as ex:
            rec.update({"ok": False, "error": f"{type(ex).__name__}: {ex}"})
        manifest[name] = rec

    bad = [k for k, v in manifest.items() if not v.get("ok")]
    with open(os.path.join(os.path.dirname(src), "decoy_prepared_manifest.json"),
              "w", encoding="utf-8") as f:
        json.dump({"n": len(manifest), "ok": ok, "failures": bad,
                   "source": pl.__file__.replace("\\", "/").rsplit("/", 1)[-1],
                   "note": "prepared with the identical pipeline as the candidate ligands",
                   "decoys": manifest}, f, indent=2, sort_keys=True)
    print(f"prepared {ok}/{len(manifest)}; coordinates lost in {zero} entries; "
          f"no polar H in {noH}")
    if bad:
        print("failures:", bad[:10])
    print(f"-> {os.path.relpath(outdir, pl.ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
