"""build_decoys.py — build a property-matched, topologically dissimilar decoy set.

Why a decoy set is needed
-------------------------
The redo reports absolute affinities, but an absolute score only means something if the protocol
can tell binders from non-binders. The sulfate self-docking control showed the geometry is
reproduced; this control asks the complementary question - do the candidate ligands rank above
molecules that are not expected to bind?

Source of decoys
----------------
`RDKit .../share/RDKit/Data/NCI/first_5K.smi`, the first 5,000 structures of the NCI Open Database
as distributed with RDKit. It is local, citable and reproducible, so the control does not depend on
a network fetch.

Method (DUD-E style)
--------------------
For each of the 20 unique candidate molecules, decoys are drawn from the NCI pool that
  1. fall inside that candidate's property windows
     (MW +-20%, Crippen logP +-1.0, HBD +-1, HBA +-1, rotatable bonds +-2, TPSA +-30 A^2), and
  2. are topologically dissimilar: maximum Morgan(radius 2, 2048-bit) Tanimoto < 0.30 against
     *every* candidate, not just the one being matched.
The second condition is what makes these decoys rather than analogues: property-matched but not
chemically related. Up to `--per-active` decoys are taken per candidate.

Writes analysis/docking/decoys/decoys.smi and decoys/decoy_manifest.json.

Usage:
    python analysis/scripts/build_decoys.py [--per-active 8] [--tanimoto 0.30]
"""
import argparse
import json
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = os.path.join(ROOT, "analysis", "docking")
NCI = os.path.join(BASE, ".pixi", "envs", "default", "Library", "share",
                   "RDKit", "Data", "NCI", "first_5K.smi")
OUTDIR = os.path.join(BASE, "decoys")
MANIFEST = os.path.join(BASE, "ligands", "ligand_manifest.json")

# per-candidate property windows
TOL = {"mw_pct": 0.20, "logp_abs": 1.0, "hbd_abs": 1, "hba_abs": 1, "rb_abs": 2, "tpsa_abs": 30.0}


def props(mol, Chem, Descriptors, rdMolDescriptors, AllChem):
    from rdkit.Chem import DataStructs
    return {"mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "rb": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol),
            "charge": Chem.GetFormalCharge(mol)}


def main():
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    ap = argparse.ArgumentParser()
    ap.add_argument("--per-active", type=int, default=8)
    ap.add_argument("--tanimoto", type=float, default=0.30)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    man = json.load(open(MANIFEST, encoding="utf-8"))
    actives = {}
    for key, v in man["ligands"].items():
        actives.setdefault(v["name"], v["smiles"])
    print(f"candidates: {len(actives)} unique molecules")

    def fp(mol):
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

    act_mols, act_fps, act_props = {}, {}, {}
    for name, smi in actives.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        act_mols[name] = m
        act_fps[name] = fp(m)
        act_props[name] = props(m, Chem, Descriptors, rdMolDescriptors, AllChem)

    # ---- the NCI pool ------------------------------------------------------------------------
    pool = []
    for ln in open(NCI, encoding="utf-8", errors="ignore"):
        parts = ln.split()
        if len(parts) < 2:
            continue
        m = Chem.MolFromSmiles(parts[0])
        if m is None or m.GetNumHeavyAtoms() < 6:
            continue
        pool.append({"nci_id": parts[1], "smiles": Chem.MolToSmiles(m),
                     "mol": m, "fp": fp(m),
                     "props": props(m, Chem, Descriptors, rdMolDescriptors, AllChem)})
    print(f"NCI pool: {len(pool)} parseable structures")

    act_canon = {Chem.MolToSmiles(m) for m in act_mols.values()}
    rng = random.Random(a.seed)

    chosen, used = {}, set()
    for name in sorted(act_mols):
        p = act_props[name]
        cands = []
        for d in pool:
            if d["nci_id"] in used or d["smiles"] in act_canon:
                continue
            q = d["props"]
            if abs(q["mw"] - p["mw"]) > TOL["mw_pct"] * max(p["mw"], 1):
                continue
            if abs(q["logp"] - p["logp"]) > TOL["logp_abs"]:
                continue
            if abs(q["hbd"] - p["hbd"]) > TOL["hbd_abs"] or abs(q["hba"] - p["hba"]) > TOL["hba_abs"]:
                continue
            if abs(q["rb"] - p["rb"]) > TOL["rb_abs"] or abs(q["tpsa"] - p["tpsa"]) > TOL["tpsa_abs"]:
                continue
            if q["charge"] != p["charge"]:
                continue
            mx = max(DataStructs.TanimotoSimilarity(d["fp"], f) for f in act_fps.values())
            if mx >= a.tanimoto:
                continue
            cands.append((d, mx))
        rng.shuffle(cands)
        take = sorted(cands[:a.per_active], key=lambda x: x[1])
        for d, mx in take:
            used.add(d["nci_id"])
            chosen[d["nci_id"]] = {"nci_id": d["nci_id"], "smiles": d["smiles"],
                                   "matched_active": name,
                                   "max_tanimoto_to_any_active": round(mx, 3),
                                   **{k: (round(v, 3) if isinstance(v, float) else v)
                                      for k, v in d["props"].items()}}
        print(f"  {name:27s} {len(take):2d} decoys "
              f"(pool candidates {len(cands)}), maxT {min([t[1] for t in take], default=float('nan')):.2f}")

    decoys = list(chosen.values())
    print(f"\ntotal decoys: {len(decoys)}")
    if not decoys:
        print("no decoys passed the filters - relax --tanimoto or widen TOL")
        return 1

    smi_path = os.path.join(OUTDIR, "decoys.smi")
    with open(smi_path, "w", encoding="utf-8", newline="\n") as f:
        for d in sorted(decoys, key=lambda x: x["nci_id"]):
            f.write(f"{d['smiles']}\tDEC_{d['nci_id']}\n")

    # summary of property match quality
    def rng_of(key, actives_only):
        vals = [act_props[n][key] for n in act_props] if actives_only else [d[key] for d in decoys]
        return (min(vals), max(vals))

    with open(os.path.join(OUTDIR, "decoy_manifest.json"), "w", encoding="utf-8") as f:
        json.dump({"source": os.path.relpath(NCI, ROOT).replace("\\", "/"),
                   "source_note": "first 5,000 structures of the NCI Open Database as shipped with RDKit",
                   "n_pool": len(pool), "n_decoys": len(decoys),
                   "per_active_cap": a.per_active, "tanimoto_cutoff": a.tanimoto,
                   "tolerance": TOL, "seed": a.seed,
                   "active_property_ranges": {k: [round(x, 2) for x in rng_of(k, True)]
                                              for k in ("mw", "logp", "hbd", "hba", "rb", "tpsa")},
                   "decoy_property_ranges": {k: [round(x, 2) for x in rng_of(k, False)]
                                             for k in ("mw", "logp", "hbd", "hba", "rb", "tpsa")},
                   "decoys": decoys}, f, indent=2)
    print(f"actives range MW {rng_of('mw', True)[0]:.0f}-{rng_of('mw', True)[1]:.0f}, "
          f"logP {rng_of('logp', True)[0]:.1f}-{rng_of('logp', True)[1]:.1f}")
    print(f"decoys  range MW {rng_of('mw', False)[0]:.0f}-{rng_of('mw', False)[1]:.0f}, "
          f"logP {rng_of('logp', False)[0]:.1f}-{rng_of('logp', False)[1]:.1f}")
    print(f"-> {os.path.relpath(smi_path, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
