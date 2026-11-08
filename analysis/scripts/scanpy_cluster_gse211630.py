# -*- coding: utf-8 -*-
"""Stage 3: scanpy clustering pipeline for GSE211630 (cell-type localisation of
the 5 MR candidate genes).

Strategy (per user's design): normal-1 is kept on a separate validation track;
the 5 OLP samples are integrated on the *union of per-sample highly variable
genes* (common-HVG integration). Two analysed objects are produced:
  1. all-samples object   (6 samples, common-gene universe = normal-1 genes) —
     TAPBP / HLA-DPA1 / EGFL8 only (GDF11 & SUOX rows absent from normal-1)
  2. OLP-only object      (NEOLP-1/2/3 + EOLP-1/2) — all 5 candidate genes present
Pipeline per object: per-sample log-normalisation -> per-sample HVG (top 2000,
flavor 'seurat') -> union subset -> PCA (arpack) -> neighbours -> Leiden
(res=0.5 & 1.0) -> UMAP. Outputs carry X = log1p(normalize_total(1e4)) counts;
raw-count objects remain the Stage-2 deliverables.
"""
import json
import os
import time

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
OUT = os.path.join("data", "processed", "geo_GSE211630")
ANN = os.path.join(OUT, "ann")
SCAN = os.path.join(OUT, "scanpy")
os.makedirs(SCAN, exist_ok=True)

CANDIDATES = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
SAMPLES = {
    "GSM6481637": ("normal-1", "normal"),
    "GSM6481638": ("NEOLP-1", "NEOLP"),
    "GSM6481639": ("NEOLP-2", "NEOLP"),
    "GSM6481640": ("NEOLP-3", "NEOLP"),
    "GSM6481641": ("EOLP-1", "EOLP"),
    "GSM6481642": ("EOLP-2", "EOLP"),
}
OLP_GSMS = [g for g, (lab, cond) in SAMPLES.items() if cond != "normal"]
N_TOP_HVG = 2000
N_PCS = 50
N_PCS_NEIGH = 30


def t0(label):
    print(f"[{time.strftime('%H:%M:%S')}] {label}", flush=True)


t0("loading 6 per-sample filtered h5ad")
parts = {}
for gsm, (lab, cond) in SAMPLES.items():
    a = sc.read_h5ad(os.path.join(ANN, f"{gsm}_filtered.h5ad"))
    a.obs_names = [f"{lab}_{i}" for i in range(a.n_obs)]  # unique, sample-prefixed
    a.obs["gsm"] = gsm
    a.obs["sample"] = lab
    a.obs["condition"] = cond
    parts[gsm] = a
    print(f"  {lab}: {a.shape}", flush=True)

# per-sample log-normalisation + HVG (common-HVG integration route)
t0("per-sample normalize + HVG")
hvg_union = set()
for gsm, a in parts.items():
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_TOP_HVG, flavor="seurat",
                                subset=False)
    n = int(a.var["highly_variable"].sum())
    hvg_union |= set(a.var_names[a.var["highly_variable"]])
    print(f"  {SAMPLES[gsm][0]}: n_hvg={n}", flush=True)

# concat objects
import anndata as ad
all_ad = ad.concat([parts[g] for g in SAMPLES], join="inner")
print(f"  all: {all_ad.shape}", flush=True)
olp_ad = ad.concat([parts[g] for g in OLP_GSMS], join="inner")
print(f"  olp: {olp_ad.shape}", flush=True)

presence = {}
for name, obj in (("all", all_ad), ("olp", olp_ad)):
    presence[name] = {g: g in obj.var_names for g in CANDIDATES}
    present = [g for g in CANDIDATES if presence[name][g]]
    missing = [g for g in CANDIDATES if not presence[name][g]]
    print(f"  {name}: candidates present={present} missing={missing}", flush=True)


def run_pipeline(obj, tag, res_list=(0.5, 1.0)):
    t0(f"{tag}: HVG union subset")
    obj.var["hvg_common"] = [g in hvg_union for g in obj.var_names]
    obj.var["highly_variable"] = obj.var["hvg_common"]
    print(f"  {tag}: common HVG in object = {int(obj.var['hvg_common'].sum())}", flush=True)
    obj = obj[:, obj.var["hvg_common"]].copy()
    print(f"  {tag}: subset shape {obj.shape}", flush=True)

    t0(f"{tag}: PCA n={N_PCS}")
    sc.pp.pca(obj, n_comps=N_PCS, svd_solver="arpack", random_state=0)

    t0(f"{tag}: neighbours + Leiden + UMAP")
    sc.pp.neighbors(obj, n_neighbors=15, n_pcs=N_PCS_NEIGH, random_state=0)
    for res in res_list:
        sc.tl.leiden(obj, resolution=res, flavor="igraph",
                     n_iterations=2, directed=False,
                     key_added=f"leiden_{res}")
    sc.tl.umap(obj, min_dist=0.3, random_state=0)
    return obj


all_ad = run_pipeline(all_ad, "all")
olp_ad = run_pipeline(olp_ad, "olp")

for obj, tag, fn in ((all_ad, "all", "gse211630_all_analyzed.h5ad"),
                     (olp_ad, "olp", "gse211630_olp_analyzed.h5ad")):
    t0(f"saving {fn}")
    obj.write_h5ad(os.path.join(SCAN, fn), compression="gzip")

# cluster x sample x condition cross-tabs
for obj, tag in ((all_ad, "all"), (olp_ad, "olp")):
    for res in (0.5, 1.0):
        ct = pd.crosstab(obj.obs[f"leiden_{res}"], obj.obs["condition"])
        ct.to_csv(os.path.join(OUT, f"cross_{tag}_leiden{res}_condition.csv"))
        cs = pd.crosstab(obj.obs[f"leiden_{res}"], obj.obs["sample"])
        cs.to_csv(os.path.join(OUT, f"cross_{tag}_leiden{res}_sample.csv"))

with open(os.path.join(SCAN, "candidate_presence.json"), "w") as f:
    json.dump({"candidates": CANDIDATES,
               "presence": presence,
               "n_hvg_union": len(hvg_union),
               "per_sample_n_hvg": {SAMPLES[g][0]: int(parts[g].var["highly_variable"].sum())
                                    for g in SAMPLES},
               "shapes": {"all": list(all_ad.shape), "olp": list(olp_ad.shape)}},
              f, indent=2)
print(f"candidate presence: {presence}", flush=True)
print(f"n_hvg_union = {len(hvg_union)}", flush=True)
for obj, tag in ((all_ad, "all"), (olp_ad, "olp")):
    print(f"{tag}: leiden_0.5 clusters = {obj.obs['leiden_0.5'].nunique()}, "
          f"leiden_1.0 = {obj.obs['leiden_1.0'].nunique()}", flush=True)
print("ALL DONE", flush=True)
