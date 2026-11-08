# -*- coding: utf-8 -*-
"""Rebuild all-genes log-normalised objects with Stage-3 embeddings/clusters.

The analysed h5ad written by scanpy_cluster_gse211630.py were HVG-subset
(4,322 / 5,874 vars) because PCA/Leiden/UMAP operate on the HVG subspace.
Candidate-gene localisation and any downstream DE must read expression over the
*whole* gene universe, so here we rebuild the full-var log-normalised objects
(identical cell order to the analysed objects) and copy obs/obsm/uns across.

Outputs: data/processed/geo_GSE211630/scanpy/
  gse211630_all_lognorm.h5ad    60,308 x 12,021 (6 samples; normal-1 universe)
  gse211630_olp_lognorm.h5ad    46,478 x 17,669 (5 OLP samples)
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
OUT = os.path.join("data", "processed", "geo_GSE211630")
ANN = os.path.join(OUT, "ann")
SCAN = os.path.join(OUT, "scanpy")

SAMPLES = {
    "GSM6481637": "normal-1", "GSM6481638": "NEOLP-1",
    "GSM6481639": "NEOLP-2", "GSM6481640": "NEOLP-3",
    "GSM6481641": "EOLP-1", "GSM6481642": "EOLP-2",
}
ALL_ORDER = list(SAMPLES)
OLP_ORDER = list(SAMPLES)[1:]  # GSM6481638..42


def rebuild(order, analyzed_path, out_path, tag):
    print(f"[{tag}] loading per-sample counts + log-normalising", flush=True)
    parts = []
    for gsm in order:
        a = sc.read_h5ad(os.path.join(ANN, f"{gsm}_filtered.h5ad"))
        a.obs_names = [f"{SAMPLES[gsm]}_{i}" for i in range(a.n_obs)]
        a.obs["gsm"] = gsm
        a.obs["sample"] = SAMPLES[gsm]
        a.obs["condition"] = "normal" if gsm == "GSM6481637" else (
            "NEOLP" if "NEOLP" in SAMPLES[gsm] else "EOLP")
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        parts.append(a)
    import anndata as ad
    obj = ad.concat(parts, join="inner")
    print(f"[{tag}] concat {obj.shape}", flush=True)

    ana = sc.read_h5ad(analyzed_path)
    same = list(obj.obs_names) == list(ana.obs_names)
    print(f"[{tag}] cell order identical to analysed object: {same}", flush=True)
    if not same:
        raise RuntimeError("cell order mismatch — embeddings cannot be copied")
    obj.obs = ana.obs.copy()  # overwrite with analysed obs (leiden cols etc.)
    obj.obsm = ana.obsm
    obj.uns = ana.uns
    obj.write_h5ad(out_path, compression="gzip")
    print(f"[{tag}] saved {out_path} ({obj.shape})", flush=True)


rebuild(ALL_ORDER,
        os.path.join(SCAN, "gse211630_all_analyzed.h5ad"),
        os.path.join(SCAN, "gse211630_all_lognorm.h5ad"), "all")
rebuild(OLP_ORDER,
        os.path.join(SCAN, "gse211630_olp_analyzed.h5ad"),
        os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"), "olp")
print("DONE", flush=True)
