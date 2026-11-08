# -*- coding: utf-8 -*-
"""
Stage D, part 1 (per condition): run LIANA CellChat + Connectome on the
GSE211630 OLP universe (46,446 annotated cells after dropping doublet-suspect
clusters 12/18) for ONE condition given as argv[1] (NEOLP | EOLP).

Grouping: 19 cell subtypes (mapped from annotation CSV via leiden_0.5).
Expression: adata.X = log1p-normalized (use_raw=False).
Resource : OmniPath 'consensus' (4,624 L-R pairs; ~61% of entities present).
CellChat: n_perms=1000, seed=1337 (default); Connectome deterministic.

Outputs per condition (results/tables/):
  geo_gse211630_liana_cellchat_{COND}.csv
  geo_gse211630_liana_connectome_{COND}.csv

Note: within-condition permutation p-values are computed over cells pooled
across the condition's samples (pseudo-replication); these tables are treated
as descriptive/mechanistic hypothesis generation, never as sample-level
evidence (see README Stage D).
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
from liana import mt

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
OBJ = os.path.join(ROOT, "data", "processed", "geo_GSE211630", "scanpy",
                   "gse211630_olp_lognorm.h5ad")
ANN = os.path.join(ROOT, "results", "tables",
                   "geo_gse211630_celltype_annotation.csv")
TBL = os.path.join(ROOT, "results", "tables")
DOUBLET = {12, 18}
RS = 1337
N_PERMS = 1000
EXPR_PROP = 0.1
MIN_CELLS = 5


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else None
    if cond not in ("NEOLP", "EOLP"):
        raise SystemExit("usage: python scrna_liana_run_cond_gse211630.py "
                         "NEOLP|EOLP")

    t0 = time.time()
    a = ad.read_h5ad(OBJ)
    ann = pd.read_csv(ANN, encoding="utf-8-sig").set_index("leiden_0.5")
    a.obs["leiden_0.5"] = a.obs["leiden_0.5"].astype(int)
    a.obs["family"] = a.obs["leiden_0.5"].map(ann["celltype_family"])
    a.obs["subtype"] = a.obs["leiden_0.5"].map(ann["celltype_sub"])
    drop = a.obs["leiden_0.5"].isin(DOUBLET)
    a = a[~drop].copy()  # 46,446 annotated cells
    assert a.n_obs == 46446, a.n_obs
    # annotate focus receivers (T/NK + Fibroblast hosts of risk-gene loss)
    t_focus = a.obs["leiden_0.5"].isin(
        {0, 1, 3, 4, 5, 6, 7, 8, 16})
    a.obs["receiver_focus"] = t_focus
    print(f"[{cond}] load+annot {time.time()-t0:.1f}s universe={a.n_obs}")

    c = a[a.obs["condition"] == cond].copy()
    print(f"[{cond}] n_cells={c.n_obs} subtypes={c.obs['subtype'].nunique()}")

    t0 = time.time()
    cc = mt.cellchat(c, groupby="subtype", use_raw=False, verbose=False,
                     inplace=False, expr_prop=EXPR_PROP, min_cells=MIN_CELLS,
                     n_perms=N_PERMS, seed=RS, n_jobs=1)
    cc.to_csv(os.path.join(TBL, f"geo_gse211630_liana_cellchat_{cond}.csv"),
              index=False, encoding="utf-8-sig")
    print(f"[{cond}] cellchat done {time.time()-t0:.0f}s rows={len(cc)}")

    t0 = time.time()
    co = mt.connectome(c, groupby="subtype", use_raw=False, verbose=False,
                       inplace=False, expr_prop=EXPR_PROP, min_cells=MIN_CELLS)
    co.to_csv(os.path.join(TBL, f"geo_gse211630_liana_connectome_{cond}.csv"),
              index=False, encoding="utf-8-sig")
    print(f"[{cond}] connectome done {time.time()-t0:.0f}s rows={len(co)}")


if __name__ == "__main__":
    main()
