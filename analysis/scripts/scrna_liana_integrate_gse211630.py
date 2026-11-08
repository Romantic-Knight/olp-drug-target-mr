# -*- coding: utf-8 -*-
"""
Stage D, part 2 (integration): merge per-condition LIANA CellChat + Connectome
results into condition-contrast tables and derive mechanistic readouts.

Inputs (results/tables/):
  geo_gse211630_liana_cellchat_{NEOLP,EOLP}.csv
  geo_gse211630_liana_connectome_{NEOLP,EOLP}.csv
  geo_gse211630_celltype_annotation.csv (cluster -> family/subtype)

Outputs (results/tables/):
  geo_gse211630_liana_signal_shift.csv      union of interactions, per-interaction
                                            contrast (sig flags, prob/scaled weight
                                            per condition, delta, rank pct, gain/loss class)
  geo_gse211630_liana_incoming_focus.csv    focus receivers (T/NK + Fibroblast):
                                            per receiver x sender-family x class counts &
                                            summed prob deltas
  geo_gse211630_liana_receptor_module_corr.csv  T/NK cells: Spearman(receptor expression,
                                            Stage-B module scores) for receptors present in
                                            significant incoming interactions
  geo_gse211630_liana_headlines.csv         top gained / top lost specific L-R pairs on
                                            focus receivers (for narrative + Fig D1)

Caveat (see README): CellChat p-values are computed within condition over
cells pooled across samples -> pseudo-replication; treat as descriptive.
Run: pixi python (needs anndata/scanpy/scipy for the module-score part).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import anndata as ad

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
TBL = os.path.join(ROOT, "results", "tables")
OBJ = os.path.join(ROOT, "data", "processed", "geo_GSE211630", "scanpy",
                   "gse211630_olp_lognorm.h5ad")
ANN = os.path.join(ROOT, "results", "tables",
                   "geo_gse211630_celltype_annotation.csv")
DOUBLET = {12, 18}
FOCUS_T = {3, 4, 5, 6, 7, 8, 16}   # T/NK receivers (Stage B subtypes)
FOCUS_F = {0, 1}                    # Fibroblast receivers
SIG = 0.05

MODULES = {
    "Stress":   ["HSPA1A", "HSPA1B", "HSP90AA1", "HSPH1"],
    "Exhaustion": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX",
                   "ENTPD1", "CXCL13"],
    "Treg":     ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "BATF"],
    "Naive":    ["IL7R", "SELL", "TCF7", "LEF1", "S100A4"],
    "Activation": ["CD69", "TNF", "IFNG", "IL32"],
    "Proliferation": ["MKI67"],
    "Cytotoxic": ["GZMB", "GZMK", "PRF1", "GNLY", "IFNG", "KLRD1"],
}
KEY = ["source", "target", "ligand", "receptor"]


def rd(name):
    p = os.path.join(TBL, name)
    return pd.read_csv(p, encoding="utf-8-sig")


def main():
    ccN, ccE = rd("geo_gse211630_liana_cellchat_NEOLP.csv"), \
        rd("geo_gse211630_liana_cellchat_EOLP.csv")
    coN, coE = rd("geo_gse211630_liana_connectome_NEOLP.csv"), \
        rd("geo_gse211630_liana_connectome_EOLP.csv")
    ann = pd.read_csv(ANN, encoding="utf-8-sig").set_index("leiden_0.5")

    for df in (ccN, ccE, coN, coE):
        for k in KEY:
            df[k] = df[k].fillna("").astype(str)

    # ---- per-interaction contrast -------------------------------------
    def prep_cc(d, tag):
        d = d.copy()
        d = d.rename(columns={"lr_probs": f"prob_{tag}",
                              "cellchat_pvals": f"p_{tag}"})
        d = d[[*KEY, f"prob_{tag}", f"p_{tag}"]]
        d[f"sig_{tag}"] = d[f"p_{tag}"] <= SIG
        return d.drop_duplicates(KEY)

    def prep_co(d, tag):
        d = d.copy()
        d = d.rename(columns={"scaled_weight": f"sw_{tag}",
                              "expr_prod": f"ep_{tag}"})
        d = d[[*KEY, f"sw_{tag}", f"ep_{tag}"]].drop_duplicates(KEY)
        # specificity percentile within condition
        d[f"rank_pct_{tag}"] = d[f"sw_{tag}"].rank(pct=True) * 100
        return d

    aN = prep_cc(ccN, "N").merge(prep_co(coN, "N"), on=KEY, how="outer")
    aE = prep_cc(ccE, "E").merge(prep_co(coE, "E"), on=KEY, how="outer")
    u = pd.concat([aN[[*KEY]], aE[[*KEY]]], ignore_index=True).drop_duplicates(KEY)
    sh = u.merge(aN, on=KEY, how="left").merge(aE, on=KEY, how="left",
                                               suffixes=("", "_y"))
    sh = sh[[c for c in sh.columns if not c.endswith("_y")]]
    sh["present_N"] = sh["prob_N"].notna()
    sh["present_E"] = sh["prob_E"].notna()
    sh["sig_N"] = sh["sig_N"].fillna(False)
    sh["sig_E"] = sh["sig_E"].fillna(False)
    sh["prob_N"] = sh["prob_N"].fillna(0.0)
    sh["prob_E"] = sh["prob_E"].fillna(0.0)
    sh["sw_N"] = sh["sw_N"].fillna(0.0)
    sh["sw_E"] = sh["sw_E"].fillna(0.0)
    sh["delta_prob"] = sh["prob_E"] - sh["prob_N"]
    sh["delta_sw"] = sh["sw_E"] - sh["sw_N"]
    # classification
    def cls(r):
        if r["sig_E"] and not r["sig_N"]:
            return "EOLP_gained"
        if r["sig_N"] and not r["sig_E"]:
            return "EOLP_lost"
        if r["sig_N"] and r["sig_E"]:
            return "sustained"
        return "none"
    sh["class"] = sh.apply(cls, axis=1)

    # annotation mapping for source/target
    ann2 = (ann.reset_index()[["leiden_0.5", "celltype_family", "celltype_sub"]]
            .dropna(subset=["celltype_sub"]))
    ann2.columns = ["cluster", "family", "subtype"]
    sub2cluster = dict(zip(ann2["subtype"], ann2["cluster"]))
    sub2family = dict(zip(ann2["subtype"], ann2["family"]))
    def map_side(df, side):
        col = side + "_cluster"
        df[col] = df[side].map(sub2cluster)
        df[side + "_family"] = df[side].map(sub2family)
        return df
    sh = map_side(sh, "source")
    sh = map_side(sh, "target")
    sh["receiver_focus"] = sh["target_cluster"].isin(FOCUS_T | FOCUS_F)

    cols = [*KEY, "source_cluster", "source_family", "target_cluster",
            "target_family", "receiver_focus", "present_N", "present_E",
            "prob_N", "prob_E", "p_N", "p_E", "sig_N", "sig_E",
            "sw_N", "sw_E", "rank_pct_N", "rank_pct_E",
            "delta_prob", "delta_sw", "class"]
    sh = sh[cols].sort_values(["receiver_focus", "class", "target_cluster",
                               "source_cluster"], ascending=[False, True, True, True])
    sh.to_csv(os.path.join(TBL, "geo_gse211630_liana_signal_shift.csv"),
              index=False, encoding="utf-8-sig")
    print("signal_shift rows:", len(sh),
          "| focus:", int(sh["receiver_focus"].sum()))

    # ---- incoming summary (focus receivers) ----------------------------
    foc = sh[sh["receiver_focus"]]
    summ = (foc.groupby(["target_cluster", "target", "source_family", "class"])
            .agg(n=("ligand", "size"),
                 sum_delta_prob=("delta_prob", "sum"),
                 sum_delta_sw=("delta_sw", "sum"))
            .reset_index())
    summ.to_csv(os.path.join(TBL, "geo_gse211630_liana_incoming_focus.csv"),
                index=False, encoding="utf-8-sig")

    # ---- headline pairs ------------------------------------------------
    def head_rows(df, cls_name, how, n=20):
        h = df[df["class"] == cls_name]
        if how == "prob":
            h = h.sort_values("delta_prob", ascending=False)
        else:
            h = h.sort_values("delta_sw", ascending=False)
        return h.head(n)
    hn = head_rows(foc, "EOLP_gained", "prob")
    hl = head_rows(foc, "EOLP_lost", "prob", n=20)
    hl = hl.sort_values("delta_prob", ascending=True)
    hl = pd.concat([hn, hl], ignore_index=True)
    hl.to_csv(os.path.join(TBL, "geo_gse211630_liana_headlines.csv"),
              index=False, encoding="utf-8-sig")

    # ---- receptor x module Spearman on T/NK cells ----------------------
    a = ad.read_h5ad(OBJ)
    a.obs["leiden_0.5"] = a.obs["leiden_0.5"].astype(int)
    tnk = a[a.obs["leiden_0.5"].isin(FOCUS_T)].copy()
    print("T/NK cells for module corr:", tnk.n_obs)
    # receptors from significant incoming interactions onto T/NK receivers
    sig_rx = sh[(sh["sig_N"] | sh["sig_E"]) &
                (sh["target_cluster"].isin(FOCUS_T))]
    rx_counts = (sig_rx["receptor"].value_counts())
    genes = [g for g in rx_counts.index if g in a.var_names][:150]
    mod_genes = {m: [g for g in gs if g in a.var_names]
                 for m, gs in MODULES.items()}
    used = sorted(set(genes) | {g for gs in mod_genes.values() for g in gs})
    X = tnk[:, used].X  # sparse
    if hasattr(X, "toarray"):
        X = X.toarray()
    Xdf = pd.DataFrame(X, columns=used, index=tnk.obs_names)
    mod_scores = pd.DataFrame({m: Xdf[gs].mean(axis=1)
                               for m, gs in mod_genes.items() if gs})
    rows = []
    for g in genes:
        if Xdf[g].std() == 0:
            continue
        for m in mod_scores.columns:
            rho, p = spearmanr(Xdf[g], mod_scores[m])
            rows.append((g, m, rho, p, int(rx_counts[g])))
    corr = pd.DataFrame(rows, columns=["receptor", "module", "rho", "p",
                                       "n_incoming_interactions"])
    corr.to_csv(os.path.join(TBL, "geo_gse211630_liana_receptor_module_corr.csv"),
                index=False, encoding="utf-8-sig")
    print("receptor-module corr rows:", len(corr))


if __name__ == "__main__":
    main()
