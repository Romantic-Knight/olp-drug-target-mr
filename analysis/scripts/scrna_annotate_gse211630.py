# -*- coding: utf-8 -*-
"""Stage 4 (option B): formal cell-type annotation of OLP Leiden clusters +
candidate-gene (GDF11/SUOX/EGFL8/TAPBP/HLA-DPA1) x cell-type enrichment tests.

Evidence basis per cluster:
  - empirical DE top genes (wilcoxon, group vs rest; 40 genes/cluster,
    analysis/scripts/_probe_annotation_gse211630.py -> stage4_probe/)
  - curated lineage markers present in the merged OLP 17,669-gene universe
    (KRT13/14/17/19, CD79A/B, COL1A1, PECAM1, CLDN5, NKG7, TRAC, CCL5 etc. are
    ABSENT because NEOLP-1/EOLP-2 lack those rows - see marker_presence.json;
    absent canonical markers were substituted by present alternatives)
  - QC probe for doublet suspicion (clusters 12/18: n=16, single sample,
    cross-lineage co-expression)
Final names (sub), coarse family (family), confidence, evidence and note are
recorded here so the annotation is fully auditable.

Outputs (results/tables + data/processed/geo_GSE211630/scanpy):
  - geo_gse211630_celltype_annotation.csv       cluster-level annotation table
  - gse211630_olp_obs_annotation.csv            per-cell obs annotation (barcode-indexed)
  - geo_gse211630_candidate_celltype_enrichment.csv   coarse family x 5 genes
  - geo_gse211630_candidate_cluster_enrichment.csv    Leiden cluster x 5 genes
  - geo_gse211630_candidate_celltype_expr.csv         frac/mean per family per gene
"""
import json
import os

import numpy as np
import pandas as pd
import scipy.stats as st
import scanpy as sc

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
SCAN = os.path.join("data", "processed", "geo_GSE211630", "scanpy")
TAB = os.path.join("results", "tables")
os.makedirs(TAB, exist_ok=True)

# ---------------------------------------------------------------- annotation
# cluster -> dict; 'family' values used as coarse cell-type groups
ANNOT = {
    0:  dict(sub="Fibroblast (SFRP2+/COL1A2+)",      family="Fibroblast",  conf="high",
             ev="COL1A2/DCN/LUM/COL3A1/COL6A3/SFRP2/VCAN/C1S/C1R/CTSK",
             note="lamina propria fibroblast, matrix-high"),
    1:  dict(sub="Fibroblast (CXCL14+/APOD+)",       family="Fibroblast",  conf="medium",
             ev="DCN/CXCL14/MEG3/APOD/SOD3/MGP/TWIST1/IGFBP5/CEBPD",
             note="CXCL14+ fibroblast subset (oral-cavity fibroblast signature); "
                  "no pericyte/SMC markers"),
    2:  dict(sub="Pericyte/SMC",                     family="Pericyte/SMC", conf="high",
             ev="TAGLN/ACTA2/CALD1/TPM1/TPM2/MCAM/MYLK/IGFBP7/RGS5/MYH11/PDGFRB",
             note="vascular smooth-muscle/pericyte"),
    3:  dict(sub="T cell (resting/recirculating, CD4/CD8)", family="T/NK", conf="medium",
             ev="CD3D/CD3E/CD2/IL7R high; CD4+CD8A both moderate; CXCR4/CREM/NR4A2/BTG1",
             note="largest T compartment; mixed CD4/CD8, recirculating-like; "
                  "minor GNLY+ cells present (NK cannot be separated at res 0.5)"),
    4:  dict(sub="T/NK (stressed, HSP-high)",        family="T/NK",       conf="low",
             ev="HSPA1A/HSPA1B/HSPA6/HSP90AA1 + CD3D/CD2/IL7R + GNLY/KLRD1/IFNG",
             note="dissociation-stress high; CD3 (T) and GNLY (NK) co-mixed -> "
                  "heterogeneous, lineage-level only"),
    5:  dict(sub="T/NK-like (low-complexity)",       family="T/NK",       conf="low",
             ev="BCL11B/NKTR/FYB1/MALAT1 + CD2; CD3 moderate",
             note="low UMI/genes (med n_counts 861); possible low-quality lymphocytes"),
    6:  dict(sub="CD8+ T (GZMK+ effector-memory)",   family="T/NK",       conf="high",
             ev="CD3D/CD3E/CD2/CD8A/CD8B/GZMK/GZMA/IL32/GIMAP7",
             note="CD8+ effector-memory/cytotoxic T"),
    7:  dict(sub="CD4+ Treg (FOXP3+)",               family="T/NK",       conf="high",
             ev="FOXP3/CTLA4/TIGIT/CD27/LTB/IL32/CD3D",
             note="regulatory T cells (Treg)"),
    8:  dict(sub="T cell (naive-like, IL7Rhigh)",    family="T/NK",       conf="medium",
             ev="IL7R/CD3D/CD3E/CD2/CD52/GIMAP7/LTB; SELL+",
             note="naive/central-memory-like; CD4/CD8 mixed (subset unresolved)"),
    9:  dict(sub="Plasmablast (GZMB+)",              family="B/Plasma",   conf="medium-high",
             ev="MZB1/JCHAIN/FCRL5/BANK1(+)/BCL11A/IRF7/PLD4/GZMB; HLA-DR+CD74+",
             note="antibody-secreting precursor; MKI67 low (not cycling); "
                  "GZMB+ plasmablast phenotype"),
    10: dict(sub="Macrophage (LYZ+/HLA-DR+)",        family="Myeloid",    conf="high",
             ev="LYZ/CD74/HLA-DRA/DPA1/DRB1/AIF1/SPI1/LST1/FCER1G; CD14+",
             note="macrophage-dominant myeloid; DC subsets not resolved"),
    11: dict(sub="Epithelial keratinocyte (basal/proliferating)", family="Epithelial",
             conf="high",
             ev="KRT5/KRT6A/KRT6B/SFN/DSP/PERP/S100A2 + MKI67/CENPF",
             note="proliferating basal keratinocyte (only epithelial cluster; "
                  "suprabasal KRT13/4 markers absent from universe)"),
    12: dict(sub=None, family=None, conf=None,
             ev="KRT5/KRT6A + LYZ/S100A8 + CD3D/CD8B co-detected",
             note="n=16, NEOLP-3 only; cross-lineage co-expression -> doublet-suspect"),
    13: dict(sub="Mast cell",                        family="Mast",       conf="high",
             ev="TPSAB1/TPSB2/CPA3/CTSG/MS4A2/HPGDS/IL1RL1/GATA2/KIT",
             note="mast cells (tryptase+)"),
    14: dict(sub="B cell (MS4A1+)",                  family="B/Plasma",   conf="high",
             ev="MS4A1/BANK1/CD19/HLA-DRA/DPA1/DQB1/CD74/SWAP70",
             note="mature naive/memory B; CD79A/B absent from universe"),
    15: dict(sub="Plasma cell (SDC1+)",              family="B/Plasma",   conf="high",
             ev="MZB1/JCHAIN/SDC1/TNFRSF17/SSR4/FKBP11/ER-machinery; FCRL5",
             note="mature antibody-secreting plasma cells"),
    16: dict(sub="CD8+ T (GZMB+/CXCL13+ cytotoxic-exhausted-like)", family="T/NK",
             conf="medium-high",
             ev="CD8A/CD8B/CD3D/GZMB/CXCL13/LAG3/LINC01871",
             note="cytotoxic/exhausted-like CD8 (NEOLP-2/NEOLP-3 enriched); "
                  "MKI67 moderate 0.17 frac"),
    17: dict(sub="Endothelial (lymphatic, CCL21+)",  family="Endothelial", conf="medium-high",
             ev="CCL21/PROX1/NRP2/MMRN1/TFPI/TM4SF1; VWF low",
             note="lymphatic endothelial cells"),
    18: dict(sub=None, family=None, conf=None,
             ev="CD3D/FOXP3 + LYZ/CXCL8/CTSB + MS4A1/BANK1 co-detected",
             note="n=16, NEOLP-2 only; cross-lineage co-expression, elevated "
                  "doublet score -> doublet-suspect"),
    19: dict(sub="Endothelial (blood, ENG+/ADGRL4+)", family="Endothelial", conf="high",
             ev="IFI27/ENG/ADGRL4/RNASE1/CALCRL/GNG11/TM4SF1/SPARCL1",
             note="blood vascular endothelial cells"),
    20: dict(sub="Myeloid (S100-high inflammatory)", family="Myeloid",   conf="low",
             ev="S100A8/S100A9/CXCL8/NAMPT/SAT1/PLEK + CD14/AIF1/FCGR3A",
             note="low transcript complexity (med n_genes 473); monocyte/"
                  "neutrophil-like inflammatory cells; CSF3R/FCGR3B/S100A12 "
                  "absent from universe so granulocyte identity unverifiable"),
}
CANDIDATES = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]

print("loading OLP lognorm object ...", flush=True)
olp = sc.read_h5ad(os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"))
leid = "leiden_0.5"
olp.obs[leid] = olp.obs[leid].astype(int)

# sanity: annotation covers all clusters
obs_cl = set(int(x) for x in olp.obs[leid].unique())
assert obs_cl == set(ANNOT.keys()), (sorted(obs_cl), sorted(ANNOT.keys()))

# ------------------------------------------------------- per-cell annotation
sub = olp.obs[leid].map(lambda c: ANNOT[int(c)]["sub"])
fam = olp.obs[leid].map(lambda c: ANNOT[int(c)]["family"])
cfd = olp.obs[leid].map(lambda c: ANNOT[int(c)]["conf"])
olp.obs["celltype_sub"] = pd.Categorical(sub)
olp.obs["celltype_family"] = pd.Categorical(fam)
olp.obs["annotation_conf"] = cfd.values
olp.obs["doublet_suspect"] = fam.isna().values
annotated = ~olp.obs["doublet_suspect"].values
print(f"cells annotated={int(annotated.sum())} doublet_suspect={int((~annotated).sum())}",
      flush=True)

# per-cell obs table (barcode-indexed) for downstream plotting/QA
ob = olp.obs[["sample", "condition", leid, "celltype_sub", "celltype_family",
              "annotation_conf", "doublet_suspect"]].copy()
ob = ob.rename(columns={leid: "leiden_0.5"})
ob.index.name = "barcode"
ob.to_csv(os.path.join(SCAN, "gse211630_olp_obs_annotation.csv"))

# -------------------------------------------------- cluster annotation table
rows = []
n_tot = olp.obs.shape[0]
for cl in sorted(ANNOT):
    m = (olp.obs[leid].values == cl)
    n = int(m.sum())
    a = ANNOT[cl]
    rows.append({
        "leiden_0.5": cl, "n_cells": n, "pct_of_OLP": round(100 * n / n_tot, 2),
        "celltype_family": a["family"], "celltype_sub": a["sub"],
        "annotation_confidence": a["conf"], "evidence_markers": a["ev"],
        "note": a["note"],
        **{f"frac_{s}": round(float((olp.obs.loc[m, 'sample'] == s).mean()), 4)
           for s in ["normal-1", "NEOLP-1", "NEOLP-2", "NEOLP-3", "EOLP-1", "EOLP-2"]},
    })
pd.DataFrame(rows).to_csv(os.path.join(TAB, "geo_gse211630_celltype_annotation.csv"),
                          index=False)

# -------------------------------------------------- expression extraction
Xc = olp[:, CANDIDATES].X  # sparse (46,478 x 5)
if not hasattr(Xc, "todense"):
    pass
Xc_d = np.asarray(Xc.todense()) if hasattr(Xc, "todense") else np.asarray(Xc)
# quick check: all 5 candidates present in OLP universe
assert all(g in olp.var_names for g in CANDIDATES), [g for g in CANDIDATES
                                                     if g not in olp.var_names]

expr_mask = Xc_d > 0

# ------------------------------------------------------------- BH correction
def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    po = p[ok]
    idx = np.argsort(po)
    m = len(po)
    raw = np.minimum(po[idx] * (m / np.arange(1, m + 1)), 1.0)
    qs = np.minimum.accumulate(raw[::-1])[::-1]
    qtmp = np.empty(m)
    qtmp[idx] = qs
    q[ok] = qtmp
    return q

# ------------------------------------------------------- coarse-family tests
famv = olp.obs["celltype_family"].values
fam_list = [f for f in FAMILY_ORDER if (famv == f).any()]
expr_fam = pd.DataFrame({  # per family per gene summary
    "gene": [], "family": [], "n_family": [], "n_expr": [], "frac_expr": [],
    "mean_log1p_all": [], "mean_log1p_expr": []})
enr_rows = []
for gi, g in enumerate(CANDIDATES):
    e = expr_mask[:, gi]
    for f in fam_list:
        mf = (famv == f) & annotated
        nf = int(mf.sum())
        ne = int((e & mf).sum())
        frac = ne / nf
        vals = Xc_d[mf, gi]
        ml_all = float(vals.mean())
        ml_expr = float(vals[e[mf]].mean()) if ne else np.nan
        # Fisher (one-sided greater: expression enriched in this family)
        nrest = int((annotated & ~(famv == f)).sum())
        nrest_e = int((e & annotated & ~(famv == f)).sum())
        table = [[ne, nrest_e], [nf - ne, nrest - nrest_e]]
        orr, p = st.fisher_exact(table, alternative="greater")
        p_log10 = float(-np.inf) if p == 0 else round(float(np.log10(p)), 2)
        enr_rows.append({"gene": g, "family": f, "n_family": nf, "n_expr": ne,
                         "frac_expr": round(frac, 4),
                         "mean_log1p_all": round(ml_all, 4),
                         "mean_log1p_expr": round(ml_expr, 4) if ne else np.nan,
                         "fold_vs_rest_frac": round(frac / max(nrest_e / nrest, 1e-9), 3),
                         "odds_ratio": round(float(orr), 3) if orr != np.inf else np.inf,
                         "fisher_p": p, "fisher_log10p": p_log10})
enr = pd.DataFrame(enr_rows)
enr["fisher_FDR_BH"] = bh(enr["fisher_p"].values)
enr = enr.sort_values(["fisher_p", "gene"])
enr.to_csv(os.path.join(TAB, "geo_gse211630_candidate_celltype_enrichment.csv"),
           index=False)

# Kruskal-Wallis across families per gene (global heterogeneity, all annotated cells)
kw_rows = []
for g in CANDIDATES:
    vals = Xc_d[annotated, CANDIDATES.index(g)]
    grp = famv[annotated]
    data = [vals[grp == f] for f in fam_list]
    H, p = st.kruskal(*data)
    kw_rows.append({"gene": g, "kruskal_H": round(float(H), 2), "kruskal_p": p,
                    "n_families": len(fam_list)})
pd.DataFrame(kw_rows).to_csv(os.path.join(TAB, "geo_gse211630_candidate_kruskal.csv"),
                             index=False)

# -------------------------------------------------------- cluster-level tests
clv = olp.obs[leid].values
cl_list = sorted(ANNOT)
cl_rows = []
for gi, g in enumerate(CANDIDATES):
    e = expr_mask[:, gi]
    for cl in cl_list:
        mc = (clv == cl) & annotated
        nc = int(mc.sum())
        ne = int((e & mc).sum())
        frac = ne / nc if nc else np.nan
        nrest = int((annotated & ~(clv == cl)).sum())
        nrest_e = int((e & annotated & ~(clv == cl)).sum())
        orr, p = st.fisher_exact([[ne, nrest_e], [nc - ne, nrest - nrest_e]],
                                 alternative="greater")
        p_log10 = float(-np.inf) if p == 0 else round(float(np.log10(p)), 2)
        cl_rows.append({"gene": g, "leiden_0.5": cl,
                        "celltype_sub": ANNOT[cl]["sub"],
                        "n_cluster": nc, "n_expr": ne,
                        "frac_expr": round(frac, 4) if nc else np.nan,
                        "odds_ratio": round(float(orr), 3) if orr != np.inf else np.inf,
                        "fisher_p": p, "fisher_log10p": p_log10})
clenr = pd.DataFrame(cl_rows)
clenr["fisher_FDR_BH"] = bh(clenr["fisher_p"].values)
clenr.to_csv(os.path.join(TAB, "geo_gse211630_candidate_cluster_enrichment.csv"),
             index=False)

# ---------------------------------------------------- headline summary json
summ = {}
for g in CANDIDATES:
    sub = enr[enr["gene"] == g].sort_values("frac_expr", ascending=False)
    top3 = [{"family": r.family, "frac_expr": r.frac_expr,
             "mean_log1p_expr": r.mean_log1p_expr,
             "fisher_FDR_BH": r.fisher_FDR_BH,
             "fisher_p": r.fisher_p} for r in sub.head(3).itertuples()]
    kw = kw_rows[CANDIDATES.index(g)]
    summ[g] = {"top3_families_by_frac": top3,
               "kruskal_H": kw["kruskal_H"], "kruskal_p": kw["kruskal_p"],
               "overall_frac_annotated": round(float(expr_mask[annotated, CANDIDATES.index(g)].mean()), 4),
               "overall_mean_log1p_all": round(float(Xc_d[annotated, CANDIDATES.index(g)].mean()), 4)}
def _json_default(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError(f"Object of type {o.__class__.__name__} is not JSON serializable")


with open(os.path.join(SCAN, "candidate_celltype_summary.json"), "w") as f:
    json.dump(summ, f, indent=1, default=_json_default)

# quick console digest
print("\n== per-family digest (frac_expr; mean_log1p_expr; FDR) ==", flush=True)
for g in CANDIDATES:
    sub = enr[enr["gene"] == g]
    s = sub.sort_values("frac_expr", ascending=False)
    line = " | ".join(f"{r.family}:{r.frac_expr:.3f}({r.mean_log1p_expr:.2f},"
                      f"FDR{r.fisher_FDR_BH:.2g})" for r in s.itertuples())
    print(f"{g}: {line}", flush=True)
print("DONE", flush=True)
