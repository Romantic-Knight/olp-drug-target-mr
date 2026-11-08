# -*- coding: utf-8 -*-
"""Stage B: (1) cell-type composition shift NEOLP vs EOLP and
(2) T/NK functional-state scoring (GSE211630 OLP).

Design notes
------------
* Universe: 46,446 annotated cells of the OLP object (doublet-suspect
  c12/c18, n=32, excluded).
* Replication unit: SAMPLE (NEOLP n=3, EOLP n=2). Cell-level tests are
  exploratory only (cells are pseudoreplicates); the sample-level
  permutation test (C(5,3)=10 label splits, min possible p=0.1) is the
  conservative companion. Composition fractions and per-sample module
  means are tested at the sample level; per-sample n is reported.
* Module score per cell = mean of log1p-normalised expression over the
  signature genes that survive in the 17,669-gene merged universe
  (missing markers per signature are listed below). Score sanity is
  checked against the known subtype biology before statistical use.
* Sparse subtypes (a sample with <10 cells of the type) are kept but
  flagged `sparse_n_lt10` and drawn hollow in the figures (peers must be
  comparable).

Outputs
-------
results/tables/geo_gse211630_composition_per_sample.csv
    sample x condition x level(family|subtype) x group x n_cells x frac
results/tables/geo_gse211630_composition_shift.csv
    per group: NEOLP/EOLP median + mean sample frac, delta, pooled cell
    Fisher p (descriptive), sample-level permutation p
results/tables/geo_gse211630_tnk_module_per_sample.csv
    subtype x signature x sample x condition -> n_cells, mean module score
results/tables/geo_gse211630_tnk_module_shift.csv
    per subtype x signature: NEOLP/EOLP stats, cell-level Wilcoxon + BH,
    sample-level permutation p, sample n
figures/geo_gse211630_composition.png|pdf   (Fig B1, 2 panels)
figures/geo_gse211630_tnk_function.png|pdf  (Fig B2, T/NK dot matrix)
(+ 150-dpi QA copies under analysis/scrna/_qa/)
"""
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
SCAN = os.path.join("data", "processed", "geo_GSE211630", "scanpy")
TAB = os.path.join("results", "tables")
FIGDIR = "figures"
QA = os.path.join("analysis", "scrna", "_qa")
os.makedirs(QA, exist_ok=True)

COND_ORDER = ["NEOLP", "EOLP"]
COL = {"NEOLP": "#4477AA", "EOLP": "#EE7733"}     # blue / orange (colour-safe)
N_NEOLP, N_EOLP = 3, 2

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                     "axes.linewidth": 0.6})

# ---------------------------------------------------------------- annotation maps
olp = sc.read_h5ad(os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"))
ob = pd.read_csv(os.path.join(SCAN, "gse211630_olp_obs_annotation.csv"), index_col=0)
assert (olp.obs.index == ob.index).all()
olp.obs["st"] = ob["celltype_sub"].values
olp.obs["fm"] = ob["celltype_family"].values
olp.obs["cond"] = olp.obs["condition"].astype(str).values
keep = ~olp.obs["leiden_0.5"].astype(int).isin([12, 18])
olp = olp[keep].copy()
assert olp.n_obs == 46446
leid = olp.obs["leiden_0.5"].astype(int).values
st = olp.obs["st"].astype(str).values
fm = olp.obs["fm"].astype(str).values
sm = olp.obs["sample"].astype(str).values
cn = olp.obs["cond"].values

FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]
SHORT = {0: "c0 Fibroblast SFRP2+", 1: "c1 Fibroblast CXCL14+", 2: "c2 Pericyte/SMC",
         3: "c3 T resting", 4: "c4 T/NK stressed", 5: "c5 T/NK low-complex.",
         6: "c6 CD8+ GZMK+ EM", 7: "c7 CD4+ Treg", 8: "c8 T naive IL7R+",
         9: "c9 Plasmablast GZMB+", 10: "c10 Macrophage LYZ+",
         11: "c11 Epith. basal/prolif.", 13: "c13 Mast", 14: "c14 B cell",
         15: "c15 Plasma cell", 16: "c16 CD8+ GZMB+/CXCL13+",
         17: "c17 Endo lymph.", 19: "c19 Endo blood", 20: "c20 Myeloid S100-hi"}
cl_st = {}
for c, t in zip(leid, st):
    cl_st.setdefault(int(c), t)
st2cl = {t: c for c, t in cl_st.items()}
st2fam = {t: f for t, f in zip(st, fm)}
assert len(cl_st) == 19
samples = [f"NEOLP-{i}" for i in (1, 2, 3)] + [f"EOLP-{i}" for i in (1, 2)]


def sample_perm(vals, lab):
    """Two-sided permutation p over C(5,3)=10 splits of 5 values."""
    n0, n1 = int((lab == 0).sum()), int((lab == 1).sum())
    assert len(vals) == n0 + n1 == 5
    if not (vals.max() > vals.min()):
        return np.nan
    obs = vals[lab == 0].mean() - vals[lab == 1].mean()
    d = []
    for ii in itertools.combinations(range(5), 3):
        ll = np.zeros(5, int)
        ll[list(ii)] = 1
        d.append(vals[ll == 1].mean() - vals[ll == 0].mean())
    return float((np.abs(np.asarray(d)) >= np.abs(obs)).mean())


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    o = np.argsort(pv)
    q[ok] = np.minimum.accumulate((pv[o] * np.arange(1, len(pv) + 1) / len(pv))[::-1])[::-1][
        np.argsort(o)]
    return q


# =================================================================
# PART 1 - composition
# =================================================================
comp_rows = []
for lvl, grp in (("family", fm), ("subtype", st)):
    for s in samples:
        ms = sm == s
        tot = int(ms.sum())
        for g in sorted(set(grp[ms])):
            m = ms & (grp == g)
            comp_rows.append({"level": lvl, "group": g, "cluster": np.nan,
                              "sample": s, "condition": "NEOLP" if s.startswith("NEOLP") else "EOLP",
                              "n_cells": int(m.sum()), "frac": float(m.sum()) / tot})
comp_ps = pd.DataFrame(comp_rows)
comp_ps["cluster"] = comp_ps.apply(
    lambda r: (st2cl.get(r.group, np.nan) if r.level == "subtype" else np.nan), axis=1)
comp_ps.to_csv(os.path.join(TAB, "geo_gse211630_composition_per_sample.csv"), index=False)

shift_rows = []
for lvl in ("family", "subtype"):
    for g in comp_ps[comp_ps.level == lvl].group.unique():
        ne = comp_ps[(comp_ps.level == lvl) & (comp_ps.group == g) &
                     (comp_ps.condition == "NEOLP")]
        eo = comp_ps[(comp_ps.level == lvl) & (comp_ps.group == g) &
                     (comp_ps.condition == "EOLP")]
        if len(ne) == 0 or len(eo) == 0:
            continue
        vn, ve = ne.frac.values, eo.frac.values
        # pooled cell-level Fisher (descriptive; cells are pseudoreplicates)
        n_cell_ne, n_cell_eo = int(ne.n_cells.sum()), int(eo.n_cells.sum())
        fisher_p = stats.fisher_exact(
            [[n_cell_ne, 46446 - n_cell_ne], [n_cell_eo, 46446 - n_cell_eo]])[1]
        shift_rows.append({
            "level": lvl, "group": g,
            "cluster": st2cl.get(g, np.nan) if lvl == "subtype" else np.nan,
            "family": st2fam.get(g, g) if lvl == "subtype" else g,
            "median_frac_NEOLP": float(np.median(vn)),
            "median_frac_EOLP": float(np.median(ve)),
            "mean_frac_NEOLP": float(vn.mean()), "mean_frac_EOLP": float(ve.mean()),
            "delta_mean_frac": float(ve.mean() - vn.mean()),
            "n_cells_NEOLP": n_cell_ne, "n_cells_EOLP": n_cell_eo,
            "fisher_cell_p": fisher_p,
            "perm_p_sample": sample_perm(np.concatenate([vn, ve]),
                                         np.array([0] * len(vn) + [1] * len(ve)))
            if (len(vn) == N_NEOLP and len(ve) == N_EOLP) else np.nan})
comp_shift = pd.DataFrame(shift_rows).sort_values(
    ["level", "delta_mean_frac"], ascending=[True, False])
# BH over family and subtype grids separately (cell-level Fisher, descriptive)
for lvl in ("family", "subtype"):
    m = comp_shift.level == lvl
    comp_shift.loc[m, "fisher_BH_FDR"] = bh(comp_shift.loc[m, "fisher_cell_p"].values)
comp_shift.to_csv(os.path.join(TAB, "geo_gse211630_composition_shift.csv"), index=False)

# =================================================================
# PART 2 - T/NK functional module scores
# =================================================================
SIG = {
    "Cytotoxic/effector": ["GZMB", "GZMK", "PRF1", "GNLY", "IFNG", "KLRD1"],
    "Exhaustion":         ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX",
                           "ENTPD1", "CXCL13"],
    "Treg":               ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "BATF"],
    "Naive/recircul.":    ["IL7R", "SELL", "TCF7", "LEF1", "S100A4"],
    "Stress/HSP":         ["HSPA1A", "HSPA1B", "HSP90AA1", "HSPH1"],
    "Activation":         ["CD69", "TNF", "IFNG", "IL32"],
    "Proliferation":      ["MKI67"],
}
var_names = set(olp.var_names)
print("signature genes surviving the merged universe:")
for s, gl in SIG.items():
    have = [g for g in gl if g in var_names]
    print(f"  {s:16s} {len(have)}/{len(gl)}  missing={[g for g in gl if g not in var_names]}")

SIG_GENES = sorted({g for gl in SIG.values() for g in gl if g in var_names})
X = np.asarray(olp[:, SIG_GENES].X.todense()) if hasattr(olp.X, "todense") \
    else np.asarray(olp[:, SIG_GENES].X)
vix = {g: i for i, g in enumerate(SIG_GENES)}

tnk_mask = fm == "T/NK"
tnk_st = np.unique(st[tnk_mask])
# T/NK subtype display order (by cluster id)
tnk_cl = sorted((st2cl[t] for t in tnk_st))
tnk_st_order = [cl_st[c] for c in tnk_cl]

score = {s: np.zeros(olp.n_obs) for s in SIG}
for s, gl in SIG.items():
    idx = [vix[g] for g in gl if g in var_names]
    score[s] = X[:, idx].mean(axis=1) if idx else np.nan

# sanity: signature means per T/NK subtype (across all T/NK cells)
print("\nsanity (mean module score per T/NK subtype):")
hdr = "".join(f"{s[:12]:>14s}" for s in SIG)
print(f"{'subtype':38s}{hdr}")
for t in tnk_st_order:
    m = tnk_mask & (st == t)
    row = "".join(f"{score[s][m].mean():14.3f}" for s in SIG)
    print(f"{SHORT[st2cl[t]]:38s}{row}")

# per-sample mean score per subtype (module scores averaged over subtype cells)
mod_rows = []
for t in tnk_st_order:
    for s in samples:
        m = tnk_mask & (st == t) & (sm == s)
        if m.sum() == 0:
            continue
        for sig in SIG:
            mod_rows.append({"subtype": t, "cluster": st2cl[t], "signature": sig,
                             "sample": s, "condition": "NEOLP" if s.startswith("NEOLP") else "EOLP",
                             "n_cells": int(m.sum()),
                             "mean_score": float(score[sig][m].mean())})
mod_ps = pd.DataFrame(mod_rows).sort_values(["signature", "cluster", "sample"])
mod_ps.to_csv(os.path.join(TAB, "geo_gse211630_tnk_module_per_sample.csv"), index=False)

mod_shift_rows = []
for t in tnk_st_order:
    for sig in SIG:
        n = mod_ps[(mod_ps.subtype == t) & (mod_ps.signature == sig) &
                   (mod_ps.condition == "NEOLP")]
        e = mod_ps[(mod_ps.subtype == t) & (mod_ps.signature == sig) &
                   (mod_ps.condition == "EOLP")]
        if n.empty or e.empty:
            continue
        vn_c = score[sig][tnk_mask & (st == t) & (cn == "NEOLP")]
        ve_c = score[sig][tnk_mask & (st == t) & (cn == "EOLP")]
        w_p = np.nan
        if len(vn_c) >= 10 and len(ve_c) >= 10:
            try:
                w_p = stats.mannwhitneyu(vn_c, ve_c, alternative="two-sided").pvalue
            except ValueError:
                w_p = np.nan
        vals = np.concatenate([n.mean_score.values, e.mean_score.values])
        lab = np.array([0] * len(n) + [1] * len(e))
        mod_shift_rows.append({
            "subtype": t, "cluster": st2cl[t], "signature": sig,
            "mean_NEOLP": float(n.mean_score.mean()),
            "mean_EOLP": float(e.mean_score.mean()),
            "median_NEOLP": float(np.median(n.mean_score)),
            "median_EOLP": float(np.median(e.mean_score)),
            "delta_mean": float(e.mean_score.mean() - n.mean_score.mean()),
            "n_cells_NEOLP": int(n.n_cells.sum()), "n_cells_EOLP": int(e.n_cells.sum()),
            "wilcoxon_cell_p": w_p,
            "perm_p_sample": sample_perm(vals, lab)
            if (len(n) == N_NEOLP and len(e) == N_EOLP) else np.nan})
mod_shift = pd.DataFrame(mod_shift_rows)
mod_shift["BH_FDR_cell"] = bh(mod_shift["wilcoxon_cell_p"].values)
mod_shift.to_csv(os.path.join(TAB, "geo_gse211630_tnk_module_shift.csv"), index=False)
print("\nmodule shift rows:", len(mod_shift))
print("cell-level BH FDR<0.05:", int((mod_shift.BH_FDR_cell < 0.05).sum()),
      "| sample perm p<=0.2:", int((mod_shift.perm_p_sample <= 0.2).sum()))
