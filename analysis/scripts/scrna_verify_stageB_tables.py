# -*- coding: utf-8 -*-
"""Stage B table verification (GSE211630 OLP scRNA-seq).

Re-derives the Stage-B statistics from the per-sample tables with an
independent reimplementation of the two tracks and compares against the
already-written shift tables:

  1. composition_per_sample : row inventory + family fractions sum to 1
     + cell-count totals; absent (cluster x sample) combinations must NOT be
     imputed (c20 / NEOLP-2 has no row).
  2. composition_shift      : 8 family + 19 subtype rows; perm_p_sample
     semantics (NaN only when a sample-combination is absent or all 5
     fractions are tied, else multiple of 0.1, min 0.1); fisher_BH_FDR is
     BH(fisher_cell_p) over all 27 rows and monotone in p.
  3. tnk_module_per_sample  : 7 subtypes x 7 signatures x 5 samples = 245
     rows; sparse (n<10) accounting.
  4. tnk_module_shift       : 49 rows; BH_FDR_cell = BH(wilcoxon_cell_p)
     monotone; FDR<0.05 and perm p<=0.2 hit counts; perm_p_sample
     re-derived from per-sample means matches the table.

Usage: run with the project pixi python from analysis/scrna.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
TAB = os.path.join(ROOT, "results", "tables")
SAMPLES = [f"NEOLP-{i}" for i in (1, 2, 3)] + [f"EOLP-{i}" for i in (1, 2)]
COND_ORDER = ["NEOLP", "EOLP"]
FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]
SHORT_KEYS = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20}


def bh_orig(p):
    """Exact BH used by scrna_composition_tnk_gse211630.py (no clip at 1)."""
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    o = np.argsort(pv)
    q[ok] = np.minimum.accumulate(
        (pv[o] * np.arange(1, len(pv) + 1) / len(pv))[::-1])[::-1][np.argsort(o)]
    return q


def sample_perm5(vals):
    """Two-sided permutation p over the C(5,3)=10 splits; NaN if no spread."""
    vals = np.asarray(vals, float)
    assert vals.size == 5
    if not (vals.max() > vals.min()):
        return np.nan
    obs = vals[:3].mean() - vals[3:].mean()       # labels: NEOLP(3) vs EOLP(2)
    d = []
    for ii in itertools.combinations(range(5), 3):
        ll = np.zeros(5, int)
        ll[list(ii)] = 1
        d.append(vals[ll == 1].mean() - vals[ll == 0].mean())
    return float((np.abs(np.asarray(d)) >= np.abs(obs)).mean())


def close(a, b, tol=1e-9):
    a, b = float(a), float(b)
    if np.isnan(a) and np.isnan(b):
        return True
    return abs(a - b) <= tol


def main():
    fails = []

    def check(name, cond, extra=""):
        flag = "PASS" if cond else "FAIL"
        print(f"[{flag}] {name} {extra}")
        if not cond:
            fails.append(name)

    # ------------------------------------------------------------- 1. per-sample composition
    cps = pd.read_csv(os.path.join(TAB, "geo_gse211630_composition_per_sample.csv"))
    check("comp per_sample: 40 family rows (8 fam x 5 samples)", cps.level.eq("family").sum() == 40)
    check("comp per_sample: 94 subtype rows (19 - 1 missing c20/NEOLP-2)",
          cps.level.eq("subtype").sum() == 94)
    check("comp per_sample: 134 total rows (missing combo NOT imputed)",
          len(cps) == 134)
    fam = cps[cps.level.eq("family")]
    tot = fam.groupby("sample")["frac"].sum()
    check("comp per_sample: family frac sums to 1 per sample",
          np.allclose(tot.reindex(SAMPLES).values, 1.0, atol=1e-9))
    ntot = fam.groupby("sample")["n_cells"].sum().reindex(SAMPLES)
    check("comp per_sample: cell totals 46,446 (27,084 NEOLP + 19,362 EOLP)",
          int(ntot.sum()) == 46446 and int(ntot[["NEOLP-1", "NEOLP-2", "NEOLP-3"]].sum()) == 27084
          and int(ntot[["EOLP-1", "EOLP-2"]].sum()) == 19362,
          f"(got NEOLP={int(ntot[['NEOLP-1','NEOLP-2','NEOLP-3']].sum())}, "
          f"EOLP={int(ntot[['EOLP-1','EOLP-2']].sum())})")
    check("comp per_sample: every listed combo has n_cells>0 and frac>0",
          (cps.n_cells > 0).all() and (cps.frac > 0).all())
    check("comp per_sample: no duplicate (level,group,sample)",
          not cps.duplicated(["level", "group", "sample"]).any())
    sub_ps = cps[cps.level.eq("subtype")]
    check("comp per_sample: subtype cluster set == annotated 19 clusters",
          set(sub_ps.cluster.astype(int).unique()) == SHORT_KEYS)

    # ------------------------------------------------------------- 2. composition shift
    csh = pd.read_csv(os.path.join(TAB, "geo_gse211630_composition_shift.csv"))
    check("comp shift: 27 rows (8 family + 19 subtype)", len(csh) == 27)
    check("comp shift: family rows have NaN cluster", csh[csh.level.eq("family")].cluster.isna().all())
    check("comp shift: subtype rows cluster == annotated set",
          set(csh[csh.level.eq("subtype")].cluster.astype(int)) == SHORT_KEYS)
    # BH over family and subtype grids SEPARATELY (as in the analysis script)
    fp = csh.fisher_cell_p.values
    fb = csh.fisher_BH_FDR.values
    q_re = np.full(len(csh), np.nan)
    for lvl in ("family", "subtype"):
        m = (csh.level == lvl).values
        q_re[m] = bh_orig(fp[m])
    o = np.argsort(fp, kind="stable")
    check("comp shift: fisher_BH_FDR is monotone in fisher_cell_p",
          np.all(np.diff(fb[o]) >= -1e-12))
    check("comp shift: fisher_BH_FDR == BH(fisher_cell_p) per level (family|subtype)",
          np.allclose(fb, q_re, atol=1e-12, equal_nan=True))
    # perm p re-derived per group from per-sample fractions
    perm_ok = True
    perm_reason = {}
    for _, r in csh.iterrows():
        key = ("subtype" if r.level == "subtype" else "family", r.group)
        rows = cps[(cps.level == r.level) & (cps.group == r.group)]
        vals = rows.set_index("sample")["frac"].reindex(SAMPLES)
        missing = vals.isna().sum()
        if missing > 0:
            expect = np.nan
            perm_reason[key] = f"missing {missing} sample-combination(s)"
        else:
            expect = sample_perm5(vals.values)
            perm_reason[key] = ""
        if not close(r.perm_p_sample, expect):
            perm_ok = False
            perm_reason[key] = (f"table={r.perm_p_sample} vs re-derived={expect} " + perm_reason[key])
    check("comp shift: perm_p_sample == re-derived permutation (NaN = absent combo or tie)",
          perm_ok)
    pv = csh.perm_p_sample.dropna().values
    check("comp shift: non-NaN perm p are multiples of 0.1 with min 0.1",
          pv.size > 0 and np.allclose(pv, np.round(pv, 6)) and
          np.allclose((pv * 10) % 1, 0, atol=1e-6) and pv.min() >= 0.1 - 1e-9,
          f"(n={pv.size}, min={pv.min():.1f})")
    nan_rows = csh[csh.perm_p_sample.isna()]
    check("comp shift: NaN perm rows = single c20 (cluster 20) subtype row",
          len(nan_rows) == 1 and nan_rows.iloc[0].level == "subtype"
          and int(nan_rows.iloc[0].cluster) == 20,
          f"(got {list(map(tuple, nan_rows[['level','cluster']].values))})")

    # ------------------------------------------------------------- 3. per-sample module scores
    mps = pd.read_csv(os.path.join(TAB, "geo_gse211630_tnk_module_per_sample.csv"))
    check("module per_sample: 245 rows (7 sub x 7 sig x 5 samples)", len(mps) == 245)
    check("module per_sample: 7 subtypes, 7 signatures, 5 samples each",
          mps.groupby(["subtype", "signature"]).size().eq(5).all() and
          mps.signature.nunique() == 7 and mps.subtype.nunique() == 7)
    check("module per_sample: no NaN mean_score, n_cells>0",
          mps.mean_score.notna().all() and (mps.n_cells > 0).all())
    sparse = mps[mps.n_cells < 10]
    check("module per_sample: sparse (<10 cells) rows flagged and present",
          len(sparse) > 0, f"(n_sparse={len(sparse)})")

    # ------------------------------------------------------------- 4. module shift
    msh = pd.read_csv(os.path.join(TAB, "geo_gse211630_tnk_module_shift.csv"))
    check("module shift: 49 rows (7 sub x 7 sig)", len(msh) == 49)
    check("module shift: unique (subtype,signature)",
          msh.duplicated(["subtype", "signature"]).sum() == 0)
    wp = msh.wilcoxon_cell_p.values
    wq = msh.BH_FDR_cell.values
    o = np.argsort(wp, kind="stable")
    check("module shift: BH_FDR_cell monotone in wilcoxon_cell_p",
          np.all(np.diff(wq[o]) >= -1e-12))
    q_re = bh_orig(wp)
    check("module shift: BH_FDR_cell == BH(wilcoxon_cell_p) over all 49 rows",
          np.allclose(wq, q_re, atol=1e-12, equal_nan=True))
    n_fdr = int((wq < 0.05).sum())
    check("module shift: cell-level BH FDR<0.05 count == 33", n_fdr == 33, f"(got {n_fdr})")
    pp = msh.perm_p_sample
    n_perm = int((pp <= 0.2).sum())
    check("module shift: sample-perm p<=0.2 count == 16 (incl. NaN-aware)", n_perm == 16,
          f"(got {n_perm})")
    check("module shift: perm p non-NaN multiples of 0.1, min 0.1",
          pp.dropna().min() >= 0.1 - 1e-9 and
          np.allclose((pp.dropna().values * 10) % 1, 0, atol=1e-6))
    # re-derive per-sample permutation from per-sample means + group means
    ok_perm = True
    for _, r in msh.iterrows():
        rows = mps[(mps.subtype == r.subtype) & (mps.signature == r.signature)]
        vals = rows.set_index("sample")["mean_score"].reindex(SAMPLES)
        if vals.isna().any():
            expect = np.nan
        else:
            expect = sample_perm5(vals.values)
        if not close(r.perm_p_sample, expect):
            ok_perm = False
            print(f"  perm mismatch {r.subtype} x {r.signature}: "
                  f"table={r.perm_p_sample} re={expect}")
    check("module shift: perm_p_sample == re-derived from per-sample means", ok_perm)
    # group mean columns = mean over per-sample means (unweighted)
    ok_mean = True
    for _, r in msh.iterrows():
        rows = mps[(mps.subtype == r.subtype) & (mps.signature == r.signature)]
        for cond in COND_ORDER:
            m = rows[rows.condition == cond].mean_score.mean()
            col = f"mean_{cond}"
            if not np.isclose(r[col], m, atol=1e-12):
                ok_mean = False
                print(f"  mean mismatch {r.subtype} x {r.signature} {cond}: "
                      f"table={r[col]} vs re={m}")
    check("module shift: mean_NEOLP/EOLP == unweighted mean of sample means", ok_mean)
    # sample n totals consistent with per-sample n_cells
    ok_n = True
    for _, r in msh.iterrows():
        rows = mps[(mps.subtype == r.subtype) & (mps.signature == r.signature)]
        for cond in COND_ORDER:
            nsum = int(rows[rows.condition == cond].n_cells.sum())
            col = f"n_cells_{cond}"
            if r[col] != nsum:
                ok_n = False
                print(f"  n mismatch {r.subtype} x {r.signature} {cond}: "
                      f"table={r[col]} vs sum={nsum}")
    check("module shift: n_cells_NEOLP/EOLP == sum of per-sample n_cells", ok_n)

    # ------------------------------------------------------------- summary
    print()
    if fails:
        print(f"VERIFICATION: {len(fails)} FAILED -> {fails}")
        sys.exit(1)
    print("VERIFICATION: all Stage-B table checks PASS")


if __name__ == "__main__":
    main()
