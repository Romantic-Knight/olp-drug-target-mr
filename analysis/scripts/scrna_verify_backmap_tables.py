# -*- coding: utf-8 -*-
"""
Programmatic verification of the Stage C (MR target reverse-mapping) tables:

  results/tables/geo_gse211630_mr_backmap_gene.csv      (5 rows)
  results/tables/geo_gse211630_mr_backmap_subtype.csv   (95 rows)

All numbers are re-derived independently from the verified Stage 4/A/B source
tables and checked against the mapping tables:
  * source tables: erosion_candidate_celltype (Stage A), composition_shift
    (Stage B), candidate_loc (Stage 4), final_target_evidence (MR final set)
  * N_NEOLP = 27,084 / N_EOLP = 19,362 pooled cell totals
  * MR direction: risk = OR>1 => expected tissue polarity +1, protective
    (OR<1) => -1

Writes results/reports/stageC_table_verification.txt and prints PASS count.
Run from project root (pandas+numpy only). exit code 0 on full PASS.
"""

import os
import sys
import numpy as np
import pandas as pd

ROOT = os.getcwd()
TAB = os.path.join(ROOT, "results", "tables")
LOG = os.path.join(ROOT, "results", "reports", "stageC_table_verification.txt")
N_NEOLP, N_EOLP = 27084, 19362

checks = []
_fails = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    if not cond:
        _fails.append(f"{name} :: {detail}")


def read(name):
    return pd.read_csv(os.path.join(TAB, name), encoding="utf-8-sig")


def main():
    g = read("geo_gse211630_mr_backmap_gene.csv")
    st = read("geo_gse211630_mr_backmap_subtype.csv")
    er = read("geo_gse211630_erosion_candidate_celltype.csv")
    comp = read("geo_gse211630_composition_shift.csv")
    loc = read("geo_gse211630_candidate_loc.csv")
    ev = read("final_target_evidence.csv")

    er_sub = er[er["level"] == "sub"].reset_index(drop=True)
    comp_sub = comp[comp["level"] == "subtype"].reset_index(drop=True)

    # ---- shape / schema ------------------------------------------------
    chk("gene table rows == 5", len(g) == 5, f"got {len(g)}")
    chk("subtype table rows == 95", len(st) == 95, f"got {len(st)}")
    req_g = {"gene", "tier", "mr_direction", "mr_or_wide", "mr_p_wide",
             "expected_polarity", "pool_expr_frac", "low_coverage",
             "tissue_mean_NEOLP", "tissue_mean_EOLP", "tissue_delta",
             "tissue_rel_pct", "observed_polarity", "consistency",
             "comp_contribution", "within_contribution",
             "interaction_contribution", "comp_share", "within_share",
             "interaction_share", "mechanism", "principal_host_subtype",
             "principal_host_family", "principal_host_expr_frac",
             "principal_host_delta_mean_all", "principal_host_FDR",
             "principal_host_perm_p", "principal_host_comp_delta_frac",
             "driver_subtype", "driver_total_contribution",
             "n_sig_within", "n_sig_within_consistent", "n_sig_within_opposite"}
    req_s = {"gene", "cluster", "celltype", "family", "n_cells_NEOLP",
             "n_cells_EOLP", "n_expr_NEOLP", "n_expr_EOLP", "frac_NEOLP",
             "frac_EOLP", "mean_all_NEOLP", "mean_all_EOLP",
             "delta_mean_all", "delta_frac", "wilcoxon_cell_p",
             "BH_FDR_cell", "perm_p_sample", "comp_delta_frac",
             "comp_perm_p_sample", "pool_frac_N", "pool_frac_E",
             "pool_frac_delta", "expr_frac_pooled", "comp_contribution",
             "within_contribution", "interaction_contribution",
             "total_contribution", "within_sign", "within_consistent",
             "sig_cell", "sig_perm", "sig_comp_perm"}
    chk("gene cols present", req_g <= set(g.columns),
        str(req_g - set(g.columns)))
    chk("subtype cols present", req_s <= set(st.columns),
        str(req_s - set(st.columns)))

    # ---- MR metadata vs final evidence --------------------------------
    order = ev["gene"].tolist()
    chk("gene order == MR evidence order",
        g["gene"].tolist() == order, str(g["gene"].tolist()))
    for _, r in g.iterrows():
        e = ev[ev["gene"] == r["gene"]].iloc[0]
        chk(f"{r['gene']} or_wide matches evidence",
            np.isclose(r["mr_or_wide"], e["or_wide"]))
        chk(f"{r['gene']} tier matches evidence", r["tier"] == e["tier"])
    dir_map = {"risk": 1, "protective": -1}
    chk("mr_direction derived from OR sign",
        all(g.loc[i, "mr_direction"] ==
            ("risk" if g.loc[i, "mr_or_wide"] > 1 else "protective")
            for i in g.index))
    chk("expected_polarity = +1 risk / -1 protective",
        all(g.loc[i, "expected_polarity"] == dir_map[g.loc[i, "mr_direction"]]
            for i in g.index))

    # ---- subtype table faithful copy of sources ------------------------
    m = st.merge(
        er_sub[["gene", "cluster", "n_cells_NEOLP", "n_cells_EOLP",
                "n_expr_NEOLP", "n_expr_EOLP", "mean_all_NEOLP",
                "mean_all_EOLP", "delta_mean_all", "BH_FDR_cell",
                "perm_p_sample"]],
        on=["gene", "cluster"], suffixes=("", "_src"), validate="one_to_one")
    for col in ["n_cells_NEOLP", "n_cells_EOLP", "n_expr_NEOLP",
                "n_expr_EOLP", "mean_all_NEOLP", "mean_all_EOLP",
                "delta_mean_all", "BH_FDR_cell", "perm_p_sample"]:
        chk(f"subtype {col} == Stage A source",
            np.allclose(m[col], m[col + "_src"], rtol=1e-9, atol=1e-12,
                        equal_nan=True))
    # Stage B composition fields aligned by cluster (direct map comparison)
    comp_map = comp_sub.set_index("cluster")[
        ["delta_mean_frac", "perm_p_sample", "family"]]
    chk("subtype family == Stage B family",
        (st["family"] == st["cluster"].map(comp_map["family"])).all())
    chk("subtype comp_delta_frac == Stage B delta_mean_frac",
        np.allclose(st["comp_delta_frac"],
                    st["cluster"].map(comp_map["delta_mean_frac"]),
                    rtol=1e-9, atol=1e-12))
    chk("subtype comp_perm_p_sample == Stage B perm_p_sample",
        np.allclose(st["comp_perm_p_sample"],
                    st["cluster"].map(comp_map["perm_p_sample"]),
                    rtol=1e-9, atol=1e-12, equal_nan=True))

    # ---- pooled totals & pooling columns ------------------------------
    tot = st.groupby("gene")[["n_cells_NEOLP", "n_cells_EOLP"]].sum()
    chk("all genes pooled NEOLP total == 27084",
        np.all(tot["n_cells_NEOLP"] == N_NEOLP),
        str(tot["n_cells_NEOLP"].tolist()))
    chk("all genes pooled EOLP total == 19362",
        np.all(tot["n_cells_EOLP"] == N_EOLP),
        str(tot["n_cells_EOLP"].tolist()))
    chk("pool_frac_N == n_cells_N/27084",
        np.allclose(st["pool_frac_N"], st["n_cells_NEOLP"] / N_NEOLP, atol=1e-12))
    chk("pool_frac_E == n_cells_E/19362",
        np.allclose(st["pool_frac_E"], st["n_cells_EOLP"] / N_EOLP, atol=1e-12))
    chk("pool_frac_delta == E - N",
        np.allclose(st["pool_frac_delta"],
                    st["pool_frac_E"] - st["pool_frac_N"], atol=1e-12))
    chk("expr_frac_pooled formula",
        np.allclose(st["expr_frac_pooled"],
                    (st["n_expr_NEOLP"] + st["n_expr_EOLP"]) /
                    (st["n_cells_NEOLP"] + st["n_cells_EOLP"]),
                    atol=1e-12))
    chk("no subtype rows without a source row (19 per gene)",
        st.groupby("gene").size().eq(19).all())

    # ---- derived per-subtype columns ----------------------------------
    chk("comp_contribution formula",
        np.allclose(st["comp_contribution"],
                    st["pool_frac_delta"] * st["mean_all_NEOLP"], atol=1e-12))
    chk("within_contribution formula",
        np.allclose(st["within_contribution"],
                    st["pool_frac_N"] *
                    (st["mean_all_EOLP"] - st["mean_all_NEOLP"]), atol=1e-12))
    chk("interaction_contribution formula",
        np.allclose(st["interaction_contribution"],
                    st["pool_frac_delta"] *
                    (st["mean_all_EOLP"] - st["mean_all_NEOLP"]), atol=1e-12))
    chk("total_contribution additive",
        np.allclose(st["total_contribution"],
                    st["comp_contribution"] + st["within_contribution"] +
                    st["interaction_contribution"], atol=1e-12))
    exp_sign = np.sign(st["mean_all_EOLP"] - st["mean_all_NEOLP"]).replace(0, np.nan)
    chk("within_sign == sign(delta_mean_all)",
        st["within_sign"].equals(exp_sign) or
        np.allclose(st["within_sign"], exp_sign, rtol=0, atol=0, equal_nan=True))
    chk("within_consistent logic",
        (st["within_consistent"] == (st["within_sign"] ==
                                     st["expected_polarity"])).all())
    chk("sig_cell logic", st["sig_cell"].equals(st["BH_FDR_cell"] < 0.05))
    chk("sig_perm logic", st["sig_perm"].equals(st["perm_p_sample"] <= 0.2))
    chk("sig_comp_perm logic",
        st["sig_comp_perm"].equals(st["comp_perm_p_sample"] <= 0.2))
    chk("expected_polarity propagated per gene",
        all(st.loc[st.gene == gg, "expected_polarity"].eq(pol).all()
            for gg, pol in
            zip(g["gene"], g["expected_polarity"])))

    # NaN semantics: perm_p_sample & comp_perm NaN appear only where the
    # source tables had NaN (sparse combinations / all ties)
    # NaN semantics: perm_p_sample NaNs match Stage A source (sparse/all-tie)
    for gg, dd in st.groupby("gene"):
        src = er_sub[er_sub["gene"] == gg]
        chk(f"{gg} perm NaN semantics",
            dd["perm_p_sample"].isna().sum() == src["perm_p_sample"].isna().sum())

    # ---- gene-level: reconstruction + decomposition + labels -----------
    for _, r in g.iterrows():
        dd = st[st["gene"] == r["gene"]]
        pN, pE = dd["pool_frac_N"].to_numpy(), dd["pool_frac_E"].to_numpy()
        mN = dd["mean_all_NEOLP"].to_numpy()
        mE = dd["mean_all_EOLP"].to_numpy()
        MN = float((pN * mN).sum())
        ME = float((pE * mE).sum())
        dM = ME - MN
        chk(f"{r['gene']} tissue_mean_NEOLP reconstruction",
            np.isclose(MN, r["tissue_mean_NEOLP"], atol=1e-9),
            f"{MN} vs {r['tissue_mean_NEOLP']}")
        chk(f"{r['gene']} tissue_mean_EOLP reconstruction",
            np.isclose(ME, r["tissue_mean_EOLP"], atol=1e-9))
        chk(f"{r['gene']} tissue_delta additive",
            np.isclose(dM, r["tissue_delta"], atol=1e-9))
        comp_t = r["comp_contribution"]
        within_t = r["within_contribution"]
        inter_t = r["interaction_contribution"]
        chk(f"{r['gene']} decomposition sums to delta",
            np.isclose(comp_t + within_t + inter_t, dM, atol=1e-9))
        rel = (dM / MN * 100.0) if MN else np.nan
        chk(f"{r['gene']} tissue_rel_pct formula",
            np.isclose(rel, r["tissue_rel_pct"], atol=1e-6))
        # whole-object loc cross-check (Stage 4 independent source)
        lr = loc[loc["gene"] == r["gene"]].iloc[0]
        chk(f"{r['gene']} matches candidate_loc whole-object NEOLP",
            np.isclose(MN, lr["NEOLP_mean_log1p"], atol=5e-3),
            f"{MN} vs {lr['NEOLP_mean_log1p']}")
        chk(f"{r['gene']} matches candidate_loc whole-object EOLP",
            np.isclose(ME, lr["EOLP_mean_log1p"], atol=5e-3))
        chk(f"{r['gene']} pool_expr_frac == candidate_loc frac_cells",
            np.isclose(r["pool_expr_frac"], lr["frac_cells"], atol=1e-9))
        chk(f"{r['gene']} consistency label",
            (r["consistency"] == "consistent") ==
            (dM != 0 and np.sign(dM) == r["expected_polarity"]),
            f"{r['consistency']} dM={dM} exp={r['expected_polarity']}")
        chk(f"{r['gene']} mechanism label",
            r["mechanism"] == ("composition_dominant"
                               if abs(comp_t) >= abs(within_t)
                               else "within_dominant"))
        host = dd.loc[dd["expr_frac_pooled"].idxmax()]
        chk(f"{r['gene']} principal host argmax",
            host["celltype"] == r["principal_host_subtype"],
            f"{host['celltype']} vs {r['principal_host_subtype']}")
        chk(f"{r['gene']} principal host frac",
            np.isclose(host["expr_frac_pooled"], r["principal_host_expr_frac"],
                       atol=1e-9))
        drv = dd.loc[dd["total_contribution"].abs().idxmax()]
        chk(f"{r['gene']} driver subtype argmax|total_contribution|",
            drv["celltype"] == r["driver_subtype"],
            f"{drv['celltype']} vs {r['driver_subtype']}")
        chk(f"{r['gene']} driver contribution value",
            np.isclose(drv["total_contribution"],
                       r["driver_total_contribution"], atol=1e-9))
        mask = dd["sig_cell"] | dd["sig_perm"]
        chk(f"{r['gene']} n_sig_within", mask.sum() == r["n_sig_within"])
        cons = dd.loc[mask & dd["within_consistent"]].shape[0]
        chk(f"{r['gene']} n_sig_within_consistent", cons == r["n_sig_within_consistent"])
        opp = dd.loc[mask & ~dd["within_consistent"] &
                     dd["within_sign"].notna()].shape[0]
        chk(f"{r['gene']} n_sig_within_opposite", opp == r["n_sig_within_opposite"])
        # low coverage flag vs pool expr fraction
        chk(f"{r['gene']} low_coverage flag", r["low_coverage"] ==
            bool(r["pool_expr_frac"] < 0.05))

    # ---- duplicates / integrity ----------------------------------------
    chk("no duplicated gene rows", g["gene"].duplicated().sum() == 0)
    chk("no duplicated subtype rows",
        st.duplicated(subset=["gene", "cluster"]).sum() == 0)
    chk("subtype cluster values == expected 19 clusters",
        set(st["cluster"]) == {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13,
                               14, 15, 16, 17, 19, 20})
    chk("no unexpected infinite values in float tables",
        not np.isinf(g.select_dtypes("number").to_numpy()).any() and
        not np.isinf(st.select_dtypes("number").to_numpy()).any())

    # ---- write log ------------------------------------------------------
    n_pass = sum(ok for _, ok in checks)
    n_fail = sum(not ok for _, ok in checks)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "w", encoding="utf-8") as f:
        f.write("Stage C table verification (geo_gse211630_mr_backmap_*)\n")
        f.write(f"checks: {len(checks)}  pass: {n_pass}  fail: {n_fail}\n")
        f.write("=" * 78 + "\n")
        for name, ok in checks:
            f.write(f"[{'PASS' if ok else 'FAIL'}] {name}\n")
        if _fails:
            f.write("-" * 78 + "\nFAIL details:\n")
            f.write("\n".join(_fails) + "\n")
    print(f"Stage C: {n_pass}/{len(checks)} checks passed, "
          f"{n_fail} failed -> {LOG}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
