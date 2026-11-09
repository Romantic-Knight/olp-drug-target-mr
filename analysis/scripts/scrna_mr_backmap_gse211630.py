# -*- coding: utf-8 -*-
"""
Stage C (MR target reverse-mapping, GSE211630): integrate the five finalized
OLP MR drug-target genes (FinnGen R13) with the scRNA-seq subtype-level
expression and composition tables produced in Stage A (erosion) and Stage B
(composition/T-NK).

Scientific question
-------------------
For each MR candidate (risk / protective based on cis-pQTL MR direction), does
the EOLP lesion tissue show the polarity predicted by the MR direction?
Two layers are inspected:
  1) tissue-level (whole 46,446-cell OLP universe) change in mean log1p
     expression, EOLP vs NEOLP, and its decomposition into a *composition*
     term (subtype fraction shifts) vs a *within-celltype* term (per-cell
     expression changes inside each subtype) vs the interaction term
     (Oaxaca-Blinder on pooled subtype fractions);
  2) subtype-level mapping: pooled expression fraction (localization) and the
     within-subtype direction of change, referenced against each subtype's
     Stage-B composition shift.

Convention (declared in README, Stage C section)
------------------------------------------------
expected tissue polarity = +1 for risk genes (OR>1: higher circulating protein
raises OLP risk -> higher tissue message expected in EOLP), -1 for protective
genes (OR<1 -> lower tissue message expected in EOLP). "consistent" means the
observed tissue delta sign equals the expected polarity; "opposite" otherwise.
Caveats: mRNA abundance, not protein; whole-tissue mixture, not circulating
pQTL compartment; no tissue pQTL is claimed.

Inputs (all previously verified, Stage 4/A/B):
  results/tables/geo_gse211630_erosion_candidate_celltype.csv   (Stage A)
  results/tables/geo_gse211630_composition_shift.csv            (Stage B)
  results/tables/geo_gse211630_candidate_loc.csv                (Stage 4)
  results/tables/final_target_evidence.csv                      (MR final)

Outputs:
  results/tables/geo_gse211630_mr_backmap_gene.csv       (5 rows, verdicts)
  results/tables/geo_gse211630_mr_backmap_subtype.csv    (5 x 19 = 95 rows)

Run: from project root with a python that has pandas+numpy (no scanpy/pixi
needed). Reproduced in the persistent wisp python kernel.
"""

import os
import numpy as np
import pandas as pd

ROOT = os.getcwd()  # project root (kernel / CLI started at project root)
TBL = os.path.join(ROOT, "results", "tables")

N_NEOLP, N_EOLP = 27084, 19362  # Stage B universe (46,446 annotated cells)


def read(name):
    p = os.path.join(TBL, name)
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    return pd.read_csv(p, encoding="utf-8-sig")


def main():
    er = read("geo_gse211630_erosion_candidate_celltype.csv")
    comp = read("geo_gse211630_composition_shift.csv")
    loc = read("geo_gse211630_candidate_loc.csv")
    ev = read("final_target_evidence.csv")

    er = er[er["level"] == "sub"].copy()
    comp_sub = comp[comp["level"] == "subtype"].copy()

    # ---- assertion 1: pooled cell totals per gene match the Stage B universe
    tot = er.groupby("gene")[["n_cells_NEOLP", "n_cells_EOLP"]].sum()
    assert np.all(tot["n_cells_NEOLP"] == N_NEOLP), "NEOLP pooled total mismatch"
    assert np.all(tot["n_cells_EOLP"] == N_EOLP), "EOLP pooled total mismatch"

    # ---- MR direction from final evidence (OR of circulating protein)
    mr = ev[["gene", "tier", "or_wide", "p_wide", "ci_wide"]].copy()
    mr["mr_direction"] = np.where(mr["or_wide"] > 1, "risk", "protective")
    mr["expected_polarity"] = np.where(mr["mr_direction"] == "risk", 1, -1)

    # ---- subtype-level table: merge Stage A erosion + Stage B composition
    st = er.merge(
        comp_sub[
            [
                "cluster", "group",
                "delta_mean_frac", "perm_p_sample",
            ]
        ].rename(
            columns={
                "group": "celltype_sub_comp",
                "delta_mean_frac": "comp_delta_frac",
                "perm_p_sample": "comp_perm_p_sample",
            }
        ),
        on="cluster", how="left", validate="many_to_one",
    )
    # sanity: subtype labels align between Stage A and Stage B
    assert (st["celltype"] == st["celltype_sub_comp"]).all(), \
        "Stage A / Stage B subtype labels mismatch"
    st.drop(columns=["celltype_sub_comp"], inplace=True)

    st["pool_frac_N"] = st["n_cells_NEOLP"] / N_NEOLP
    st["pool_frac_E"] = st["n_cells_EOLP"] / N_EOLP
    st["pool_frac_delta"] = st["pool_frac_E"] - st["pool_frac_N"]
    st["expr_frac_pooled"] = (st["n_expr_NEOLP"] + st["n_expr_EOLP"]) / (
        st["n_cells_NEOLP"] + st["n_cells_EOLP"]
    )

    # ---- Oaxaca-Blinder decomposition (base = NEOLP)
    st["comp_contribution"] = st["pool_frac_delta"] * st["mean_all_NEOLP"]
    st["within_contribution"] = st["pool_frac_N"] * (
        st["mean_all_EOLP"] - st["mean_all_NEOLP"]
    )
    st["interaction_contribution"] = st["pool_frac_delta"] * (
        st["mean_all_EOLP"] - st["mean_all_NEOLP"]
    )
    st["total_contribution"] = (
        st["comp_contribution"]
        + st["within_contribution"]
        + st["interaction_contribution"]
    )
    st = st.merge(mr[["gene", "mr_direction", "expected_polarity"]], on="gene")

    # subtype-level within-consistency vs MR per-cell expectation
    d = st["mean_all_EOLP"] - st["mean_all_NEOLP"]
    st["within_sign"] = np.sign(d).replace(0, np.nan)
    st["within_consistent"] = st["within_sign"] == st["expected_polarity"]
    st["sig_cell"] = st["BH_FDR_cell"] < 0.05
    st["sig_perm"] = st["perm_p_sample"] <= 0.2
    st["sig_comp_perm"] = st["comp_perm_p_sample"] <= 0.2

    # ---- gene-level table
    genes = []
    for g, dd in st.groupby("gene"):
        pN, pE = dd["pool_frac_N"].to_numpy(), dd["pool_frac_E"].to_numpy()
        mN = dd["mean_all_NEOLP"].to_numpy()
        mE = dd["mean_all_EOLP"].to_numpy()
        M_N, M_E = float((pN * mN).sum()), float((pE * mE).sum())
        dM = M_E - M_N
        comp_t = float(((pE - pN) * mN).sum())
        within_t = float((pN * (mE - mN)).sum())
        inter_t = float(((pE - pN) * (mE - mN)).sum())
        rel = (dM / M_N * 100.0) if M_N else np.nan

        # principal host = subtype with highest pooled expression fraction
        host_row = dd.loc[dd["expr_frac_pooled"].idxmax()]
        host_frac = host_row["expr_frac_pooled"]
        # Stage-B composition shift of that host subtype
        host_comp_delta = host_row["comp_delta_frac"]

        # principal driver = subtype with the largest |total contribution|
        drv_row = dd.loc[dd["total_contribution"].abs().idxmax()]

        # counts among statistically-supported within-subtype comparisons
        sigmask = dd["sig_cell"] | dd["sig_perm"]
        sig_consistent = int(
            dd.loc[sigmask & dd["within_consistent"]].shape[0]
        )
        sig_opposite = int(
            dd.loc[sigmask & ~dd["within_consistent"] & dd["within_sign"].notna()].shape[0]
        )

        ev_row = mr[mr["gene"] == g].iloc[0]
        loc_row = loc[loc["gene"] == g].iloc[0]

        # The dominant term is chosen across all three components. Comparing only
        # |composition| with |within-cell-type| ignored the interaction term, which is
        # the largest of the three for EGFL8, and mislabelled that gene.
        _terms = {"composition_dominant": comp_t,
                  "within_dominant": within_t,
                  "interaction_dominant": inter_t}
        mech = max(_terms, key=lambda k: abs(_terms[k]))
        pol = "up" if dM > 0 else ("down" if dM < 0 else "flat")
        consistency = (
            "consistent"
            if (dM != 0 and np.sign(dM) == ev_row["expected_polarity"])
            else "opposite"
        )
        genes.append(
            dict(
                gene=g,
                tier=ev_row["tier"],
                mr_direction=ev_row["mr_direction"],
                mr_or_wide=ev_row["or_wide"],
                mr_ci_wide=ev_row["ci_wide"],
                mr_p_wide=ev_row["p_wide"],
                mr_effect_wide=f"OR {ev_row['or_wide']:.2f} "
                               f"({ev_row['ci_wide']})",
                expected_polarity=ev_row["expected_polarity"],
                pool_expr_frac=loc_row["frac_cells"],
                low_coverage=bool(loc_row["frac_cells"] < 0.05),
                tissue_mean_NEOLP=M_N,
                tissue_mean_EOLP=M_E,
                tissue_delta=dM,
                tissue_rel_pct=rel,
                observed_polarity=pol,
                consistency=consistency,
                comp_contribution=comp_t,
                within_contribution=within_t,
                interaction_contribution=inter_t,
                comp_share=(comp_t / dM) if dM else np.nan,
                within_share=(within_t / dM) if dM else np.nan,
                interaction_share=(inter_t / dM) if dM else np.nan,
                mechanism=mech,
                principal_host_subtype=host_row["celltype"],
                principal_host_family=host_row["family"],
                principal_host_expr_frac=host_frac,
                principal_host_delta_mean_all=host_row["delta_mean_all"],
                principal_host_FDR=host_row["BH_FDR_cell"],
                principal_host_perm_p=host_row["perm_p_sample"],
                principal_host_comp_delta_frac=host_comp_delta,
                driver_subtype=drv_row["celltype"],
                driver_total_contribution=drv_row["total_contribution"],
                n_sig_within=sigmask.sum(),
                n_sig_within_consistent=sig_consistent,
                n_sig_within_opposite=sig_opposite,
            )
        )
    gene_df = pd.DataFrame(genes)

    # stable ordering: MR tier (final evidence order) then cluster
    order = ev["gene"].tolist()
    gene_df["_ord"] = gene_df["gene"].map({g: i for i, g in enumerate(order)})
    gene_df = gene_df.sort_values("_ord").drop(columns="_ord").reset_index(drop=True)
    st = st.sort_values(
        ["gene", "cluster"],
        key=lambda s: s.map({g: i for i, g in enumerate(order)})
        if s.name == "gene" else s,
    ).reset_index(drop=True)

    gene_df.to_csv(
        os.path.join(TBL, "geo_gse211630_mr_backmap_gene.csv"),
        index=False, encoding="utf-8-sig",
    )
    st.to_csv(
        os.path.join(TBL, "geo_gse211630_mr_backmap_subtype.csv"),
        index=False, encoding="utf-8-sig",
    )

    # ---- console summary
    cols = [
        "gene", "tier", "mr_direction", "mr_effect_wide", "pool_expr_frac",
        "tissue_mean_NEOLP", "tissue_mean_EOLP", "tissue_delta",
        "tissue_rel_pct", "observed_polarity", "consistency",
        "mechanism", "comp_share", "within_share", "interaction_share",
        "driver_subtype", "driver_total_contribution",
        "principal_host_subtype", "principal_host_expr_frac",
        "principal_host_comp_delta_frac", "n_sig_within_consistent",
        "n_sig_within_opposite",
    ]
    print(gene_df[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\nsubtype table rows: {st.shape[0]}  (expect 95)")
    print("written: geo_gse211630_mr_backmap_gene.csv / _subtype.csv")


if __name__ == "__main__":
    main()
