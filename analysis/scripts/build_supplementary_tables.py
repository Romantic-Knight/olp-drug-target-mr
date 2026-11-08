"""Build the curated supplementary tables for the OLP drug-target MR manuscript.

Inputs: results/tables/*.csv (existing analysis outputs)
Outputs:
  results/tables/supplementary/suppl_table_01..05_*.csv   - the five numbered on-request
      tables, in first-citation order (author decision 2026-09-23)
  results/tables/supplementary/strobe_mr_checklist.csv    - numberless, submitted separately
  results/tables/supplementary/_retired/*.csv             - four former tables whose content a
      figure already displays; kept because build_supplementary_figures.R reads one of them
All numbers are carried over from the primary analysis tables without re-computation.

2026-09-23 numbering (was -> now):
  S2 -> Table 1 (full 5,040-gene screen)      S1 -> Table 2 (22 WIDE hits)
  S5 -> Table 3 (SuSiE credible sets)         S6 -> Table 4 (PheWAS genome-wide hits)
  S3 -> Table 5 (cross-ancestry pQTL)         S10 -> strobe_mr_checklist.csv (no number)
  S4, S7, S8, S9 -> _retired/
"""
import pandas as pd
import numpy as np
import os

OUT = "results/tables/supplementary"
RETIRED = os.path.join(OUT, "_retired")
os.makedirs(RETIRED, exist_ok=True)
T = "results/tables"

def save(df, name):
    path = os.path.join(OUT, name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)
    print(name, df.shape)

# ----------------------------------------------------------------------
# Table 2 — 22 WIDE Bonferroni hits: MR, instrument strength, replication, coloc
# ----------------------------------------------------------------------
# Guard (2026-09-23): results/tables/susie_coloc_results.csv was overwritten by a failed
# SuSiE run and no longer carries pp_h4_max / n_pairs / conv_eqtl / conv_fg. Those values
# survive in results/tables/supplementary/_retired/coloc_abf_susie_44runs.csv (pp_h4_max and
# n_pairs are also recomputable from results/tables/susie_coloc_pairs.csv). Fail loudly here
# instead of raising a bare KeyError, and never silently drop a column from the table.
_sus_cols = pd.read_csv(f"{T}/susie_coloc_results.csv", nrows=0).columns
_missing = [c for c in ("pp_h4_max", "n_pairs", "conv_eqtl", "conv_fg") if c not in _sus_cols]
if _missing:
    raise SystemExit(
        f"susie_coloc_results.csv no longer carries {_missing}; it was overwritten by a "
        f"failed SuSiE run. The values survive in "
        f"results/tables/supplementary/_retired/coloc_abf_susie_44runs.csv - recalibrate from "
        f"results/tables/susie_coloc_pairs.csv (pp_h4_max, n_pairs) and re-run "
        f"run_susie_coloc.R for the convergence flags before regenerating this table. "
        f"See results/tables/supplementary/_retired/README.md.")
summ = pd.read_csv(f"{T}/mr_coloc_summary.csv")            # 22 rows
iv   = pd.read_csv(f"{T}/iv_strength_F.csv")               # 22 rows (WIDE hits)
xco  = pd.read_csv(f"{T}/mr_olp_hits_crosscohort.csv")
sus  = pd.read_csv(f"{T}/susie_coloc_results.csv")[["gene", "cohort", "pp_h4_max"]]
sus_wide = sus[sus.cohort == "WIDE"].rename(columns={"pp_h4_max": "susie_PP.H4_WIDE"})[["gene", "susie_PP.H4_WIDE"]]
sus_incl = sus[sus.cohort == "INCLAVO"].rename(columns={"pp_h4_max": "susie_PP.H4_INCLAVO"})[["gene", "susie_PP.H4_INCLAVO"]]

st1 = summ.merge(iv[["gene", "iv", "effect_allele", "other_allele", "allele_state",
                     "effect_allele_freq", "maf", "z_exp", "n_samples_eqtlgen", "r2_exp", "F"]],
                 on="gene", how="left")
st1 = st1.merge(xco[["gene", "dir_consistent", "repl_nominal", "repl_bonf"]], on="gene", how="left")
st1 = st1.merge(sus_wide, on="gene", how="left").merge(sus_incl, on="gene", how="left")

# CI columns are not in mr_coloc_summary -> bring from the raw MR tables
mrw = pd.read_csv(f"{T}/mr_olp_wide_results.csv")[["gene", "OR", "OR_lo", "OR_hi", "p_MR"]].rename(
    columns={"OR": "MR_OR_WIDE", "OR_lo": "ORlo_WIDE", "OR_hi": "ORhi_WIDE", "p_MR": "MR_P_WIDE"})
mri = pd.read_csv(f"{T}/mr_olp_inclavo_results.csv")[["gene", "OR", "OR_lo", "OR_hi", "p_MR"]].rename(
    columns={"OR": "MR_OR_INCLAVO", "OR_lo": "ORlo_INCLAVO", "OR_hi": "ORhi_INCLAVO", "p_MR": "MR_P_INCLAVO"})
summ = summ.drop(columns=["OR_WIDE", "P_WIDE", "OR_INCLAVO", "P_INCLAVO"]).rename(
    columns={"PPH4_WIDE": "abf_PP.H4_WIDE", "PPH4_INCLAVO": "abf_PP.H4_INCLAVO"})
st1 = summ.merge(iv[["gene", "iv", "effect_allele", "other_allele", "allele_state",
                     "effect_allele_freq", "maf", "z_exp", "n_samples_eqtlgen", "r2_exp", "F"]],
                 on="gene", how="left")
st1 = st1.merge(xco[["gene", "dir_consistent", "repl_nominal", "repl_bonf"]], on="gene", how="left")
st1 = st1.merge(sus_wide, on="gene", how="left").merge(sus_incl, on="gene", how="left")
st1 = st1.merge(mrw, on="gene", how="left").merge(mri, on="gene", how="left")
# variant_grch38 exists only in final_target_evidence (5 genes) -> fill from IV table if absent
save(st1, "suppl_table_02_wide_hits_mr_replication_coloc.csv")

# ----------------------------------------------------------------------
# Table 1 — full 5,040-gene MR results, both definitions
# ----------------------------------------------------------------------
F_N = 5040.0
st2 = mrw.merge(mri, on="gene", how="outer").rename(columns={"OR_WIDE": "OR_WIDE", "OR_INCLAVO": "OR_INCLAVO"})
mrw_full = pd.read_csv(f"{T}/mr_olp_wide_results.csv")
mri_full = pd.read_csv(f"{T}/mr_olp_inclavo_results.csv")
keep_w = ["gene", "snp", "chr", "pos", "effect_allele", "other_allele", "zscore", "pvalue",
          "NrSamples", "MAF", "Z_MR", "p_MR", "OR", "OR_lo", "OR_hi", "p_bonf", "p_fdr", "significant_bonf"]
keep_i = ["gene", "Z_MR", "p_MR", "OR", "OR_lo", "OR_hi", "p_bonf", "p_fdr", "significant_bonf"]
w = mrw_full[keep_w].rename(columns={"snp": "IV_rsID", "chr": "chr37", "pos": "pos37",
                                     "effect_allele": "EA_eQTL", "other_allele": "OA_eQTL",
                                     "zscore": "Z_eQTL", "pvalue": "P_eQTL", "NrSamples": "N_eQTL",
                                     "MAF": "MAF_EUR", "Z_MR": "Z_WIDE", "p_MR": "P_WIDE",
                                     "OR": "OR_WIDE", "OR_lo": "ORlo_WIDE", "OR_hi": "ORhi_WIDE",
                                     "p_bonf": "Pbonf_WIDE", "p_fdr": "PFDR_WIDE", "significant_bonf": "SigBonf_WIDE"})
i = mri_full[keep_i].rename(columns={"Z_MR": "Z_INCLAVO", "p_MR": "P_INCLAVO", "OR": "OR_INCLAVO",
                                     "OR_lo": "ORlo_INCLAVO", "OR_hi": "ORhi_INCLAVO",
                                     "p_bonf": "Pbonf_INCLAVO", "p_fdr": "PFDR_INCLAVO", "significant_bonf": "SigBonf_INCLAVO"})
st2 = w.merge(i, on="gene", how="outer")
# F-statistic approximation, same formula as iv_strength_F: R2 = Z^2/N ; F = (N-2)R2/(1-R2)
st2["F_approx"] = (st2["N_eQTL"] - 2) * (st2["Z_eQTL"] ** 2 / st2["N_eQTL"]) / (1 - st2["Z_eQTL"] ** 2 / st2["N_eQTL"])
st2 = st2[["gene", "IV_rsID", "chr37", "pos37", "EA_eQTL", "OA_eQTL", "MAF_EUR", "N_eQTL", "Z_eQTL", "P_eQTL",
           "F_approx", "OR_WIDE", "ORlo_WIDE", "ORhi_WIDE", "Z_WIDE", "P_WIDE", "Pbonf_WIDE", "PFDR_WIDE", "SigBonf_WIDE",
           "OR_INCLAVO", "ORlo_INCLAVO", "ORhi_INCLAVO", "Z_INCLAVO", "P_INCLAVO", "Pbonf_INCLAVO", "PFDR_INCLAVO", "SigBonf_INCLAVO"]]
assert len(st2) == int(F_N), f"expected 5040 rows, got {len(st2)}"
save(st2, "suppl_table_01_full_mr_screen_5040.csv")

# ----------------------------------------------------------------------
# Table 5 — pQTL assessment (as-is; 4 rows)
# ----------------------------------------------------------------------
st3 = pd.read_csv(f"{T}/pqtl_assessment.csv")
save(st3, "suppl_table_05_pqtl_assessment.csv")

# ----------------------------------------------------------------------
# RETIRED Table 4 — coloc.abf full (44 runs) + SuSiE summary (source data for Supplementary Figure 3)
# ----------------------------------------------------------------------
abf = pd.read_csv(f"{T}/coloc_olp_results.csv")
sus_cols = ["gene", "cohort", "n_snps", "n_cs_eqtl", "n_cs_fg", "n_cs_eqtl_full",
            "n_pairs", "pp_h4_max", "conv_eqtl", "conv_fg"]
sus44 = pd.read_csv(f"{T}/susie_coloc_results.csv")[sus_cols]

# Carry forward what only the retired copy still holds (2026-09-23).
# `susie_coloc_results.csv` was overwritten by a failed SuSiE run; `rebuild_susie_coloc_results.py`
# restored the pair-derived columns, but the convergence diagnostics (and some credible-set
# counts) are not derivable from the pair table and survive only in this retired copy. Never
# blank a value that exists there: fill nulls from it and say so. New non-null values always win.
RETIRED_44 = os.path.join(OUT, "_retired", "coloc_abf_susie_44runs.csv")
if os.path.exists(RETIRED_44):
    prev = pd.read_csv(RETIRED_44)
    src = {"pp_h4_max": "susie_PP.H4_max"}
    carry_cols = [c for c in sus44.columns if c not in ("gene", "cohort")]
    prev = prev.rename(columns=src)
    carried = {}
    for c in carry_cols:
        if c not in prev.columns:
            continue
        filled = sus44[["gene", "cohort", c]].merge(
            prev[["gene", "cohort", c]], on=["gene", "cohort"], how="left",
            suffixes=("", "_prev"))[f"{c}_prev"]
        need = sus44[c].isna() & filled.notna()
        if need.any():
            # cast through object so filling a boolean/NaN column does not warn about
            # assigning into an incompatible dtype
            sus44[c] = sus44[c].astype(object).where(sus44[c].notna(), filled.astype(object))
            carried[c] = int(need.sum())
    if carried:
        print("WARNING: carried forward from the retired copy of this table "
              "(their source was destroyed by the 2026-09-23 overwrite):",
              ", ".join(f"{k} x{v}" for k, v in sorted(carried.items())))
        print("         see results/tables/_susie_coloc_results_recovery_note.md")
    # keep integer counts integer-valued in the CSV (empty where unknown) so regenerating this
    # table reproduces the pre-overwrite file instead of printing 12.0. The fill above turns
    # these columns into `object`, so test for anything that is not already integer.
    # `n_pairs` is deliberately NOT cast: it has no value for the six non-estimable runs, so the
    # pre-overwrite file wrote it as a float (2.0) and a cast would change every row.
    for c in ("n_snps", "n_cs_eqtl", "n_cs_fg", "n_cs_eqtl_full"):
        if c in sus44.columns and sus44[c].dtype.kind not in ("i", "u"):
            sus44[c] = pd.to_numeric(sus44[c], errors="coerce").astype("Int64")
st4 = abf.merge(sus44, on=["gene", "cohort"], how="left", suffixes=("_abf", "_sus"))
st4["coloc_call"] = np.select([st4["PP.H4"] > 0.8, st4["PP.H4"] > 0.5],
                              ["strong", "moderate"], default="weak")
st4 = st4.rename(columns={"n_snps_abf": "n_snps", "pp_h4_max": "susie_PP.H4_max"})
st4 = st4.sort_values(["cohort", "PP.H4"], ascending=[True, False])
save(st4, "_retired/coloc_abf_susie_44runs.csv")

# ----------------------------------------------------------------------
# Table 3 — SuSiE credible sets (756 rows, as-is)
# ----------------------------------------------------------------------
st5 = pd.read_csv(f"{T}/susie_coloc_cs.csv")
save(st5, "suppl_table_03_susie_credible_sets.csv")

# ----------------------------------------------------------------------
# Table 4 — PheWAS genome-wide hits (121 rows, as-is)
# ----------------------------------------------------------------------
st6 = pd.read_csv(f"{T}/phewas_olp_gw_hits.csv")
save(st6, "suppl_table_04_phewas_gw_hits.csv")

# ----------------------------------------------------------------------
# RETIRED Table 7 — pleiotropy verdict: 5 candidates x 9 phenotypes (Figure 5C/D display this)
# ----------------------------------------------------------------------
pl = pd.read_csv(f"{T}/pleiotropy_coloc_results.csv")
sp = pd.read_csv(f"{T}/susie_pleiotropy_results.csv")
cand5 = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
pl5 = pl[pl.gene.isin(cand5)].copy()
# phenotypes used in Fig7B (OLP included as internal control)
phenos = sorted(pl5["pheno"].unique())
print("Fig7B phenotypes (n=%d):" % len(phenos), phenos)
pl5 = pl5[pl5.pheno.isin(phenos)].merge(
    sp[["gene", "pheno", "pp_h4_susie", "n_cs_fg", "hit1", "hit2"]],
    on=["gene", "pheno"], how="left")
pl5["susie_verdict"] = np.where(
    pl5["pp_h4_susie"].isna(), "no-susie-run",
    np.where(pl5["pp_h4_susie"] >= 0.8, "susie>=0.8", "susie<0.8"))
pl5["abf_verdict"] = np.where(pl5["PP.H4"] >= 0.8, "abf>=0.8", "abf<0.8")
pl5["final_verdict"] = np.select(
    [(pl5["PP.H4"] >= 0.8) & (pl5["pp_h4_susie"] >= 0.8), (pl5["PP.H4"] < 0.8) & (pl5["pp_h4_susie"] >= 0.8),
     (pl5["PP.H4"] >= 0.8) & (pl5["pp_h4_susie"] < 0.8) & pl5["pp_h4_susie"].notna()],
    ["confirmed", "susie-only", "abf-only-overturned"], default="neither")
pl5 = pl5[["gene", "pheno", "nsnp", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "abf_verdict",
           "pp_h4_susie", "susie_verdict", "final_verdict", "n_cs_fg", "hit1", "hit2"]]
save(pl5, "_retired/pleiotropy_verdict_5x9.csv")

# ----------------------------------------------------------------------
# RETIRED Table 8 — molecular docking full results (Figure 3B labels every passing ligand)
# ----------------------------------------------------------------------
st8 = pd.read_csv(f"{T}/docking_results.csv")
rec = {"GDF11": "5E4G (apo, chain A)", "SUOX": "2A9A (apo, chain A)"}
st8["receptor_PDB"] = st8["target"].map(rec)
st8["good_docking_lt_-7"] = st8["best_affinity_kcal_mol"] < -7
st8["affinity_rank"] = st8.groupby("target")["best_affinity_kcal_mol"].rank(ascending=False).astype(int)
st8 = st8.sort_values(["target", "best_affinity_kcal_mol"]).reset_index(drop=True)
save(st8, "_retired/docking_results.csv")

# ----------------------------------------------------------------------
# RETIRED Table 9 — Steiger direction test, all (50 rows; Figure 5B covers OLP WIDE)
# ----------------------------------------------------------------------
st9 = pd.read_csv(f"{T}/steiger_direction_test.csv")
save(st9, "_retired/steiger_direction.csv")

# ----------------------------------------------------------------------
# STROBE-MR checklist — no table number; submitted as a separate file
# ----------------------------------------------------------------------
strobe = [
    ["Title/abstract", 1, "Title and abstract indicate MR design", "Partial", "Manuscript stage: structured abstract planned (eQTLGen/FinnGen R13 dual-cohort, main finding GDF11)"],
    ["Introduction", 2, "Background and scientific rationale", "Yes", "literature_analysis_1_2.md; research_plan_oral_disease_v2.md"],
    ["Introduction", 3, "Objectives / pre-specified hypotheses", "Yes", "Main report §1 (dual-cohort cis-eQTL -> OLP druggable-genome screen); no prospective registry (author decision)"],
    ["Methods", "4a", "Data sources: population/time/design", "Yes", "eQTLGen 31,684 European whole blood; FinnGen R13 500,186 Finnish; OLP WIDE 9,526 / INCLAVO 2,491 cases"],
    ["Methods", "4b", "Data availability statement", "Partial", "URLs available; deCODE controlled-access + Eldjarn 2023 citation to be consolidated into manuscript Data Availability section"],
    ["Methods", "5a", "Instrument (IV) selection", "Yes", "Single top cis-eQTL per gene (eQTLGen, P<5e-8), cis window +/-1 Mb, min-P SNP; single-IV Wald ratio pre-specified"],
    ["Methods", "5b", "Instrument strength (F statistic)", "Yes", "iv_strength_F.csv; all F>10 (5 candidates F=61.4-2,928.2)"],
    ["Methods", 6, "Outcome GWAS source and definition", "Yes", "FinnGen R13 two OLP definitions (GRCh38); 3,724 allele flips, 26 mismatches excluded"],
    ["Methods", 7, "Sample overlap / two-sample assumption", "Partial", "eQTLGen and FinnGen partially overlapping European ancestry not quantified; disclaimer to be added in Methods"],
    ["Methods", "8a", "Main statistical method", "Yes", "Wald ratio identity beta_MR=beta_out/beta_exp, Z_MR=Z_out*sign(Z_exp); Bonferroni alpha=0.05/5,040=9.92e-6; dual-cohort replication"],
    ["Methods", "8b", "Assumption 1 (IV-exposure) checked", "Yes", "eQTLGen P<5e-8 (Z 4.5-47); F>10 for all hits"],
    ["Methods", "8c", "Assumption 2 (no confounding) checked", "Yes", "Steiger direction test (5 candidates expression->outcome; GDF11 P=6.7e-3, SUOX/TAPBP/HLA-DPA1 P~0, EGFL8 P=0.082) + PheWAS reverse screen"],
    ["Methods", "8d", "Assumption 3 (no horizontal pleiotropy) checked", "Yes", "coloc.abf (44 runs) + coloc.susie multi-signal re-analysis (44 runs Step 4 + chr12 18 runs Step 5c)"],
    ["Methods", "8e", "Sensitivity / robustness", "Yes", "Dual-cohort replication (19/22 nominal, 4/22 Bonferroni); SuSiE verdict spectrum; internal positive control PP.H4 delta=0"],
    ["Methods", "8f", "Sample size / power", "Partial", "N, s, case counts reported; formal single-IV power calculation not performed (optional)"],
    ["Methods", 9, "Software and versions", "Yes", "R 4.4.2, coloc 5.2.3, susieR 0.12.42, ggplot2 3.5.1, Python 3.14.7, AutoDock Vina 1.2.7"],
    ["Methods", 10, "Missing data / variant exclusions", "Yes", "26 rsID mismatches removed; 12 genes without NrSamples removed; indels/multiallelic and MAF<1% excluded (pQTL)"],
    ["Results", 11, "Descriptive data (participants, IVs)", "Yes", "5,459 genes with cis-eQTL / 5,040 tested; 1 IV per gene; 9,526 / 2,491 cases"],
    ["Results", 12, "IV-assumption assessment", "Partial", "Steiger + coloc reported; F now provided (see 5b)"],
    ["Results", 13, "Main results (OR/CI/P/sample size)", "Yes", "mr_olp_*_results.csv (5,040 rows); 22 WIDE hits; final_target_evidence.csv (5 candidates)"],
    ["Results", 14, "Sensitivity/direction/pleiotropy results", "Yes", "coloc PP.H4, SuSiE verdict, Steiger, PheWAS spectrum"],
    ["Results", 15, "Other analyses (subgroups etc.)", "N/A", "No subgroup/mediation analyses (not pre-specified); chr12 GDF11-SUOX LD non-independence documented via SuSiE/side-effect narrowing"],
    ["Discussion", 16, "Key results", "Yes", "Main report §5/§8; target_cards.md (GDF11 TIER1; SUOX/EGFL8 TIER2; TAPBP/HLA-DPA1 TIER3)"],
    ["Discussion", 17, "Limitations", "Yes", "Main report §7 (8 items: single IV, beta approximation, LD ancestry mismatch, no usable pQTL, apo docking, MHC pleiotropy)"],
    ["Discussion", 18, "Interpretation / mechanisms", "Yes", "MHC enrichment 41%; TAPBP antigen-processing mechanism; GDF11 side-effect profile (T1D/hypothyroid/nasal polyp)"],
    ["Discussion", 19, "Generalisability", "Partial", "To state: European-ancestry derivation; cross-ancestry replication with pQTL warranted"],
    ["Other", "20a", "Funding / conflicts of interest", "No", "Author input required (None if unfunded)"],
    ["Other", "20b", "Data availability + access", "Partial", "URLs in audit §1; deCODE controlled-access statement needed"],
    ["Other", "20c", "Code / reproducibility scripts", "Partial", "analysis/scripts (60+ scripts); environment lock + README input-output mapping recommended (Zenodo/scripts appendix)"],
    ["Other", "20d", "Ethics / informed consent", "No", "Author-confirmed wording required for FinnGen/eQTLGen/INTERVAL/Knight_ADRC/1000G; public/controlled summary statistics typically exempt"],
]
st10 = pd.DataFrame(strobe, columns=["domain", "item", "checklist_item", "status", "evidence_or_location"])
save(st10, "strobe_mr_checklist.csv")

print("Done. Numbered tables + checklist ->", OUT, "; retired ->", RETIRED)
