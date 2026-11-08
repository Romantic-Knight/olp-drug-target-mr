"""Steiger directionality test for the 5 candidate-gene IVs.

For each candidate IV we compare its correlation with the eQTL exposure
(r_exp, from eQTLGen z-score & sample size) against its correlation with each
outcome (r_out, from FinnGen beta/se and effective sample size):
  r = sqrt(z^2 / (z^2 + N))
A significantly larger r_exp than r_out supports expression -> disease
direction (the IV primarily tags the exposure), the assumption underlying
cis-eQTL MR. We test the difference with the Fisher z transformation.

Output: results/tables/steiger_direction_test.csv
"""
import pandas as pd, numpy as np, os, math

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
OLP_CASES = {"K11_LICHEN_PLANUS_WIDE": (9526, 490660), "K11_ORAL_LICHEN_PLANUS_INCLAVO": (2491, 497695)}
PHENO_META = {
    "K11_COELIAC": (5308, 476554), "T1D_WIDE": (11197, 396409),
    "E4_HYTHY_AI_STRICT": (56098, 358456), "J10_NASALPOLYP": (8793, 369050),
    "L12_LICHENSCLERATROPH": (4181, 462830), "E4_GRAVES_STRICT": (4247, 495939),
    "RHEUMA_SEROPOS_STRICT": (6358, 488366), "AUTOIMMUNE": (131428, 368758),
}
IVS = {"GDF11": "rs7302200", "SUOX": "rs773107", "EGFL8": "rs41268896",
       "TAPBP": "rs3106190", "HLA-DPA1": "rs114390018"}

# eQTLGen side: z and N for the IV from the coloc windows file
eq = pd.read_csv(os.path.join(ROOT, "data", "processed", "coloc", "coloc_eqtl_windows.csv"))

rows = []
for g, rsid in IVS.items():
    e = eq[(eq.gene == g) & (eq.snp == rsid)]
    if e.empty:
        print(f"WARN: IV {rsid} not in eQTL windows for {g}")
        continue
    z_exp = e.iloc[0]["zscore"]; n_exp = e.iloc[0]["NrSamples"]
    r_exp = math.sqrt(z_exp**2 / (z_exp**2 + n_exp))
    # FinnGen side from pleiotropy region files (all phenotypes incl. OLP)
    for pheno, (cases, ctrls) in {**OLP_CASES, **PHENO_META}.items():
        # load IV beta/se for this pheno
        p = os.path.join(ROOT, "data", "processed", "phewas", "coloc", f"fg_{pheno}_{g}.csv")
        if not os.path.exists(p):
            if pheno.startswith("K11_LICHEN_PLANUS_WIDE"):
                p = os.path.join(ROOT, "data", "processed", "coloc", "coloc_fg_wide.csv")
            elif pheno.startswith("K11_ORAL"):
                p = os.path.join(ROOT, "data", "processed", "coloc", "coloc_fg_inclavo.csv")
            else:
                continue
        fg = pd.read_csv(p)
        hit = fg[fg.snp == rsid]
        if hit.empty:
            continue
        b, se = hit.iloc[0]["beta"], hit.iloc[0]["sebeta"]
        n_eff = 4 * cases * ctrls / (cases + ctrls)
        z_out = b / se
        r_out = math.sqrt(z_out**2 / (z_out**2 + n_eff))
        # Fisher z comparison (Steiger 1980), two-sided p from |z|
        zf1 = 0.5 * math.log((1 + r_exp) / (1 - r_exp))
        zf2 = 0.5 * math.log((1 + r_out) / (1 - r_out))
        se_diff = math.sqrt(1 / (n_exp - 3) + 1 / (n_eff - 3))
        z_diff = (zf1 - zf2) / se_diff
        p_steiger = 1 - math.erf(abs(z_diff) / math.sqrt(2))
        direction = "expression->outcome" if r_exp > r_out else "outcome->expression"
        rows.append({"gene": g, "iv": rsid, "pheno": pheno, "n_cases": cases,
                     "r_exp": r_exp, "r_out": r_out, "direction": direction,
                     "p_steiger": p_steiger,
                     "z_exp": z_exp, "z_out": z_out})
res = pd.DataFrame(rows)
res = res[["gene", "iv", "pheno", "n_cases", "z_exp", "z_out", "r_exp", "r_out", "direction", "p_steiger"]]
res.to_csv(os.path.join(ROOT, "results", "tables", "steiger_direction_test.csv"), index=False)
print(res.to_string(index=False))
print("\nrows:", len(res))
# summary: how many IV-outcome pairs pass Steiger (p<0.05, expression->outcome)
passed = (res.p_steiger < 0.05) & (res.direction == "expression->outcome")
print("\nSteiger-supported (p<0.05, expression->outcome):", passed.sum(), "/", len(res))
print(res[passed].gene.value_counts().to_string())
