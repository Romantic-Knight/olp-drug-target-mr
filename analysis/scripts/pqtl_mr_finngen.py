# SUOX/GDF11 pQTL-MR: deCODE SomaScan pQTL (GRCh38, N=35,652) -> FinnGen R13 OLP (WIDE + INCLAVO).
# Strict MR identity as in the rest of the pipeline:
#   beta_MR = beta_out / beta_exp ; Z_MR = Z_out * sign(Z_exp) ; p_MR = 2*pnorm(-|Z_MR|)
# Effect allele alignment: deCODE effectAllele vs FinnGen ref/alt (beta is for alt in FinnGen).
import os, sys, math, csv
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from fetch_finngen_region import fetch_region

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
os.chdir(ROOT)

OUT = "results/tables/pqtl_mr_olp_results.csv"
PHENOS = {
    "K11_LICHEN_PLANUS_WIDE": "WIDE",
    "K11_ORAL_LICHEN_PLANUS_INCLAVO": "INCLAVO",
}
# deCODE pQTL sig SNPs (p<=1e-5) extracted by pqtl_decode_extract.py
PQTL_FILES = {
    "SUOX_SeqId18324": "results/tables/pqtl_decode_suox_sig.csv",
    "GDF11_SeqId14587": "results/tables/pqtl_decode_gdf11_sig.csv",
}
WINDOW = (55_500_000, 56_500_000)  # GRCh38 chr12 covering GDF11 3'/SUOX gene + pQTL cluster

def pnorm2(z):
    # two-sided normal tail p
    import math
    return math.erfc(abs(z) / math.sqrt(2))

def load_pqtl():
    rows = []
    for label, fn in PQTL_FILES.items():
        with open(fn) as f:
            for r in csv.DictReader(f):
                if r["rsids"].startswith("rs"):
                    rows.append({
                        "pqtl": label, "rsid": r["rsids"], "pos": int(r["pos_int"]),
                        "exp_effect": r["effectAllele"], "exp_other": r["otherAllele"],
                        "exp_beta": float(r["Beta"]), "exp_se": float(r["SE"]),
                        "exp_p": float(r["Pval"]), "exp_maf": float(r["ImpMAF"]),
                    })
    return rows

def main():
    pqtl = load_pqtl()
    print(f"{len(pqtl)} pQTL sig SNPs loaded")
    fg_cache = {}
    results = []
    for pheno, cohort in PHENOS.items():
        if pheno not in fg_cache:
            reg = fetch_region(pheno, "12", *WINDOW)
            fg_cache[pheno] = reg.set_index("snp") if not reg.empty else None
            print(f"[{pheno}] fetched {len(reg)} variants chr12:{WINDOW[0]:,}-{WINDOW[1]:,}")
        fg = fg_cache[pheno]
        if fg is None:
            print(f"[{pheno}] EMPTY"); continue
        for q in pqtl:
            rs = q["rsid"]
            if rs not in fg.index:
                # try comma-stripped first rsid
                print(f"  [{cohort}] {rs} ({q['pqtl']}) not in FinnGen region -> skip")
                continue
            f = fg.loc[rs]
            # allele alignment
            if f["alt"] == q["exp_effect"]:
                b_out, sign_align = float(f["beta"]), 1.0
            elif f["ref"] == q["exp_effect"]:
                b_out, sign_align = -float(f["beta"]), -1.0
            else:
                print(f"  [{cohort}] {rs}: allele mismatch fg {f['ref']}/{f['alt']} vs exp {q['exp_effect']}/{q['exp_other']} -> skip")
                continue
            se_out = float(f["sebeta"]); p_out = float(f["pval"])
            # strict MR identity; beta_out already aligned to exp effect allele
            b_exp, se_exp = q["exp_beta"], q["exp_se"]
            if b_exp == 0 or se_exp == 0:
                continue
            z_out = b_out / se_out if se_out > 0 else float("nan")
            b_mr = b_out / b_exp
            z_mr = z_out * (1 if b_exp > 0 else -1) * sign_align
            p_mr = pnorm2(z_mr)
            results.append({
                "cohort": cohort, "pheno": pheno, "gene": q["pqtl"].split("_")[0],
                "pqtl": q["pqtl"], "rsid": rs, "pos_grch38": q["pos"],
                "exp_effect": q["exp_effect"], "exp_other": q["exp_other"],
                "exp_beta": b_exp, "exp_se": se_exp, "exp_p": q["exp_p"], "exp_maf": q["exp_maf"],
                "fg_ref": f["ref"], "fg_alt": f["alt"], "fg_af_alt": f["af_alt"],
                "out_beta": b_out, "out_se": se_out, "out_p": p_out,
                "beta_mr": b_mr, "z_mr": z_mr, "p_mr": p_mr,
                "or_olp_per_sd": math.exp(b_mr),
            })
            print(f"  [{cohort}] {q['pqtl']} {rs}: exp_beta={b_exp:+.4f} out_beta={b_out:+.4f} "
                  f"-> beta_MR={b_mr:+.4f} p_MR={p_mr:.3e} OR={math.exp(b_mr):.3f}")
    with open(OUT, "w", newline="") as f:
        cols = list(results[0].keys())
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(results)
    print("wrote", OUT, f"({len(results)} rows)")

if __name__ == "__main__":
    main()
