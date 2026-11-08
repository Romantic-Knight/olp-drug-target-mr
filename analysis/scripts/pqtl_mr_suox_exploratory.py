# Exploratory SUOX pQTL-MR (deCODE SomaScan -> FinnGen R13 OLP WIDE/INCLAVO).
# Context: no biallelic common cis-pQTL reaches strict cis significance (all P>3e-6);
# these single-SNP Wald-ratio results are directional/exploratory only, NOT tier-changing.
# Matching: GRCh38 pos + allele concordance (direct or complement flip). Strict MR identity used.
import os, sys, math, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch_finngen_region import parse_tbi, reg2bins, fetch_bytes, decompress_bgzf

FG_DIR = r"D:\AI-Wisp-Science\Drug Target MR\data\raw\finngen_r13"
ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
os.chdir(ROOT)
OUT = "results/tables/pqtl_mr_suox_exploratory.csv"
WINDOW = (55_500_000, 56_800_000)

PHENOS = {
    "K11_LICHEN_PLANUS_WIDE": "WIDE",
    "K11_ORAL_LICHEN_PLANUS_INCLAVO": "INCLAVO",
}

def fetch_raw_variants(pheno, beg38, end38):
    gz_url = f"https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_{pheno}.gz"
    tbi_path = os.path.join(FG_DIR, f"finngen_R13_{pheno}.gz.tbi")
    refs = parse_tbi(tbi_path)
    bins_idx = refs["12"]
    chunks = []
    for b in reg2bins(beg38, end38):
        chunks.extend(bins_idx.get(b, []))
    chunks = sorted(set(chunks), key=lambda c: (c[0] >> 16, c[1] >> 16))
    merged = []
    for c in chunks:
        rs0, rs1 = c[0] >> 16, c[1] >> 16
        if merged and rs0 <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], rs1))
        else:
            merged.append([rs0, rs1])
    rows = {}
    for (s0, s1) in merged:
        raw = fetch_bytes(gz_url, s0, s1)
        txt = decompress_bgzf(raw)
        for line in txt.decode("utf-8", "replace").split("\n"):
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 13:
                continue
            pos = int(f[1])
            if beg38 <= pos <= end38:
                rows[pos] = {"pos": pos, "ref": f[2], "alt": f[3], "rsid": f[4],
                             "pval": float(f[6]), "beta": float(f[8]), "se": float(f[9]),
                             "af_alt": float(f[10])}
    return rows


def load_pqtl():
    rows = []
    with open("results/tables/pqtl_decode_SUOX_SeqId18324_biallelic.csv") as f:
        for r in csv.DictReader(f):
            rows.append({"gene": "SUOX", "rsid": r["rsids"], "pos": int(r["pos_int"]),
                         "eff": r["effectAllele"], "oth": r["otherAllele"],
                         "beta": float(r["Beta"]), "se": float(r["SE"]),
                         "p": float(r["Pval"]), "maf": float(r["ImpMAF"])})
    return rows


def align(fg, q):
    """Return (out_beta_aligned, sign_align) or None if no allele concordance."""
    eff, oth = q["eff"], q["oth"]
    for candidate_ref, candidate_alt in [(fg["ref"], fg["alt"]), (fg["alt"], fg["ref"])]:
        if candidate_alt == eff and candidate_ref == oth:
            return fg["beta"], 1.0
        if candidate_alt == oth and candidate_ref == eff:
            return -fg["beta"], -1.0
    return None


def main():
    pqtl = load_pqtl()
    print(f"{len(pqtl)} SUOX biallelic pQTL candidates")
    results = []
    for pheno, cohort in PHENOS.items():
        fg = fetch_raw_variants(pheno, *WINDOW)
        print(f"[{cohort}] FinnGen variants {len(fg)} in window")
        for q in pqtl:
            if q["pos"] not in fg:
                print(f"  {q['rsid']} pos {q['pos']:,}: NOT in FinnGen -> skip")
                continue
            f = fg[q["pos"]]
            aligned = align(f, q)
            if aligned is None:
                print(f"  {q['rsid']} pos {q['pos']:,}: allele mismatch fg {f['ref']}/{f['alt']} "
                      f"vs {q['eff']}/{q['oth']} -> skip")
                continue
            b_out, sa = aligned
            b_exp, se_exp = q["beta"], q["se"]
            z_out = b_out / f["se"]
            b_mr = b_out / b_exp
            z_mr = z_out * (1 if b_exp > 0 else -1) * sa
            p_mr = math.erfc(abs(z_mr) / math.sqrt(2))
            results.append({
                "cohort": cohort, "gene": "SUOX", "rsid": q["rsid"], "pos_grch38": q["pos"],
                "exp_eff": q["eff"], "exp_oth": q["oth"], "exp_beta": b_exp,
                "exp_se": se_exp, "exp_p": q["p"], "exp_maf": q["maf"],
                "fg_ref": f["ref"], "fg_alt": f["alt"], "fg_af_alt": f["af_alt"],
                "fg_rsid": f["rsid"], "out_beta": b_out, "out_se": f["se"], "out_p": f["pval"],
                "beta_mr": b_mr, "z_mr": z_mr, "p_mr": p_mr, "or_olp": math.exp(b_mr),
            })
            print(f"  [{cohort}] {q['rsid']} (fg_rs={f['rsid']}): exp_b={b_exp:+.4f} out_b={b_out:+.4f} "
                  f"-> beta_MR={b_mr:+.4f} p_MR={p_mr:.3e}")
    if results:
        with open(OUT, "w", newline="") as f:
            cols = list(results[0].keys())
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(results)
        print("wrote", OUT, f"({len(results)} rows)")
    else:
        print("no analyzable pairs")

if __name__ == "__main__":
    main()
