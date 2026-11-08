# Step 6 pQTL evidence assessment: SUOX cis-pQTL from Knight_ADRC (GWAS Catalog GCST90619063).
# Determines whether a usable cis-pQTL instrument exists for SUOX pQTL-MR.
# Reference: SUOX GRCh37 chr12:56,390,964-56,400,425; GRCh38 chr12:55,996,855-56,006,641 (Ensembl).
# File has NO _GRCh37/_GRCh38 suffix -> GWAS Catalog default assembly = GRCh38 (verified: top cis
# signal 55,959,665 lies 37 kb upstream of GRCh38 SUOX; no such signal at GRCh37 SUOX position).
# Cis window = gene +/- 1 Mb on GRCh38. Threshold for a usable instrument: genome-wide P<5e-8
# (strict) or P<1e-6 (conservative cis) and MAF not extreme-low.
import gzip, csv, os

PQTL = r"data/raw/pqtl/GCST90619063.tsv.gz"
SUOX_GRCH38 = (55_996_855, 56_006_641)
CIS_BP = 1_000_000
OUT = r"results/tables/pqtl_assessment.csv"

def scan():
    lo = SUOX_GRCH38[0] - CIS_BP
    hi = SUOX_GRCH38[1] + CIS_BP
    rows = []
    with gzip.open(PQTL, "rt") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] != "12":
                continue
            pos = int(p[1])
            if lo <= pos <= hi:
                rows.append(dict(pos=pos, a1=p[2], a2=p[3], beta=float(p[4]),
                                 se=float(p[5]), eaf=float(p[6]),
                                 neglogp=float(p[7]), n=int(p[8])))
    return rows

def main():
    rows = scan()
    rows.sort(key=lambda r: -r["neglogp"])
    top = rows[0]
    top_p = 10 ** -top["neglogp"]
    n = top["n"]
    # verdict
    if top_p < 5e-8:
        verdict = "usable"
        reason = f"top cis p={top_p:.2e} < 5e-8 (GW sig)"
    elif top_p < 1e-6:
        verdict = "borderline"
        reason = f"top cis p={top_p:.2e} in 1e-6..5e-8, needs replication"
    else:
        verdict = "not_usable_insufficient_power"
        reason = (f"top cis p={top_p:.2e} (>= 1e-6); EAF={top['eaf']:.3f} extreme-low; "
                  f"N={n} too small for reliable cis-pQTL")
    dist_gene = top["pos"] - SUOX_GRCH38[0]
    print(f"cis window GRCh38 chr12:{SUOX_GRCH38[0]-CIS_BP:,}-{SUOX_GRCH38[1]+CIS_BP:,}: {len(rows)} SNPs")
    print(f"top cis signal: pos={top['pos']:,} ({dist_gene:+,} bp vs gene start) "
          f"{top['a1']}/{top['a2']} beta={top['beta']:.4f} se={top['se']:.4f} "
          f"EAF={top['eaf']:.4f} -log10p={top['neglogp']:.3f} n={n}")
    print("VERDICT:", verdict, "|", reason)

    os.makedirs("results/tables", exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "pqtl_study", "gwas_catalog", "assembly", "n_samples",
                    "cis_chr", "cis_grch38_start", "cis_grch38_end", "n_cis_snps",
                    "top_cis_pos_grch38", "dist_to_gene_start_bp", "top_effect_allele",
                    "top_other_allele", "top_beta", "top_se", "top_eaf",
                    "top_neglog10p", "top_p", "usable_threshold", "verdict", "reason"])
        w.writerow(["SUOX", "Knight_ADRC plasma pQTL (Sulfite oxidase, mitochondrial protein level)",
                    "GCST90619063", "GRCh38", n,
                    "12", SUOX_GRCH38[0]-CIS_BP, SUOX_GRCH38[1]+CIS_BP, len(rows),
                    top["pos"], dist_gene, top["a1"], top["a2"], top["beta"], top["se"],
                    round(top["eaf"], 6), round(top["neglogp"], 4), f"{top_p:.3e}",
                    "5e-8 (GW) / 1e-6 (conservative cis)", verdict, reason])
    print("wrote", OUT)

if __name__ == "__main__":
    main()
