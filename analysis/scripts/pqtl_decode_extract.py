# Extract ALL cis pQTL SNPs with P < 1e-5 in GDF11/SUOX merged window (chr12 GRCh38).
# SUOX top p=1.3e-7 (rs747931940) is cis-Bonferroni significant; SUOX pQTL-MR is feasible.
import gzip, os, csv

BASE = r"remote/deCODE-largescaleplasma-2023"
OUT_SUOX = r"results/tables/pqtl_decode_suox_sig.csv"
OUT_GDF11 = r"results/tables/pqtl_decode_gdf11_sig.csv"

FILE_SUOX = os.path.join(BASE, "SomaScan_SUOX_SeqId18324.txt.gz")
FILE_GDF11 = os.path.join(BASE, "SomaScan_GDF11_SeqId14587.txt.gz")
FILE_GDF11_MSTN = os.path.join(BASE, "SomaScan_GDF11_SeqId2765_MSTN.txt.gz")

LO, HI = 54_700_000, 57_100_000
THRESH = 1e-5


def scan(path, thresh):
    hits = []
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] != "chr12":
                continue
            pos = int(p[1])
            if not (LO <= pos <= HI):
                continue
            pval = float(p[7])
            if pval <= thresh:
                rec = dict(zip(hdr, p))
                rec["pos_int"] = pos
                hits.append(rec)
    return hits


def write(hits, out, label, thresh):
    hits.sort(key=lambda r: float(r["Pval"]))
    with open(out, "w", newline="") as f:
        cols = ["pqtl", "Chrom", "pos_int", "rsids", "effectAllele", "otherAllele",
                "Beta", "SE", "Pval", "N", "ImpMAF"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for h in hits:
            w.writerow({c: h.get(c, "") for c in cols} | {"pqtl": label})
    print(f"{label}: {len(hits)} SNPs with P<={thresh:g} -> {out}")


def main():
    suox = scan(FILE_SUOX, THRESH)
    write(suox, OUT_SUOX, "SUOX_SeqId18324", THRESH)
    g1 = scan(FILE_GDF11, THRESH)
    write(g1, OUT_GDF11, "GDF11_SeqId14587", THRESH)
    g2 = scan(FILE_GDF11_MSTN, THRESH)
    write(g2, OUT_GDF11.replace("14587", "2765"), "GDF11_SeqId2765_MSTN", THRESH)
    # overlap with FinnGen coloc region? report matched rsids
    import fnmatch
    for fn in ["data/processed/coloc/coloc_fg_wide.csv", "data/processed/coloc/coloc_fg_inclavo.csv"]:
        with open(fn) as f:
            fg_rs = {r["snp"] for r in csv.DictReader(f)}
        for out, label in [(OUT_SUOX, "SUOX"), (OUT_GDF11, "GDF11"), (OUT_GDF11.replace("14587", "2765"), "GDF11_MSTN")]:
            with open(out) as f:
                rs = {r["rsids"] for r in csv.DictReader(f) if r["rsids"].startswith("rs")}
            common = rs & fg_rs
            print(f"{label} sig rsids in {fn.split('/')[-1]}: {len(common)} -> {sorted(common)[:10]}")


if __name__ == "__main__":
    main()
