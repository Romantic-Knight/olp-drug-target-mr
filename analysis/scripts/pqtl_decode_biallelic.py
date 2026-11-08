# Extract biallelic common (MAF>0.01) cis pQTL SNPs (P<1e-4) for SUOX / GDF11 from deCODE SomaScan.
# Excludes indels (allele length > 1) which are unreliable as instruments (poly-A/T repeats).
import gzip, os, csv

BASE = r"remote/deCODE-largescaleplasma-2023"
OUT_DIR = r"results/tables"
LO, HI = 54_700_000, 57_100_000
THRESH = 1e-4
MIN_MAF = 0.01
FILES = {
    "SUOX_SeqId18324": os.path.join(BASE, "SomaScan_SUOX_SeqId18324.txt.gz"),
    "GDF11_SeqId14587": os.path.join(BASE, "SomaScan_GDF11_SeqId14587.txt.gz"),
    "GDF11_SeqId2765_MSTN": os.path.join(BASE, "SomaScan_GDF11_SeqId2765_MSTN.txt.gz"),
}


def is_biallelic_snp(eff, oth):
    if not eff or not oth:
        return False
    if len(eff) != 1 or len(oth) != 1:
        return False
    return {eff, oth} <= {"A", "C", "G", "T"} and eff != oth


def scan(path):
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
            if pval > THRESH:
                continue
            eff, oth = p[4], p[5]
            if not is_biallelic_snp(eff, oth):
                continue
            if float(p[11]) < MIN_MAF:
                continue
            rec = dict(zip(hdr, p))
            rec["pos_int"] = pos
            hits.append(rec)
    return hits


def write(hits, out, label):
    hits.sort(key=lambda r: float(r["Pval"]))
    with open(out, "w", newline="") as f:
        cols = ["pqtl", "Chrom", "pos_int", "rsids", "effectAllele", "otherAllele",
                "Beta", "SE", "Pval", "N", "ImpMAF"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for h in hits:
            w.writerow({c: h.get(c, "") for c in cols} | {"pqtl": label})
    print(f"{label}: {len(hits)} biallelic common cis SNPs (P<={THRESH:g}, MAF>{MIN_MAF}) -> {out}")


def main():
    for label, path in FILES.items():
        hits = scan(path)
        tag = label.split("_")[0] + "_" + label.split("_")[1]
        write(hits, os.path.join(OUT_DIR, f"pqtl_decode_{tag}_biallelic.csv"), label)


if __name__ == "__main__":
    main()
