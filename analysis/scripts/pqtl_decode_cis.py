# Scan deCODE SomaScan pQTL files (GRCh38) for cis signals at GDF11 / SUOX.
# Genes (GRCh38): GDF11 chr12:55,743,122-55,757,264 ; SUOX chr12:55,996,855-56,006,641
# Merged scan window covers both genes +/- 1 Mb: chr12:54,700,000-57,100,000.
# Usage: python pqtl_decode_cis.py
import gzip, os, csv, sys

BASE = r"remote/deCODE-largescaleplasma-2023"
OUT = r"results/tables/pqtl_decode_cis_top.csv"

FILES = {
    "GDF11_SeqId14587": os.path.join(BASE, "SomaScan_GDF11_SeqId14587.txt.gz"),
    "GDF11_SeqId2765_MSTN": os.path.join(BASE, "SomaScan_GDF11_SeqId2765_MSTN.txt.gz"),
    "SUOX_SeqId18324": os.path.join(BASE, "SomaScan_SUOX_SeqId18324.txt.gz"),
}

GENES = {  # GRCh38
    "GDF11": (55_743_122, 55_757_264),
    "SUOX": (55_996_855, 56_006_641),
}
LO, HI = 54_700_000, 57_100_000
GW_SIG = 5e-8


def assign_gene(pos):
    hits = []
    for g, (s, e) in GENES.items():
        if s - 1_000_000 <= pos <= e + 1_000_000:
            hits.append(g)
    return "|".join(hits) if hits else "outside"


def scan(path):
    top = []  # (pval, dict)
    n_chr12 = 0
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] != "chr12":
                continue
            n_chr12 += 1
            pos = int(p[1])
            if not (LO <= pos <= HI):
                continue
            pval = float(p[7])
            if pval > 1e-4:  # keep only suggestive cis candidates
                continue
            rec = dict(zip(hdr, p))
            top.append((pval, pos, rec))
    top.sort(key=lambda x: x[0])
    return top, n_chr12


def main():
    os.makedirs("results/tables", exist_ok=True)
    rows = []
    summary = []
    for label, path in FILES.items():
        print(f"scanning {label} ...", flush=True)
        top, n_chr12 = scan(path)
        summary.append((label, n_chr12, len(top)))
        for pval, pos, rec in top[:15]:
            g = assign_gene(pos)
            rows.append({
                "pqtl": label, "file": os.path.basename(path),
                "rsid": rec.get("rsids"), "chrom": rec.get("Chrom"), "pos_grch38": pos,
                "gene_cis": g, "effectAllele": rec.get("effectAllele"),
                "otherAllele": rec.get("otherAllele"), "Beta": rec.get("Beta"),
                "SE": rec.get("SE"), "Pval": rec.get("Pval"), "N": rec.get("N"),
                "ImpMAF": rec.get("ImpMAF"),
            })
    for label, n, k in summary:
        print(f"{label}: chr12 lines={n}, suggestive(p<1e-4) in window={k}")
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT, f"({len(rows)} rows)")


if __name__ == "__main__":
    main()
