"""Download 1000G Phase 3 VCF headers (HTTP Range + multi-member gzip) and
extract per-chromosome sample lists (VCF column order).
Usage: python fetch_1000g_vcf_headers.py"""
import urllib.request, gzip, os, sys

outdir = r"D:\AI-Wisp-Science\Drug Target MR\data\raw\1000g"
os.makedirs(outdir, exist_ok=True)
base = ("https://1000genomes.s3.amazonaws.com/release/20130502/"
        "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

chrs = [2, 5, 6, 11, 12, 17]
for c in chrs:
    url = base.format(c)
    fn = os.path.join(outdir, f"chr{c}_samples.txt")
    if os.path.exists(fn):
        print(f"chr{c}: cached"); continue
    req = urllib.request.Request(url, headers={"Range": "bytes=0-1572863", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        raw = r.read()
    txt = gzip.decompress(raw).decode("utf-8", errors="replace")
    chrom_lines = [l for l in txt.split("\n") if l.startswith("#CHROM")]
    if not chrom_lines:
        print(f"chr{c}: #CHROM not found in first 1.5MB; FAIL"); continue
    samples = chrom_lines[0].split("\t")[9:]
    with open(fn, "w") as f:
        f.write("\n".join(samples))
    print(f"chr{c}: {len(samples)} samples -> {fn}")
print("Done.")
