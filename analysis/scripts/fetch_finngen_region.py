"""Fetch FinnGen R13 summary-stats regions for pleiotropy colocalisation.

FinnGen R13 summary stats (Google Cloud Storage, ~800 MB per phenotype) are
bgzip-compressed with a tabix index (.tbi). This script reads only the byte
spans overlapping the requested region via HTTP Range, using the same
tabix/BGZF machinery as build_ld_1000g.py, so no full-file download is needed.

For each (phenotype x gene) pair it writes a FinnGen-side coloc input file
(with the same schema as coloc_fg_wide.csv) restricted to rsids present in the
gene's 1000G LD panel:

  data/processed/phewas/coloc/fg_<pheno>_<gene>.csv   (snp, ref, alt, pval, beta, sebeta, af_alt)

Windows are derived from the LD SNP positions (GRCh37) lifted to GRCh38 using
the empirically verified liftoff at each locus (chr12: -393,784 bp;
chr6 MHC: +32,223 bp) plus a 200 kb safety margin; rsid matching between the
GRCh37 eQTL side and GRCh38 FinnGen side absorbs any remaining offset.

Usage: python fetch_finngen_region.py [pheno ...]   (default: all in PHENOS)
"""
import struct, gzip, os, sys, io, time
import pandas as pd
import urllib.request

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
FG_DIR = os.path.join(ROOT, "data", "raw", "finngen_r13")
OUT_DIR = os.path.join(ROOT, "data", "processed", "phewas", "coloc")
os.makedirs(OUT_DIR, exist_ok=True)

BASE = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_{}.gz"

# phenotype -> (num_cases, num_controls) from finngen_R13_manifest.tsv
MANIFEST = os.path.join(FG_DIR, "finngen_R13_manifest.tsv")
PHENOS = {
    "K11_COELIAC": (5308, 476554),
    "T1D_WIDE": (11197, 396409),
    "E4_HYTHY_AI_STRICT": (56098, 358456),
    "J10_NASALPOLYP": (8793, 369050),
    "L12_LICHENSCLERATROPH": (4181, 462830),
    "E4_GRAVES_STRICT": (4247, 495939),
    "RHEUMA_SEROPOS_STRICT": (6358, 488366),
    "AUTOIMMUNE": (131428, 368758),
}

# gene -> (chrom, GRCh37 pos min/max from LD panel, GRCh38 shift)
GENES = {
    "GDF11":    ("12", 56131052, 56508409, -393784),
    "SUOX":     ("12", 55419667, 57394635, -393784),
    "EGFL8":    ("6", 31321211, 32666960, 32223),
    "TAPBP":    ("6", 32577973, 33833963, 32223),
    "HLA-DPA1": ("6", 32045226, 33933375, 32223),
}
MARGIN = 200_000


# ---------------- tabix index parsing (same as build_ld_1000g.py) ----------------
def parse_tbi(path):
    with open(path, "rb") as f:
        d = f.read()
    if d[:2] == b"\x1f\x8b":
        d = gzip.decompress(d)
    assert d[:4] == b"TBI\x01", "not a tabix index"
    off = 4
    n_ref = struct.unpack_from("<i", d, off)[0]; off += 4
    off += 24
    l_nm = struct.unpack_from("<i", d, off)[0]; off += 4
    names = d[off:off + l_nm].decode("ascii", "replace"); off += l_nm
    ref_names = [x for x in names.replace("\x00", "\t").split("\t") if x]
    refs = {}
    for rn in ref_names:
        n_bin = struct.unpack_from("<i", d, off)[0]; off += 4
        bins = {}
        for _ in range(n_bin):
            bin_, n_chunk = struct.unpack_from("<Ii", d, off); off += 8
            chunks = []
            for _ in range(n_chunk):
                cbeg, cend = struct.unpack_from("<QQ", d, off); off += 16
                chunks.append((cbeg, cend))
            bins[bin_] = chunks
        n_intv = struct.unpack_from("<i", d, off)[0]; off += 4
        off += 8 * n_intv
        refs[rn] = bins
    return refs


def reg2bins(beg, end):
    bins = [0]
    e = end - 1
    bins += list(range(1 + (beg >> 26), 1 + (e >> 26) + 1))
    bins += list(range(9 + (beg >> 23), 9 + (e >> 23) + 1))
    bins += list(range(73 + (beg >> 20), 73 + (e >> 20) + 1))
    bins += list(range(585 + (beg >> 17), 585 + (e >> 17) + 1))
    bins += list(range(4681 + (beg >> 14), 4681 + (e >> 14) + 1))
    bins += list(range(37449 + (beg >> 11), 37449 + (e >> 11) + 1))
    return bins


def decompress_bgzf(data, maxout=1_000_000_000):
    out = bytearray()
    off = 0
    while off + 18 <= len(data) and len(out) < maxout:
        if data[off:off + 2] != b"\x1f\x8b":
            break
        bsize = int.from_bytes(data[off + 16:off + 18], "little")
        block = data[off:off + bsize + 1]
        if len(block) < bsize + 1:
            break
        try:
            out += gzip.decompress(block)
        except Exception:
            break
        off += bsize + 1
    return bytes(out)


def fetch_bytes(url, start, end):
    req = urllib.request.Request(url, headers={
        "Range": f"bytes={start}-{end}", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return r.read()


# ---------------- main ----------------
def fetch_region(pheno, chrom, beg38, end38):
    """Return DataFrame of variants in [beg38, end38] on chrom (GRCh38)."""
    gz_url = BASE.format(pheno)
    tbi_path = os.path.join(FG_DIR, f"finngen_R13_{pheno}.gz.tbi")
    if not os.path.exists(tbi_path):
        req = urllib.request.Request(gz_url + ".tbi", headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=300) as r:
            tbi = r.read()
        with open(tbi_path, "wb") as f:
            f.write(tbi)
        print(f"  [{pheno}] downloaded tbi ({len(tbi)} bytes)")
    refs = parse_tbi(tbi_path)
    ref_name = chrom if chrom in refs else None
    if ref_name is None:
        print(f"  [{pheno}] chrom {chrom} not in index"); return pd.DataFrame()
    bins_idx = refs[ref_name]
    chunks = []
    for b in reg2bins(beg38, end38):
        chunks.extend(bins_idx.get(b, []))
    chunks = sorted(chunks, key=lambda c: (c[0] >> 16, c[1] >> 16))
    merged = []
    for c in chunks:
        rs0, rs1 = c[0] >> 16, c[1] >> 16
        if merged and rs0 <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], rs1))
        else:
            merged.append([rs0, rs1])
    rows = []
    for (s0, s1) in merged:
        raw = fetch_bytes(gz_url, s0, s1)
        txt = decompress_bgzf(raw)
        if not txt:
            continue
        lines = txt.decode("utf-8", "replace").split("\n")
        for line in lines:
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 13:
                continue
            pos = int(f[1])
            if not (beg38 <= pos <= end38):
                continue
            rsids = f[4]
            if rsids == "." or "," in rsids or "," in f[3]:
                continue  # require a single rsid and biallelic
            try:
                pval = float(f[6]); beta = float(f[8]); se = float(f[9]); maf = float(f[10])
            except ValueError:
                continue
            rows.append((rsids, f[2], f[3], pval, beta, se, maf))
    return pd.DataFrame(rows, columns=["snp", "ref", "alt", "pval", "beta", "sebeta", "af_alt"])


def main():
    args = sys.argv[1:]
    phenos = args if args else list(PHENOS.keys())
    # pre-load LD panel rsid sets
    ld_snps = {}
    for g, (chrom, p0, p1, shift) in GENES.items():
        p = os.path.join(ROOT, "data", "processed", "coloc", "ld_1000g", f"{g}_snps.csv")
        ld_snps[g] = set(pd.read_csv(p)["rsid"])

    for pheno in phenos:
        if pheno not in PHENOS:
            print("unknown pheno:", pheno); continue
        for g, (chrom, p0, p1, shift) in GENES.items():
            beg38 = (p0 + shift) - MARGIN
            end38 = (p1 + shift) + MARGIN
            out = os.path.join(OUT_DIR, f"fg_{pheno}_{g}.csv")
            if os.path.exists(out):
                print(f"[{pheno}][{g}] cached"); continue
            t0 = time.time()
            reg = fetch_region(pheno, chrom, beg38, end38)
            if reg.empty:
                print(f"[{pheno}][{g}] EMPTY region"); continue
            keep = reg["snp"].isin(ld_snps[g])
            reg = reg[keep].drop_duplicates("snp")
            n_case, n_ctl = PHENOS[pheno]
            reg.to_csv(out, index=False)
            print(f"[{pheno}][{g}] chr{chrom}:{beg38:,}-{end38:,} "
                  f"-> {len(reg)} LD-panel SNPs ({keep.sum()}/{len(keep)}) in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
