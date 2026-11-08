"""Query LOCAL FinnGen R13 bgzip+tabix summary stats for a region (GRCh38).

Reuses the tabix/BGZF machinery from fetch_finngen_region.py but reads the
full local .gz files (data/raw/finngen_r13/) instead of remote HTTP Range.
Writes a CSV per phenotype:

  data/processed/pqtl/fg_<pheno>_<region>.csv
  columns: snp, chrom, pos, ref, alt, pval, beta, sebeta, af_alt

Usage: python query_finngen_local.py <pheno> <chrom> <beg38> <end38> [<tag>]
"""
import os, sys, struct, gzip
import pandas as pd

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
FG_DIR = os.path.join(ROOT, "data", "raw", "finngen_r13")
OUT_DIR = os.path.join(ROOT, "data", "processed", "pqtl")
os.makedirs(OUT_DIR, exist_ok=True)


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


def decompress_bgzf(data, maxout=500_000_000):
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


def query_local(pheno, chrom, beg38, end38):
    gz_path = os.path.join(FG_DIR, f"finngen_R13_{pheno}.gz")
    tbi_path = gz_path + ".tbi"
    assert os.path.exists(gz_path), gz_path
    assert os.path.exists(tbi_path), tbi_path
    refs = parse_tbi(tbi_path)
    if chrom not in refs:
        print(f"[{pheno}] chrom {chrom} not in index"); return pd.DataFrame()
    chunks = []
    for b in reg2bins(beg38, end38):
        chunks.extend(refs[chrom].get(b, []))
    chunks = sorted(chunks, key=lambda c: (c[0] >> 16, c[1] >> 16))
    merged = []
    for c in chunks:
        rs0, rs1 = c[0] >> 16, c[1] >> 16
        if merged and rs0 <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], rs1))
        else:
            merged.append([rs0, rs1])
    rows = []
    with open(gz_path, "rb") as fh:
        for (s0, s1) in merged:
            fh.seek(s0)
            raw = fh.read(s1 - s0 + 1)
            txt = decompress_bgzf(raw)
            if not txt:
                continue
            for line in txt.decode("utf-8", "replace").split("\n"):
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) < 13:
                    continue
                pos = int(f[1])
                if not (beg38 <= pos <= end38):
                    continue
                if f[4] == "." or "," in f[4] or "," in f[3]:
                    continue  # require single rsid and biallelic
                try:
                    pval = float(f[6]); beta = float(f[8]); se = float(f[9]); maf = float(f[10])
                except ValueError:
                    continue
                rows.append((f[4], f[0].lstrip("chr"), pos, f[2], f[3], pval, beta, se, maf))
    df = pd.DataFrame(rows, columns=["snp", "chrom", "pos", "ref", "alt", "pval", "beta", "sebeta", "af_alt"])
    df = df.drop_duplicates("snp")
    return df


def main():
    pheno, chrom, b, e = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    tag = sys.argv[5] if len(sys.argv) > 5 else f"chr{chrom}_{b}_{e}"
    out = os.path.join(OUT_DIR, f"fg_{pheno}_{tag}.csv")
    df = query_local(pheno, chrom, b, e)
    df.to_csv(out, index=False)
    print(f"[{pheno}] {len(df)} variants -> {out}")


if __name__ == "__main__":
    main()
