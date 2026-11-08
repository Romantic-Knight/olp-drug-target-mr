"""Build 1000G Phase 3 (GRCh37, EUR) LD matrices for SuSiE coloc genes.

For each gene: wanted SNP set (min-P < 1e-4 subset from susie_gene_subsets.csv) ->
locate the tabix bins overlapping the region -> fetch the bgzf byte spans via
HTTP Range -> decompress -> extract EUR genotype dosages -> MAF filter -> LD (Pearson r).

Outputs per gene:
  data/processed/coloc/ld_1000g/<gene>_snps.csv   (rsid,pos,ref,alt,maf)
  data/processed/coloc/ld_1000g/<gene>_ld.csv.gz  (N x N matrix, rownames = rsid)

Usage: python build_ld_1000g.py [gene ...]   (default: all genes in susie_gene_subsets.csv)
"""
import struct, gzip, os, sys, io
import numpy as np
import pandas as pd
import urllib.request

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
TBI_DIR = os.path.join(ROOT, "data", "raw", "1000g")
OUT_DIR = os.path.join(ROOT, "data", "processed", "coloc", "ld_1000g")
PANEL = os.path.join(TBI_DIR, "integrated_call_samples_v3.20130502.ALL.panel")
VCF_BASE = ("https://1000genomes.s3.amazonaws.com/release/20130502/"
            "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------- tabix index parsing ----------------
def parse_tbi(path):
    with open(path, "rb") as f:
        d = f.read()
    # tabix index files are gzip-compressed themselves (single member)
    if d[:2] == b"\x1f\x8b":
        d = gzip.decompress(d)
    assert d[:4] == b"TBI\x01", "not a tabix index"
    off = 4
    n_ref = struct.unpack_from("<i", d, off)[0]; off += 4
    off += 24  # format, col_seq, col_beg, col_end, meta, skip
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

# ---------------- bgzf multi-block decompression ----------------
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
def main():
    subset = pd.read_csv(os.path.join(ROOT, "data", "processed", "coloc", "susie_gene_subsets.csv"))
    args = sys.argv[1:]
    genes = args if args else subset["gene"].tolist()

    panel = pd.read_csv(PANEL, sep="\t")
    eur = set(panel.loc[panel["super_pop"] == "EUR", "sample"])

    # per-chromosome: VCF sample order (from locally extracted headers) + tbi
    sample_cache = {}
    tbi_cache = {}
    for c in sorted(subset["chr"].unique()):
        cs = str(int(c))
        sp = pd.read_csv(os.path.join(TBI_DIR, f"chr{cs}_samples.txt"), header=None)[0].tolist()
        sample_cache[cs] = [i for i, s in enumerate(sp) if s in eur]
        tbi_cache[cs] = parse_tbi(os.path.join(TBI_DIR, f"chr{cs}.tbi"))

    for g in genes:
        sub = subset[subset["gene"] == g]
        chr_s = str(int(sub["chr"].iloc[0]))
        want = pd.read_csv(os.path.join(ROOT, "data", "processed", "coloc", "coloc_eqtl_windows.csv"))
        want = want[want["gene"] == g]
        # wanted rsids: from the subset (minP<1e-4) — stored by merging on snp
        # (susie_gene_subsets.csv only has window coords; re-derive the SNP set here)
        fgw = pd.read_csv(os.path.join(ROOT, "data", "processed", "coloc", "coloc_fg_wide.csv"))
        fgi = pd.read_csv(os.path.join(ROOT, "data", "processed", "coloc", "coloc_fg_inclavo.csv"))
        pw = fgw.set_index("snp")["pval"].to_dict()
        pi = fgi.set_index("snp")["pval"].to_dict()
        want["pw"] = want["snp"].map(pw).fillna(1.0)
        want["pi"] = want["snp"].map(pi).fillna(1.0)
        want["minp"] = np.minimum(want["pvalue"], np.minimum(want["pw"], want["pi"]))
        want = want[want["minp"] < 1e-4].sort_values("minp")
        if len(want) > 2000:
            want = want.head(2000)
        if want.empty:
            print(f"[{g}] no wanted SNPs; skip"); continue

        beg = int(want["pos"].min()) - 2000
        end = int(want["pos"].max()) + 2000
        refs = tbi_cache[chr_s]
        ref_name = chr_s if chr_s in refs else ("chr" + chr_s if ("chr" + chr_s) in refs else None)
        if ref_name is None:
            print(f"[{g}] ref '{chr_s}' not in tbi index"); continue
        bins_idx = refs[ref_name]
        chunks = []
        for b in reg2bins(beg, end):
            chunks.extend(bins_idx.get(b, []))
        if not chunks:
            print(f"[{g}] no index chunks; skip"); continue
        # merge in real-file-offset space (voffset >> 16)
        chunks = sorted(chunks, key=lambda c: (c[0] >> 16, c[1] >> 16))
        merged = []
        for c in chunks:
            rs0, rs1 = c[0] >> 16, c[1] >> 16
            if merged and rs0 <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], rs1))
            else:
                merged.append([rs0, rs1])
        url = VCF_BASE.format(chr_s)
        eur_idx = sample_cache[chr_s]
        wanted_rs = set(want["snp"])
        wanted_pos = dict(zip(want["snp"], want["pos"]))

        rows = []
        # chunk offsets are BGZF virtual offsets: real file offset = voffset >> 16
        for (s0, s1) in merged:
            raw = fetch_bytes(url, s0, s1)
            txt = decompress_bgzf(raw)
            if not txt:
                continue
            for line in txt.decode("utf-8", "replace").split("\n"):
                if not line.startswith("#"):
                    f = line.split("\t")
                    if len(f) < 12:
                        continue
                    rid = f[2]
                    if rid == "." or rid not in wanted_rs:
                        continue
                    if "," in f[4]:
                        continue  # biallelic only
                    pos_i = int(f[1])
                    gt_field = [f[9 + i] for i in eur_idx]
                    gt = []
                    bad = False
                    for x in gt_field:
                        a = x.split(":", 1)[0]
                        if a == "." or a == "./.":
                            gt.append(np.nan); continue
                        if "|" in a:
                            p, q = a.split("|")
                        elif "/" in a:
                            p, q = a.split("/")
                        else:
                            bad = True; break
                        try:
                            gt.append(int(p) + int(q))
                        except ValueError:
                            bad = True; break
                    if bad:
                        continue
                    rows.append((rid, pos_i, f[3], f[4], np.array(gt, dtype=float)))
        if len(rows) < 10:
            print(f"[{g}] too few variants ({len(rows)}); skip")
            continue
        rs = np.array([r[0] for r in rows])
        pos = np.array([r[1] for r in rows])
        ref = np.array([r[2] for r in rows])
        alt = np.array([r[3] for r in rows])
        GT = np.vstack([r[4] for r in rows])
        maf = np.nanmean(GT, axis=1) / 2.0
        keep = (~np.isnan(maf)) & (maf > 0.01) & (maf < 0.99)
        rs, pos, ref, alt, maf, GT = rs[keep], pos[keep], ref[keep], alt[keep], maf[keep], GT[keep]
        if GT.shape[0] < 10:
            print(f"[{g}] too few after MAF filter; skip")
            continue
        # impute missing dosages with per-SNP mean so the LD matrix is complete
        for i in range(GT.shape[0]):
            m = np.isnan(GT[i])
            if m.any():
                GT[i, m] = np.nanmean(GT[i])
        # LD: Pearson correlation of dosages
        GTc = GT - np.mean(GT, axis=1, keepdims=True)
        denom = np.sqrt(np.sum(GTc ** 2, axis=1))
        norm = GTc / denom[:, None]
        LD = norm @ norm.T
        np.fill_diagonal(LD, 1.0)
        snps_out = os.path.join(OUT_DIR, f"{g}_snps.csv")
        ld_out = os.path.join(OUT_DIR, f"{g}_ld.csv.gz")
        pd.DataFrame({"rsid": rs, "pos": pos, "ref": ref, "alt": alt, "maf": maf}).to_csv(snps_out, index=False)
        pd.DataFrame(LD, index=rs, columns=rs).to_csv(ld_out, compression="gzip")
        print(f"[{g}] chr{chr_s}: {len(rows)} fetched -> {GT.shape[0]} SNPs, LD {LD.shape}, saved")

if __name__ == "__main__":
    main()
