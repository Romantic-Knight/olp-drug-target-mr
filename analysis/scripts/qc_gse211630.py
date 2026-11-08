# -*- coding: utf-8 -*-
"""QC of GSE211630 per-sample count matrices.

For each sample matrix (gene x cell, gzipped TSV):
  - dimensions (n_genes, n_cells from header)
  - total UMI / cells; UMI per cell; genes detected per cell (approx via row scan)
  - gene symbol sanity: fraction of human-style symbols vs mouse-style symbols;
    presence of human marker genes (HLA-DRA, CD3D, MS4A1, GDF11, SUOX, TAPBP,
    EGFL8, KRT14, EPCAM...) ; presence of mouse-only patterns (GmNNNN, mt-*, 
    lowerCamel mouse symbols e.g. Actb) ; mitochondrial prefix MT- vs mt-
Key question: does the matrix reflect human (GRCh38) or mouse (mm10) annotation?
Output: data/processed/geo_GSE211630/qc_report.json (+ console summary)
"""
import gzip
import json
import os
import re

SRC = os.path.join("data", "raw", "geo", "GSE211630", "matrix")
OUT = os.path.join("data", "processed", "geo_GSE211630")
os.makedirs(OUT, exist_ok=True)

HUMAN_MARKERS = [
    "HLA-A", "HLA-B", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "TAPBP", "CD3D", "CD3E",
    "CD3G", "CD4", "CD8A", "MS4A1", "CD79A", "CD14", "LYZ", "NKG7", "GZMB", "GDF11",
    "SUOX", "EGFL8", "KRT5", "KRT14", "EPCAM", "PTPRC", "VIM", "COL1A1", "MUC1",
]
MOUSE_ONLY_PATTERNS = [
    re.compile(r"^Gm\d+$"),          # mouse predicted genes
    re.compile(r"^mt-[A-Z]"),        # mouse mt genes lowercase 'mt-'
    re.compile(r"^[A-Z][a-z]{2,}\d{1,2}[a-z]$"),  # mouse CamelCase symbol like Actb? (Actb no digit)
]
MOUSE_CAMEL = re.compile(r"^[A-Z][a-z]+[a-z0-9]*$")  # e.g. Actb, Gapdh, Cd3d
HUMAN_MT = re.compile(r"^MT-")
HUMAN_LINC = re.compile(r"^LINC\d+")


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def qc_file(path):
    res = {}
    per_cell = []            # cell barcodes
    col_umi = []             # per-cell total UMI (int)
    col_genes = []           # per-cell gene count (approx: nonzero entries)
    n_genes = 0
    detected_markers = {}
    mouse_camel_hits = 0
    gm_hits = 0
    mt_prefix = {"human": 0, "mouse": 0}
    human_style = 0          # all-caps symbols with digits/dashes (human-like)
    ambiguous = 0
    total_umi = 0
    nz_cells_total = 0
    rows_nonzero = []        # per gene: #cells expressing (nonzero count)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        per_cell = [h.strip('"') for h in header if h.strip() and not h.startswith('"')]
        # header may be quoted; unquote all
        per_cell = [h.strip('"').strip() for h in header]
        per_cell = [h for h in per_cell if h]
        col_umi = [0] * len(per_cell)
        col_genes = [0] * len(per_cell)
        mt_umi = [0] * len(per_cell)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0].strip('"').strip()
            if not gene:
                continue
            n_genes += 1
            # count nonzero per cell + UMI
            nz = 0
            mt_row = bool(HUMAN_MT.match(gene))
            for j, tok in enumerate(parts[1:]):
                if j >= len(col_umi):
                    break
                if not tok or tok == "0":
                    continue
                try:
                    v = int(tok)
                except ValueError:
                    continue
                if v > 0:
                    col_umi[j] += v
                    col_genes[j] += 1
                    if mt_row:
                        mt_umi[j] += v
                    nz += 1
            rows_nonzero.append(nz)
            nz_cells_total += nz
            # symbol classification
            if gene in HUMAN_MARKERS:
                detected_markers[gene] = detected_markers.get(gene, 0) + 1
            if MOUSE_ONLY_PATTERNS[0].match(gene):
                gm_hits += 1
            if MOUSE_ONLY_PATTERNS[1].match(gene):
                mt_prefix["mouse"] += 1
            if HUMAN_MT.match(gene):
                mt_prefix["human"] += 1
            if MOUSE_CAMEL.match(gene) and gene not in HUMAN_MARKERS and \
               not HUMAN_LINC.match(gene):
                mouse_camel_hits += 1
            if re.match(r"^[A-Z0-9][A-Z0-9_.-]+$", gene) or HUMAN_LINC.match(gene):
                human_style += 1
            else:
                ambiguous += 1
    res["n_cells"] = len(per_cell)
    res["n_genes"] = n_genes
    res["total_umi"] = sum(col_umi)
    res["mean_umi_per_cell"] = round(sum(col_umi) / max(1, len(per_cell)), 1)
    res["median_umi_per_cell"] = sorted(col_umi)[len(col_umi) // 2] if col_umi else 0
    res["cells_with_umi0"] = sum(1 for u in col_umi if u == 0)
    res["median_genes_per_cell"] = sorted(col_genes)[len(col_genes) // 2] if col_genes else 0
    res["genes_expressed_in_gt0_cells"] = sum(1 for g in rows_nonzero if g > 0)
    res["genes_expressed_in_gt1_cell"] = sum(1 for g in rows_nonzero if g > 1)
    res["markers_detected"] = detected_markers
    res["mouse_gm_genes"] = gm_hits
    res["mouse_mt_genes"] = mt_prefix["mouse"]
    res["human_mt_genes"] = mt_prefix["human"]
    res["mouse_camelcase_hits"] = mouse_camel_hits
    res["human_style_genes"] = human_style
    res["ambiguous_symbols"] = ambiguous
    tot = [u for u in col_umi if u > 0]
    mtf = [mt_umi[j] / col_umi[j] for j in range(len(col_umi)) if col_umi[j] > 0]
    res["cells_with_mt_umi"] = sum(1 for m in mtf if m > 0)
    res["mt_fraction_mean"] = round(sum(mtf) / len(mtf), 4) if mtf else None
    res["mt_fraction_median"] = sorted(mtf)[len(mtf) // 2] if mtf else None
    return res


def main():
    files = sorted(f for f in os.listdir(SRC) if f.endswith(".txt.gz"))
    report = {"samples": {}}
    for fn in files:
        label = fn.replace("_count_data-", "|").replace(".txt.gz", "")
        print("=" * 80)
        print("QC:", fn)
        r = qc_file(os.path.join(SRC, fn))
        report["samples"][fn] = r
        md = [k for k in HUMAN_MARKERS if k in r["markers_detected"]]
        print(f"  cells={r['n_cells']} genes={r['n_genes']} "
              f"totalUMI={r['total_umi']} meanUMI/cell={r['mean_umi_per_cell']}")
        print(f"  zeroUMIcells={r['cells_with_umi0']} humanMT={r['human_mt_genes']} "
              f"mouseMT={r['mouse_mt_genes']} GmNNNN={r['mouse_gm_genes']} "
              f"mouseCamel={r['mouse_camelcase_hits']}")
        print(f"  human-style symbols={r['human_style_genes']} "
              f"ambiguous={r['ambiguous_symbols']}")
        print(f"  human markers found ({len(md)}): {md}")
        if len(md) < 10:
            missing = [k for k in HUMAN_MARKERS if k not in r["markers_detected"]]
            print("  !! few markers; missing e.g.:", missing[:15])
    with open(os.path.join(OUT, "qc_report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)
    print("\nreport ->", os.path.join(OUT, "qc_report.json"))


if __name__ == "__main__":
    main()
