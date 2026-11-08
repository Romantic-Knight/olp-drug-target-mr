# -*- coding: utf-8 -*-
"""
GSE211630: read 6 count matrices -> per-sample AnnData -> cell-level QC filtering
(mt%<20% for the 3 MT-containing samples; MAD outlier + scrublet doublet removal
for the 3 samples without MT rows) -> per-sample filtered AnnData + per-sample
10x MTX export for Seurat + cell-level QC metrics.

Usage (project root):
    <scrna pixi python> analysis/scripts/scrna_preprocess_gse211630.py

Outputs (project-relative):
  data/processed/geo_GSE211630/ann/<GSM>_filtered.h5ad     (per-sample filtered AnnData)
  data/processed/geo_GSE211630/seurat_export/<GSM>/...   (matrix.mtx, features.tsv, barcodes.tsv, metadata.csv)
  data/processed/geo_GSE211630/qc_cell_metrics.csv.gz
  data/processed/geo_GSE211630/preprocess_summary.json

Rationale for per-sample objects (no full union matrix here): sample gene
universes differ strongly (normal-1 has only ~12k of ~27k genes; 3 samples lack
MT-* rows). A union-wide merged matrix would materialise large explicit-zero
blocks for normal-1 and is exactly the filtering-inconsistency artefact we
documented in the first-round QC report. Cross-sample integration is deferred to
the clustering stage on per-sample common highly variable genes.

File-format note (important): the GEO count-matrix header row consists ONLY of
cell barcodes (no leading gene-column placeholder), while data rows are
"<gene>\t<count per cell>". Cell count must therefore be taken from the data
rows (header may be either barcodes-only or barcodes-prefixed-by-a-gene cell);
the loader calibrates against the first data row.
"""
import gzip
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Anchor all IO to the project root (this file lives in <root>/analysis/scripts,
# three levels up from the file), so the pipeline does not depend on the
# invoking working directory.
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)

MATRIX_DIR = os.path.join("data", "raw", "geo", "GSE211630", "matrix")

# GSM -> (file, condition, label)
SAMPLES = [
    ("GSM6481637", "GSM6481637_count_data-normal-1.txt.gz", "normal", "normal-1"),
    ("GSM6481638", "GSM6481638_count_data-NEOLP-1.txt.gz", "NEOLP", "NEOLP-1"),
    ("GSM6481639", "GSM6481639_count_data-NEOLP-2.txt.gz", "NEOLP", "NEOLP-2"),
    ("GSM6481640", "GSM6481640_count_data-NEOLP-3.txt.gz", "NEOLP", "NEOLP-3"),
    ("GSM6481641", "GSM6481641_count_data-EOLP-1.txt.gz", "EOLP", "EOLP-1"),
    ("GSM6481642", "GSM6481642_count_data-EOLP-2.txt.gz", "EOLP", "EOLP-2"),
]
MT_SAMPLES = {"GSM6481639", "GSM6481640", "GSM6481641"}  # samples carrying MT-* rows
MT_FRAC_CUT = 0.20
MT_FRAC_REPORT = 0.10
MAD_K = 3.0
DOUBLET_RATE = 0.06

OUTDIR = os.path.join("data", "processed", "geo_GSE211630")


def _load(path):
    """Two-pass parser -> scipy CSR in cells x genes orientation (no transpose).

    Pass 1 collects gene names and counts nonzeros, calibrating the cell count
    from the first data row. Pass 2 fills COO arrays; COO -> CSR conversion
    (stable path; avoids the csc->csr transpose that corrupted earlier runs).
    """
    t0 = time.time()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        genes = []
        nnz = 0
        expected = None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if expected is None:
                expected = len(parts) - 1
            genes.append(parts[0].strip('"'))
            for tok in parts[1:]:
                if tok != "0" and tok:
                    nnz += 1
    if expected is None:
        raise ValueError(f"empty matrix file: {path}")
    if len(header) == expected:
        barcodes = [t.strip('"') for t in header]
    elif len(header) == expected + 1:
        barcodes = [t.strip('"') for t in header[1:]]
    else:
        raise ValueError(f"header cells {len(header)} vs data cells {expected}: {path}")
    if len(barcodes) != expected:
        raise ValueError(f"barcode count {len(barcodes)} != data cells {expected}: {path}")
    if len(genes) == 0:
        raise ValueError(f"no gene rows: {path}")

    row = np.empty(nnz, dtype=np.int32)
    colg = np.empty(nnz, dtype=np.int32)
    data = np.empty(nnz, dtype=np.int32)
    pos = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()
        for gi, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            for ci in range(1, len(parts)):
                tok = parts[ci]
                if tok == "0" or not tok:
                    continue
                row[pos] = ci - 1
                colg[pos] = gi
                data[pos] = int(tok)
                pos += 1
    if pos != nnz:
        raise ValueError(f"nnz mismatch {pos} != {nnz}: {path}")

    coo = sp.coo_matrix((data, (row, colg)), shape=(expected, len(genes)))
    csr = coo.tocsr()
    del row, colg, data, coo
    print(f"  loaded {os.path.basename(path)}: cells x genes {csr.shape} "
          f"nnz={nnz} in {time.time()-t0:.0f}s", flush=True)
    return csr, barcodes, genes


def mad_limits(vals_log10):
    med = np.median(vals_log10)
    mad = 1.4826 * np.median(np.abs(vals_log10 - med))
    return med - MAD_K * mad, med + MAD_K * mad


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    summary = {}
    cell_rows = []
    adata_parts = []

    import anndata as ad

    for gsm, fname, condition, label in SAMPLES:
        t0 = time.time()
        print(f"\n=== {gsm} ({label}, {condition}) ===", flush=True)
        csr, barcodes, genes = _load(os.path.join(MATRIX_DIR, fname))
        X = csr.astype(np.float32)
        del csr
        if len(genes) != X.shape[1]:
            raise ValueError(f"gene count mismatch {len(genes)} != {X.shape[1]}: {fname}")
        adata = ad.AnnData(X, obs=pd.DataFrame(index=barcodes))
        adata.var_names = genes
        adata.obs["sample"] = gsm
        adata.obs["condition"] = condition
        adata.obs["label"] = label

        n_counts = np.asarray(X.sum(axis=1)).ravel()
        n_genes = np.diff(X.indptr)
        mt_idx = np.array([i for i, g in enumerate(genes) if g.startswith("MT-")])
        mt_frac = np.zeros(X.shape[0], dtype=np.float64)
        if len(mt_idx):
            mt_umi = np.asarray(X[:, mt_idx].sum(axis=1)).ravel()
            with np.errstate(divide="ignore", invalid="ignore"):
                mt_frac = np.where(n_counts > 0, mt_umi / np.maximum(n_counts, 1), 0.0)
        else:
            mt_frac[:] = np.nan
        adata.obs["n_counts"] = n_counts
        adata.obs["n_genes"] = n_genes
        adata.obs["mt_frac"] = mt_frac

        doublet_score = np.full(X.shape[0], np.nan)
        predicted_doublet = np.zeros(X.shape[0], dtype=bool)
        threshold = None
        if gsm not in MT_SAMPLES:  # scrublet doublet prediction only for MT-free samples
            try:
                from scrublet import Scrublet
                scr = Scrublet(X, expected_doublet_rate=DOUBLET_RATE)
                ds, pr = scr.scrub_doublets(verbose=False)
                doublet_score = np.asarray(ds)
                predicted_doublet = np.asarray(pr).astype(bool)
                threshold = float(scr.threshold_)
                print(f"  scrublet done: threshold={threshold:.3f} "
                      f"predicted={int(predicted_doublet.sum())}/{len(predicted_doublet)}",
                      flush=True)
            except Exception as e:
                print(f"  scrublet FAILED: {e!r}", flush=True)
        else:
            print("  scrublet skipped (MT-sample filter set: mt%<cut only)", flush=True)
        adata.obs["doublet_score"] = doublet_score
        adata.obs["predicted_doublet"] = predicted_doublet

        kept_mt = np.ones(X.shape[0], dtype=bool)
        kept_mad = np.ones(X.shape[0], dtype=bool)
        kept_doublet = np.ones(X.shape[0], dtype=bool)
        n_mt10 = n_mt20 = n_mad_lo = n_mad_hi = n_dbl = 0

        if gsm in MT_SAMPLES:
            kept_mt = mt_frac < MT_FRAC_CUT
            n_mt20 = int((~kept_mt).sum())
            n_mt10 = int((mt_frac >= MT_FRAC_REPORT).sum())
        else:
            lg_nc = np.log10(np.maximum(n_counts, 1.0))
            lg_ng = np.log10(np.maximum(n_genes, 1.0))
            lo_c, hi_c = mad_limits(lg_nc)
            lo_g, hi_g = mad_limits(lg_ng)
            mad_lo = (lg_nc < lo_c) | (lg_ng < lo_g)
            mad_hi = (lg_nc > hi_c) | (lg_ng > hi_g)
            kept_mad = ~(mad_lo | mad_hi)
            n_mad_lo = int(mad_lo.sum())
            n_mad_hi = int(mad_hi.sum())
            kept_doublet = ~predicted_doublet
            n_dbl = int(predicted_doublet.sum())

        keep = kept_mt & kept_mad & kept_doublet
        n_keep = int(keep.sum())

        print(f"  filter: n_cells={X.shape[0]} -> kept={n_keep} "
              f"(mt>20%:{n_mt20} mt>=10%:{n_mt10} mad_lo:{n_mad_lo} "
              f"mad_hi:{n_mad_hi} dbl:{n_dbl})", flush=True)

        adata.obs["kept_mt"] = kept_mt
        adata.obs["kept_mad"] = kept_mad
        adata.obs["kept_doublet"] = kept_doublet
        adata.obs["kept_final"] = keep
        adata.obs["scrublet_threshold"] = threshold if threshold is not None else np.nan

        for i in range(X.shape[0]):
            cell_rows.append({
                "cell": barcodes[i], "sample": gsm, "condition": condition,
                "label": label,
                "n_genes": int(n_genes[i]), "n_counts": int(n_counts[i]),
                "mt_frac": (None if np.isnan(mt_frac[i]) else float(mt_frac[i])),
                "doublet_score": (None if np.isnan(doublet_score[i]) else float(doublet_score[i])),
                "predicted_doublet": bool(predicted_doublet[i]),
                "kept_mt": bool(kept_mt[i]), "kept_mad": bool(kept_mad[i]),
                "kept_doublet": bool(kept_doublet[i]), "kept_final": bool(keep[i]),
            })

        summary[gsm] = {
            "label": label, "condition": condition,
            "genes_raw": int(len(genes)),
            "cells_raw": int(X.shape[0]),
            "cells_mt20": n_mt20, "cells_mt10": n_mt10,
            "cells_mad_lo": n_mad_lo, "cells_mad_hi": n_mad_hi,
            "cells_doublet": n_dbl,
            "doublet_threshold": threshold,
            "cells_kept": n_keep,
            "pct_kept": round(100.0 * n_keep / X.shape[0], 2),
        }

        adata_parts.append(adata[keep].copy())
        print(f"  [{label}] {time.time()-t0:.0f}s elapsed", flush=True)

    ann_dir = os.path.join(OUTDIR, "ann")
    os.makedirs(ann_dir, exist_ok=True)
    for adata in adata_parts:
        gsm = adata.obs["sample"].iloc[0]
        ap = os.path.join(ann_dir, f"{gsm}_filtered.h5ad")
        adata.write_h5ad(ap, compression="gzip")
        print(f"  wrote {ap} ({adata.shape}, "
              f"{os.path.getsize(ap)/1e6:.1f} MB)", flush=True)

    # Merged object restricted to the common-gene universe (join="inner"): no
    # explicit-zero blocks are materialised for samples with narrower gene
    # universes (normal-1 ~12k of ~27k). This is the object downstream
    # clustering/integration should consume; per-sample h5ad above remain the
    # single-sample-validated originals.
    merged = ad.concat(adata_parts, join="inner")
    m_path = os.path.join(OUTDIR, "gse211630_filtered.h5ad")
    merged.write_h5ad(m_path, compression="gzip")
    print(f"  wrote {m_path} ({merged.shape}, "
          f"{os.path.getsize(m_path)/1e6:.1f} MB, common genes={merged.n_vars})",
          flush=True)

    export_dir = os.path.join(OUTDIR, "seurat_export")
    os.makedirs(export_dir, exist_ok=True)
    for adata in adata_parts:
        gsm = adata.obs["sample"].iloc[0]
        d = os.path.join(export_dir, gsm)
        os.makedirs(d, exist_ok=True)
        Xc = adata.X
        if sp.issparse(Xc):
            Xc = Xc.tocoo()
        else:
            Xc = sp.coo_matrix(Xc)
        with open(os.path.join(d, "matrix.mtx"), "w") as f:
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write(f"{Xc.shape[1]} {Xc.shape[0]} {Xc.nnz}\n")
            for r, c, v in zip(Xc.row, Xc.col, Xc.data):
                f.write(f"{c+1} {r+1} {int(v)}\n")
        with open(os.path.join(d, "features.tsv"), "w") as f:
            for g in adata.var_names:
                f.write(f"{g}\t{g}\tGene Expression\n")
        with open(os.path.join(d, "barcodes.tsv"), "w") as f:
            for b in adata.obs_names:
                f.write(f"{b}\n")
        adata.obs.to_csv(os.path.join(d, "metadata.csv"))
        print(f"  MTX export {gsm}: {adata.shape}", flush=True)

    df = pd.DataFrame(cell_rows)
    metrics_path = os.path.join(OUTDIR, "qc_cell_metrics.csv.gz")
    df.to_csv(metrics_path, index=False, compression="gzip")
    print(f"  wrote {metrics_path} ({df.shape[0]} rows)", flush=True)

    sum_path = os.path.join(OUTDIR, "preprocess_summary.json")
    totals = {
        "total_cells_raw": int(sum(v["cells_raw"] for v in summary.values())),
        "total_cells_kept": int(sum(v["cells_kept"] for v in summary.values())),
        "genes_common_merged": int(merged.n_vars),
    }
    with open(sum_path, "w", encoding="utf-8") as f:
        json.dump({"samples": summary, "totals": totals, "params": {
            "mt_frac_cut": MT_FRAC_CUT, "mt_frac_report": MT_FRAC_REPORT,
            "mad_k": MAD_K, "doublet_rate": DOUBLET_RATE,
            "mt_samples": sorted(MT_SAMPLES),
        }}, f, indent=2)
    print(f"  wrote {sum_path}", flush=True)
    print("\nALL DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
