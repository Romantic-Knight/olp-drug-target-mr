# -*- coding: utf-8 -*-
"""QC summary figure for GSE211630 (per-sample cells/genes/UMI) + summary CSV."""
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REP = json.load(open(os.path.join("data", "processed", "geo_GSE211630", "qc_report.json"),
                     encoding="utf-8"))["samples"]
CAND = json.load(open(os.path.join("data", "processed", "geo_GSE211630",
                                   "candidate_genes_check.json"), encoding="utf-8"))

ORDER = ["GSM6481637_count_data-normal-1.txt.gz",
         "GSM6481638_count_data-NEOLP-1.txt.gz",
         "GSM6481639_count_data-NEOLP-2.txt.gz",
         "GSM6481640_count_data-NEOLP-3.txt.gz",
         "GSM6481641_count_data-EOLP-1.txt.gz",
         "GSM6481642_count_data-EOLP-2.txt.gz"]
LAB = {"GSM6481637_count_data-normal-1.txt.gz": "normal-1",
       "GSM6481638_count_data-NEOLP-1.txt.gz": "NEOLP-1",
       "GSM6481639_count_data-NEOLP-2.txt.gz": "NEOLP-2",
       "GSM6481640_count_data-NEOLP-3.txt.gz": "NEOLP-3",
       "GSM6481641_count_data-EOLP-1.txt.gz": "EOLP-1",
       "GSM6481642_count_data-EOLP-2.txt.gz": "EOLP-2"}
COL = {"normal": "#9aa0a6", "NEOLP": "#4c8bf5", "EOLP": "#e8710a"}

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))
labels = [LAB[f] for f in ORDER]
n_cells = [REP[f]["n_cells"] for f in ORDER]
n_genes = [REP[f]["n_genes"] for f in ORDER]
umi = [REP[f]["mean_umi_per_cell"] for f in ORDER]
grp_col = [COL["normal" if l.startswith("normal") else "NEOLP" if "NEOLP" in l else "EOLP"]
           for l in labels]

ax = axes[0]
ax.bar(labels, n_cells, color=grp_col)
ax.set_title("Cells per sample"); ax.set_ylabel("cells"); ax.tick_params(axis="x", rotation=45)
ax.set_ylim(0, max(n_cells) * 1.12)
for i, v in enumerate(n_cells):
    ax.text(i, v + max(n_cells) * 0.015, f"{v:,}", ha="center", va="bottom", fontsize=8)
ax = axes[1]
ax.bar(labels, n_genes, color=grp_col)
ax.set_title("Genes (rows) per sample"); ax.set_ylabel("genes"); ax.tick_params(axis="x", rotation=45)
ax.set_ylim(0, max(n_genes) * 1.12)
for i, v in enumerate(n_genes):
    ax.text(i, v + max(n_genes) * 0.015, f"{v:,}", ha="center", va="bottom", fontsize=8)
ax = axes[2]
ax.bar(labels, umi, color=grp_col)
ax.set_title("Mean UMI per cell"); ax.set_ylabel("UMI"); ax.tick_params(axis="x", rotation=45)
ax.set_ylim(0, max(umi) * 1.12)
for i, v in enumerate(umi):
    ax.text(i, v + max(umi) * 0.015, f"{v:,.0f}", ha="center", va="bottom", fontsize=8)
from matplotlib.patches import Patch
handles = [Patch(facecolor=COL["normal"], label="normal oral mucosa (n=1)"),
           Patch(facecolor=COL["NEOLP"], label="NEOLP (n=3)"),
           Patch(facecolor=COL["EOLP"], label="EOLP (n=2)")]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.955),
           ncol=3, fontsize=9, frameon=False)
fig.subplots_adjust(top=0.80, bottom=0.13)
fig.suptitle("GSE211630 per-sample QC (scRNA-seq, oral lichen planus)", y=0.99,
             fontsize=13)
os.makedirs("figures", exist_ok=True)
fig.savefig(os.path.join("figures", "geo_gse211630_qc.png"), dpi=220, bbox_inches="tight")
print("figure saved: figures/geo_gse211630_qc.png")

# ---- summary CSV ----
os.makedirs(os.path.join("results", "tables"), exist_ok=True)
with open(os.path.join("results", "tables", "geo_gse211630_qc_summary.csv"), "w",
          encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample_file", "label", "group", "n_cells", "n_genes", "total_umi",
                "mean_umi_per_cell", "zero_umi_cells", "human_mt_rows",
                "mouse_mt_rows", "Gm_genes", "mouse_camel_hits",
                "GDF11_cells", "SUOX_cells", "EGFL8_cells", "TAPBP_cells",
                "HLA-DPA1_cells"])
    for fn in ORDER:
        r = REP[fn]
        c = CAND[fn]["genes"]
        g = "normal" if LAB[fn].startswith("normal") else ("NEOLP" if "NEOLP" in LAB[fn] else "EOLP")
        w.writerow([fn, LAB[fn], g, r["n_cells"], r["n_genes"], r["total_umi"],
                    r["mean_umi_per_cell"], r["cells_with_umi0"],
                    r["human_mt_genes"], r["mouse_mt_genes"], r["mouse_gm_genes"],
                    r["mouse_camelcase_hits"],
                    c.get("GDF11", {}).get("cells_expressing", "row_absent"),
                    c.get("SUOX", {}).get("cells_expressing", "row_absent"),
                    c.get("EGFL8", {}).get("cells_expressing", "row_absent"),
                    c.get("TAPBP", {}).get("cells_expressing", "row_absent"),
                    c.get("HLA-DPA1", {}).get("cells_expressing", "row_absent")])
print("csv saved: results/tables/geo_gse211630_qc_summary.csv")
