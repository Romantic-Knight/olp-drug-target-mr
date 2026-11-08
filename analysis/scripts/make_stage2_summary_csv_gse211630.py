# -*- coding: utf-8 -*-
"""One-row-per-sample Stage-2 filtering summary table (results/tables/).

Derived from preprocess_summary.json — no recomputation.
"""
import json
import os

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
SRC = os.path.join("data", "processed", "geo_GSE211630", "preprocess_summary.json")
DST = os.path.join("results", "tables", "geo_gse211630_cellqc_summary.csv")

with open(SRC, encoding="utf-8") as f:
    data = json.load(f)

rows = []
for gsm, s in data["samples"].items():
    rows.append({
        "gsm": gsm, "label": s["label"], "condition": s["condition"],
        "genes_raw": s["genes_raw"], "cells_raw": s["cells_raw"],
        "cells_mt20_excluded": s["cells_mt20"], "cells_mt10_reference": s["cells_mt10"],
        "cells_mad_low": s["cells_mad_lo"], "cells_mad_high": s["cells_mad_hi"],
        "cells_scrublet_doublet": s["cells_doublet"],
        "doublet_threshold": s["doublet_threshold"],
        "cells_kept": s["cells_kept"], "pct_kept": round(s["pct_kept"], 2),
    })
df = pd.DataFrame(rows)
tot = data["totals"]
df.loc[len(df)] = {
    "gsm": "TOTAL", "label": "all samples", "condition": "normal+NEOLP+EOLP",
    "genes_raw": "—", "cells_raw": tot["total_cells_raw"],
    "cells_mt20_excluded": "", "cells_mt10_reference": "",
    "cells_mad_low": "", "cells_mad_high": "", "cells_scrublet_doublet": "",
    "doublet_threshold": "", "cells_kept": tot["total_cells_kept"],
    "pct_kept": round(100.0 * tot["total_cells_kept"] / tot["total_cells_raw"], 2),
}
os.makedirs("results/tables", exist_ok=True)
df.to_csv(DST, index=False)
print(f"wrote {DST}\n{df.to_string(index=False)}")
