"""Rebuild results/tables/susie_coloc_results.csv after the failed-run overwrite.

What happened (found 2026-09-23): a failed `run_susie_coloc.R` overwrote the file with its own
short error-path schema, so the file lost `n_pairs`, `pp_h4_max`, `pp_h4_2nd`, `hit1`, `hit2`,
`idx1`, `idx2`, `niter_eqtl`, `conv_eqtl`, `niter_fg` and `conv_fg`, and 27 of its 44 rows lost
their credible-set counts.

This script restores the file **from `results/tables/susie_coloc_pairs.csv` only** (author
decision 2026-09-23):

  restore   n_pairs, pp_h4_max, pp_h4_2nd, hit1, hit2, idx1, idx2
            (all are exact derivations of the per-credible-set-pair table)
  keep      gene, cohort, n_snps, n_cs_eqtl, n_cs_fg, n_cs_eqtl_full
            (the values the damaged file still holds; every non-null value was checked against
            the pre-overwrite copy and agrees)
  blank     niter_eqtl, conv_eqtl, niter_fg, conv_fg
            (SuSiE convergence diagnostics are NOT derivable from the pair table)
  set       err = "" for the 38 runs that produced pairs, and the production label
            "no CS in one/both datasets" for the six runs whose FinnGen credible-set count is 0
            (production `run_susie_coloc.R` sets err = "" on its success path)

The output uses the exact 18-column schema of the production success path, so a future authentic
re-run of `run_susie_coloc.R` produces a drop-in compatible file. The blanked columns and the
provenance of every restored value are recorded in
`results/tables/_susie_coloc_results_recovery_note.md`.

The pre-recovery (damaged) file is archived to
`results/reports/_pre_renumber_archive/susie_coloc_results_pre_recovery_20260923.csv`.

Usage:  python analysis/scripts/rebuild_susie_coloc_results.py [--apply]
Dry run (default) prints the planned change and the verification report only.
"""
import hashlib
import os
import shutil
import sys

import pandas as pd

APPLY = "--apply" in sys.argv
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TARGET = os.path.join(ROOT, "results", "tables", "susie_coloc_results.csv")
PAIRS = os.path.join(ROOT, "results", "tables", "susie_coloc_pairs.csv")
RETIRED = os.path.join(ROOT, "results", "tables", "supplementary", "_retired",
                       "coloc_abf_susie_44runs.csv")
ARCHIVE = os.path.join(ROOT, "results", "reports", "_pre_renumber_archive",
                       "susie_coloc_results_pre_recovery_20260923.csv")
NOTE = os.path.join(ROOT, "results", "tables", "_susie_coloc_results_recovery_note.md")

KEY = ["gene", "cohort"]
# production schema of the success path in run_susie_coloc.R (lines 117-124)
SCHEMA = ["gene", "cohort", "n_snps", "n_cs_eqtl", "n_cs_fg", "n_cs_eqtl_full",
          "n_pairs", "pp_h4_max", "pp_h4_2nd", "hit1", "hit2", "idx1", "idx2",
          "niter_eqtl", "conv_eqtl", "niter_fg", "conv_fg", "err"]
KEEP = ["n_snps", "n_cs_eqtl", "n_cs_fg", "n_cs_eqtl_full"]
BLANK = ["niter_eqtl", "conv_eqtl", "niter_fg", "conv_fg"]
NO_CS = "no CS in one/both datasets"


def md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest().upper()[:16]


def build():
    cur = pd.read_csv(TARGET)
    pairs = pd.read_csv(PAIRS)
    retired = pd.read_csv(RETIRED)
    old_keys = list(map(tuple, cur[KEY].values))

    # ---- derive the pair-based columns -------------------------------------
    df = cur[KEY + KEEP].copy()
    grp = pairs.groupby(KEY, sort=False)
    derived = []
    for k, sub in grp:
        s = sub["PP.H4.abf"].sort_values(ascending=False).reset_index(drop=True)
        best = sub.loc[sub["PP.H4.abf"].idxmax()]
        derived.append({
            "gene": k[0], "cohort": k[1], "n_pairs": len(sub),
            "pp_h4_max": round(float(s.iloc[0]), 4),
            "pp_h4_2nd": round(float(s.iloc[1]), 4) if len(s) > 1 else float("nan"),
            "hit1": best["hit1"], "hit2": best["hit2"],
            "idx1": best["idx1"], "idx2": best["idx2"],
        })
    df = df.merge(pd.DataFrame(derived), on=KEY, how="left")
    # the six runs with no credible-set pair stay empty (the pre-overwrite file, and the
    # production failure path, leave n_pairs unset rather than writing 0)
    has_pairs = df["n_pairs"].notna()

    # ---- err: "" where a pair table exists, production label otherwise -----
    df["err"] = ""
    df.loc[~has_pairs, "err"] = NO_CS

    # ---- blank the non-derivable convergence diagnostics -------------------
    for c in BLANK:
        df[c] = pd.NA

    df = df[SCHEMA]
    return cur, pairs, retired, df, old_keys


def verify(cur, pairs, retired, new, old_keys):
    ok = True
    print(f"\nrows: {len(new)} (expected 44) | schema: {len(new.columns)} columns")
    if len(new) != 44:
        ok = False
        print("  FAIL row count")
    if list(new.columns) != SCHEMA:
        ok = False
        print("  FAIL schema differs from the production success path")
    if list(map(tuple, new[KEY].values)) != old_keys or len(set(old_keys)) != 44:
        ok = False
        print("  FAIL gene/cohort key set or order changed")
    else:
        print("  OK  gene x cohort keys unchanged (44 unique, same order)")

    # preserved columns must be untouched
    kept = cur[KEY + KEEP].merge(new[KEY + KEEP], on=KEY, suffixes=("_old", "_new"))
    for c in KEEP:
        same = kept[f"{c}_old"].fillna(-999) == kept[f"{c}_new"].fillna(-999)
        if not same.all():
            ok = False
            print(f"  FAIL preserved column {c} changed in {(~same).sum()} row(s)")
    print(f"  OK  preserved columns untouched ({', '.join(KEEP)}); "
          f"non-null after recovery: " +
          ", ".join(f"{c}={int(new[c].notna().sum())}/44" for c in KEEP))

    # independent recomputation of the pair-derived columns
    grp = pairs.groupby(KEY, sort=False)
    bad = 0
    for _, r in new.iterrows():
        k = (r["gene"], r["cohort"])
        if k in grp.groups:
            s = grp.get_group(k)["PP.H4.abf"].sort_values(ascending=False)
            if abs(r["pp_h4_max"] - round(float(s.iloc[0]), 4)) > 1e-12 or r["n_pairs"] != len(s):
                bad += 1
        elif not (pd.isna(r["pp_h4_max"]) and pd.isna(r["n_pairs"])):
            bad += 1
    if bad:
        ok = False
        print(f"  FAIL pair-derived values wrong in {bad} row(s)")
    else:
        print("  OK  n_pairs / pp_h4_max reproduce independently from susie_coloc_pairs.csv")

    # cross-check against the pre-overwrite copy
    m = new.merge(retired[KEY + ["n_pairs", "susie_PP.H4_max"]], on=KEY, how="left",
                  suffixes=("", "_ret"))
    d1 = (m["pp_h4_max"].round(4).fillna(-1) - m["susie_PP.H4_max"].round(4).fillna(-1)).abs()
    d2 = (m["n_pairs"].fillna(-1) - m["n_pairs_ret"].fillna(-1)).abs()
    if d1.max() > 1e-9 or d2.max() > 1e-9:
        ok = False
        print(f"  FAIL vs pre-overwrite copy: pp_h4_max max diff {d1.max():.2e}, "
              f"n_pairs max diff {d2.max():.2e}")
    else:
        print("  OK  pp_h4_max (4 dp) and n_pairs match the pre-overwrite copy in 44/44 rows")
    print(f"  non-estimable runs (pp_h4_max NA, n_pairs NA): "
          f"{int(new.pp_h4_max.isna().sum())}")
    print(f"  convergence columns blank: {', '.join(BLANK)}")
    return ok


def note_text(new, pairs, cur, retired):
    gain = {}
    m = new.merge(retired, on=KEY, how="left", suffixes=("", "_ret"))
    for c in ["n_cs_eqtl", "n_cs_fg", "n_cs_eqtl_full", "n_snps", "conv_eqtl", "conv_fg"]:
        gain[c] = int((m[c].isna() & m[f"{c}_ret"].notna()).sum())
    return f"""# Recovery note — `results/tables/susie_coloc_results.csv`

Written 2026-09-23 by `analysis/scripts/rebuild_susie_coloc_results.py` (author decision:
rebuild the 44 rows **from `susie_coloc_pairs.csv` only**).

## Why

A failed `run_susie_coloc.R` run overwrote this file with its own short error-path schema
(7 columns: `gene, cohort, err, n_snps, n_cs_fg, n_cs_eqtl, n_cs_eqtl_full`). Compared with the
production success-path schema (18 columns) it had lost `n_pairs`, `pp_h4_max`, `pp_h4_2nd`,
`hit1`, `hit2`, `idx1`, `idx2`, `niter_eqtl`, `conv_eqtl`, `niter_fg`, `conv_fg`, and 27 of the
44 rows had lost their credible-set counts. The damaged state is preserved at
`results/reports/_pre_renumber_archive/susie_coloc_results_pre_recovery_20260923.csv`.

## What each column now holds

| Column(s) | Provenance |
|---|---|
| `gene`, `cohort` | unchanged keys (22 WIDE Bonferroni hits × 2 OLP definitions) |
| `n_pairs`, `pp_h4_max`, `pp_h4_2nd`, `hit1`, `hit2`, `idx1`, `idx2` | **re-derived from `susie_coloc_pairs.csv`** (exact functions of the credible-set-pair table: pair count, largest / second-largest `PP.H4.abf`, and the identifiers of the best pair). `pp_h4_max` is rounded to 4 dp to match `run_susie_coloc.R` |
| `n_snps`, `n_cs_eqtl`, `n_cs_fg`, `n_cs_eqtl_full` | **kept** from the damaged file; every non-null value was checked against the pre-overwrite copy and agrees. Rows whose counts the failed run lost are left empty (`n_snps` present in 17/44) |
| `err` | `""` for the 38 runs that produced a pair table (the production success path sets `err = ""`); `"no CS in one/both datasets"` for the 6 runs whose FinnGen credible-set count is 0 — the production label for that condition |
| `niter_eqtl`, `conv_eqtl`, `niter_fg`, `conv_fg` | **blank / unavailable** — SuSiE iteration and convergence diagnostics are not derivable from the pair table |

Rows with `pp_h4_max` empty and `n_pairs` empty: **{int(new.pp_h4_max.isna().sum())}** — the
"no FinnGen credible set / run not estimable" cases, the same six runs the figure legends call
non-estimable. Those rows carry no pair table at all, so `n_pairs` is left empty rather than 0,
exactly as the pre-overwrite file (and the production failure path) had it.

## Still recoverable, not restored here

`results/tables/supplementary/_retired/coloc_abf_susie_44runs.csv` is a pre-overwrite derived
copy of this file (built as the former Supplementary Table 4) and still holds, for all 44 rows,
values that are blank or partial above. Rows that would gain data from it:

| Column | rows still blank here | fillable from the retired copy |
|---|---|---|
| `n_cs_eqtl` | {int(new.n_cs_eqtl.isna().sum())} | {gain['n_cs_eqtl']} |
| `n_cs_fg` | {int(new.n_cs_fg.isna().sum())} | {gain['n_cs_fg']} |
| `n_cs_eqtl_full` | {int(new.n_cs_eqtl_full.isna().sum())} | {gain['n_cs_eqtl_full']} |
| `n_snps` | {int(new.n_snps.isna().sum())} | {gain['n_snps']} |
| `conv_eqtl` | {int(new.conv_eqtl.isna().sum())} | {gain['conv_eqtl']} |
| `conv_fg` | {int(new.conv_fg.isna().sum())} | {gain['conv_fg']} |

They were deliberately **not** copied back here, because the author decision restricted this
rebuild to the pair table. If a full restoration is wanted, take those columns from that retired
copy (or re-run `run_susie_coloc.R`, which regenerates every column authentically).

The SuSiE convergence diagnostics are still not derivable from this file. To stop a
regeneration from silently blanking them, `analysis/scripts/build_supplementary_tables.py`
**carries those values forward** from the retired copy (printing a
`carried forward from the retired copy` warning) instead of overwriting them with nulls. With
that in place the table pipeline runs end to end again and its ten outputs are byte-identical to
the files on disk — verified by `analysis/scripts/_test_supp_table_rebuild.py` (10/10).

## Downstream scripts unblocked by this rebuild

`plot_susie_coloc.R` (Figure 4B), `journal_fig4.R`, `qa_susie_consistency.R`,
`qa_susie_pleiotropy.R`, `run_susie_pleiotropy.R` and `build_supplementary_tables.py` all read
`pp_h4_max`, which now exists again — Figure 4B can be re-rendered from this file once more.
(This rebuild has **not** re-rendered any figure; the existing figure files are unchanged.)

Independent confirmation that the restored values are the pre-overwrite ones:
`analysis/scripts/qa_susie_consistency.R` cross-checks `pp_h4_max` against
`results/tables/susie_abf_comparison.csv` and reproduces all 44 verdicts, including the values
quoted in the manuscript — CD44 0.0302, TDP2 0.1768, EGFL8 (INCLAVO) 0.9585, TAPBP 0.9879/0.9827
and HLA-DPA1 1.0000 / 0.9999.
"""


def main():
    cur, pairs, retired, new, old_keys = build()
    print("=" * 78)
    print("RECOVERY PLAN (from susie_coloc_pairs.csv only)")
    print("=" * 78)
    print(f"target            : results/tables/susie_coloc_results.csv")
    print(f"  before          : {cur.shape[1]} columns, MD5 {md5(TARGET)}")
    print(f"  after           : {len(SCHEMA)} columns (production success-path schema)")
    print(f"  restored        : n_pairs, pp_h4_max, pp_h4_2nd, hit1, hit2, idx1, idx2")
    print(f"  blanked         : {', '.join(BLANK)}")
    print(f"  archive damaged : results/reports/_pre_renumber_archive/"
          f"susie_coloc_results_pre_recovery_20260923.csv")

    ok = verify(cur, pairs, retired, new, old_keys)

    if not APPLY:
        print("\nDRY-RUN — nothing written (pass --apply)")
        sys.exit(0 if ok else 1)
    if not ok:
        print("\nVERIFICATION FAILED — nothing written")
        sys.exit(1)

    if not os.path.exists(ARCHIVE):
        shutil.copy2(TARGET, ARCHIVE)
        print(f"\narchived damaged file -> {os.path.relpath(ARCHIVE, ROOT)}")
    new.to_csv(TARGET, index=False)
    with open(NOTE, "w", encoding="utf-8", newline="") as f:
        f.write(note_text(new, pairs, cur, retired))
    print(f"wrote results/tables/susie_coloc_results.csv ({len(new)} rows, "
          f"MD5 {md5(TARGET)})")
    print("wrote results/tables/_susie_coloc_results_recovery_note.md")


if __name__ == "__main__":
    main()
