"""analyse_decoys.py — does the docking actually separate candidates from non-binders?

The decoy set is property-matched (MW/logP/HBD/HBA/RB/TPSA inside each candidate's windows) and
topologically dissimilar (max Morgan-2 Tanimoto < 0.30 against every candidate), so it is a fair
stand-in for molecules not expected to bind. If the candidate scores were not better than the
decoys', then an absolute affinity alone would carry no evidence of binding.

Reports, per receptor and overall:
  * the pass rate at the text's -7 kcal/mol threshold, candidates vs decoys
  * a two-sided Fisher exact test on that 2x2 table
  * the rank-based AUC (P[a random candidate scores better than a random decoy]) with Mann-Whitney U
  * the rank of the best candidate among all molecules docked, and the decoy best score

Writes results/tables/docking_decoy_control.csv and prints the summary.

Usage:
    python analysis/scripts/analyse_decoys.py
"""
import csv
import json
import os
from collections import defaultdict
from statistics import median

from scipy.stats import fisher_exact, mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REDO = os.path.join(ROOT, "analysis", "docking", "redo")
CAND = os.path.join(REDO, "docking_redo_summary.json")
DEC = os.path.join(REDO, "docking_redo_decoys_summary.json")
OUT = os.path.join(ROOT, "results", "tables", "docking_decoy_control.csv")
THRESH = -7.0


def load(path, is_decoy):
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}")
    out = []
    for r in json.load(open(path, encoding="utf-8")):
        out.append({"target": r["target"], "name": r["ligand"],
                    "score": r["best_affinity"], "kind": "decoy" if is_decoy else "candidate"})
    return out


def main():
    cands = load(CAND, False)
    decoys = load(DEC, True)
    print(f"candidates {len(cands)}, decoys {len(decoys)}\n")

    rows = []
    for scope in ("GDF11", "SUOX", "all"):
        c = [r for r in cands if scope == "all" or r["target"] == scope]
        d = [r for r in decoys if scope == "all" or r["target"] == scope]
        if not c or not d:
            continue
        cp = sum(1 for r in c if r["score"] < THRESH)
        dp = sum(1 for r in d if r["score"] < THRESH)
        table = [[cp, len(c) - cp], [dp, len(d) - dp]]
        odds, p_fish = fisher_exact(table, alternative="two-sided")
        # AUC from U: higher score = better = more negative affinity, so use -score as the ranking
        u, p_mw = mannwhitneyu([-r["score"] for r in c], [-r["score"] for r in d],
                               alternative="greater")
        auc = u / (len(c) * len(d))
        allm = sorted(c + d, key=lambda r: r["score"])
        best_rank = next(i for i, r in enumerate(allm, 1) if r["kind"] == "candidate")
        rows.append({"scope": scope, "n_candidates": len(c), "n_decoys": len(d),
                     "cand_pass_-7": cp, "cand_pass_pct": round(100 * cp / len(c), 1),
                     "decoy_pass_-7": dp, "decoy_pass_pct": round(100 * dp / len(d), 1),
                     "cand_median": round(median(r["score"] for r in c), 3),
                     "decoy_median": round(median(r["score"] for r in d), 3),
                     "best_cand": round(min(r["score"] for r in c), 3),
                     "best_decoy": round(min(r["score"] for r in d), 3),
                     "fisher_p": round(p_fish, 4), "auc": round(auc, 3),
                     "auc_p": round(p_mw, 4), "best_cand_rank_of_all": best_rank,
                     "n_all": len(allm)})

    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        print(f"[{r['scope']}]  candidates {r['n_candidates']}, decoys {r['n_decoys']}")
        print(f"   pass < {THRESH}: candidates {r['cand_pass_-7']} "
              f"({r['cand_pass_pct']}%) vs decoys {r['decoy_pass_-7']} ({r['decoy_pass_pct']}%)")
        print(f"   median: candidates {r['cand_median']} vs decoys {r['decoy_median']}")
        print(f"   best:   candidates {r['best_cand']} vs decoys {r['best_decoy']}")
        print(f"   Fisher exact p = {r['fisher_p']}; AUC = {r['auc']} (p = {r['auc_p']})")
        print(f"   best candidate ranks {r['best_cand_rank_of_all']} of {r['n_all']} molecules")
        print()

    # interpretation, stated against fixed thresholds so it does not drift
    overall = [r for r in rows if r["scope"] == "all"][0]
    print("INTERPRETATION")
    print(f"  AUC {overall['auc']} - 0.5 is no discrimination, 1.0 is perfect separation.")
    if overall["auc"] < 0.6:
        print("  -> the protocol barely separates candidates from property-matched non-binders, so an")
        print("     absolute affinity is NOT evidence of binding for this target set.")
    elif overall["auc"] < 0.7:
        print("  -> weak separation only; docking should be reported as a coarse screen.")
    else:
        print("  -> the candidates rank above the decoys; docking carries some discriminating power.")
    if overall["fisher_p"] > 0.05:
        print("  -> the difference in pass rate is not significant (p > 0.05).")
    print(f"\ntable -> {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
