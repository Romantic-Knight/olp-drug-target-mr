"""compare_docking_redo.py — published vs redone docking scores, side by side.

Produces the table the manuscript's docking section has to be rewritten from: per ligand, the
published single-seed affinity against the redone best-of-three-seeds, the difference, and whether
each crosses the good-docking threshold (< -7 kcal/mol) that the text quotes.

Writes results/tables/docking_redo_comparison.csv and prints a summary.

Usage:
    python analysis/scripts/compare_docking_redo.py
"""
import csv
import json
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PUB = os.path.join(ROOT, "results", "tables", "docking_results.csv")
NEW = os.path.join(ROOT, "analysis", "docking", "redo", "docking_redo_summary.json")
OUT = os.path.join(ROOT, "results", "tables", "docking_redo_comparison.csv")
THRESH = -7.0


def main():
    published = {}
    with open(PUB, encoding="utf-8") as f:
        for row in csv.reader(f):
            if len(row) >= 3 and row[0] in ("GDF11", "SUOX"):
                published[(row[0], row[1].strip())] = float(row[2])

    new = json.load(open(NEW, encoding="utf-8"))
    newmap = {(r["target"], r["ligand"]): r for r in new}

    keys = sorted(set(published) | set(newmap))
    rows = []
    for k in keys:
        t, lig = k
        p = published.get(k)
        n = newmap.get(k)
        rows.append({
            "target": t,
            "ligand": lig,
            "published_kcal": "" if p is None else round(p, 3),
            "redo_best_kcal": "" if n is None else n["best_affinity"],
            "redo_worst_seed_kcal": "" if n is None else n["worst_seed_best"],
            "seed_spread_kcal": "" if n is None else n["seed_spread"],
            "delta_kcal": "" if (p is None or n is None) else round(n["best_affinity"] - p, 3),
            "published_pass": "" if p is None else int(p < THRESH),
            "redo_pass": "" if n is None else int(n["best_affinity"] < THRESH),
        })

    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    both = [r for r in rows if r["published_kcal"] != "" and r["redo_best_kcal"] != ""]
    pub_pass = [r for r in both if r["published_pass"] == 1]
    new_pass = [r for r in both if r["redo_pass"] == 1]
    deltas = [r["delta_kcal"] for r in both]
    spreads = [r["seed_spread_kcal"] for r in both]

    print(f"target-ligand pairs compared: {len(both)}")
    print(f"  pass < {THRESH} kcal/mol: published {len(pub_pass)} -> redo {len(new_pass)}")
    print(f"  score change: min {min(deltas):+.2f}, median "
          f"{sorted(deltas)[len(deltas) // 2]:+.2f}, max {max(deltas):+.2f} kcal/mol")
    print(f"  seed spread: max {max(spreads):.2f} kcal/mol")
    print()
    print("still passing after redo:")
    for r in sorted(new_pass, key=lambda x: x["redo_best_kcal"]):
        print(f"  {r['target']:5s} {r['ligand']:18s} published {r['published_kcal']:>7} "
              f"-> redo {r['redo_best_kcal']:>7}  ({r['delta_kcal']:+.2f})")
    lost = [r for r in pub_pass if r["redo_pass"] == 0]
    print(f"\nno longer passing ({len(lost)}):")
    for r in sorted(lost, key=lambda x: x["delta_kcal"]):
        print(f"  {r['target']:5s} {r['ligand']:18s} published {r['published_kcal']:>7} "
              f"-> redo {r['redo_best_kcal']:>7}  ({r['delta_kcal']:+.2f})")
    gained = [r for r in new_pass if r["published_pass"] == 0]
    if gained:
        print(f"\nnewly passing ({len(gained)}):")
        for r in gained:
            print(f"  {r['target']:5s} {r['ligand']:18s} published {r['published_kcal']:>7} "
                  f"-> redo {r['redo_best_kcal']:>7}  ({r['delta_kcal']:+.2f})")
    print(f"\ntable -> {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
