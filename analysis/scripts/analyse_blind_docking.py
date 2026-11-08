"""analyse_blind_docking.py — where did the blind-scan poses land?

The blind scan gave some ligands a better score than the focused active-site box did, which implies
their best poses are not in the MTE/SO4 site. This locates the poses relative to the reference
sites so the claim can be stated precisely:

  * MTE/SO4 active site  - centre of the published grid box (midpoint of Mo and the catalytic SO4)
  * catalytic Mo         - the metal position in 2A9A chain A

For each ligand the top modes are read from the blind output and the minimum distance from any
ligand atom to those two reference points is reported. A pose within a few Angstrom of the active
site is "at" it; anything tens of Angstrom away is a different surface pocket.

Usage:
    python analysis/scripts/analyse_blind_docking.py
"""
import glob
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BLIND = os.path.join(ROOT, "analysis", "docking", "redo", "blind_SUOX")
OUT = os.path.join(ROOT, "analysis", "docking", "redo", "blind_SUOX_analysis.json")

# the published SUOX grid centre = midpoint of Mo (1.784, 1.336, 11.699) and catalytic SO4
# S (-2.172, 2.256, 8.631); recomputed below for the record
MO = (1.784, 1.336, 11.699)
SO4_S = (-2.172, 2.256, 8.631)
ACTIVE = tuple((MO[i] + SO4_S[i]) / 2.0 for i in range(3))


def read_modes(path, max_modes=3):
    """[(mode, affinity, [(x,y,z), ...])] for the first `max_modes` MODEL blocks."""
    out, cur, aff, mode = [], None, None, 0
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("MODEL"):
            if cur is not None:
                out.append((mode, aff, cur))
                if len(out) >= max_modes:
                    break
            mode += 1
            cur, aff = [], None
        elif ln.startswith("REMARK") and "VINA RESULT" in ln:
            try:
                aff = float(ln.split()[3])
            except (IndexError, ValueError):
                pass
        elif cur is not None and ln.startswith(("ATOM", "HETATM")):
            cur.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    if cur is not None and len(out) < max_modes:
        out.append((mode, aff, cur))
    return out


def main():
    files = sorted(glob.glob(os.path.join(BLIND, "*.pdbqt")))
    if not files:
        print("no blind outputs found - run `redock.py --blind SUOX` first")
        return 1
    print(f"active-site reference (Mo/SO4 midpoint): "
          f"({ACTIVE[0]:.2f}, {ACTIVE[1]:.2f}, {ACTIVE[2]:.2f})\n")
    results = []
    for f in files:
        name = os.path.basename(f)[:-6]
        modes = read_modes(f)
        if not modes:
            continue
        best = modes[0]
        # closest approach of any atom of the top mode to each reference
        d_active = min(math.dist(a, ACTIVE) for a in best[2])
        d_mo = min(math.dist(a, MO) for a in best[2])
        # centroids of each of the top modes, to see whether they agree with each other
        cents = [tuple(sum(a[i] for a in m[2]) / len(m[2]) for i in range(3)) for m in modes]
        spread = max(math.dist(cents[0], c) for c in cents) if len(cents) > 1 else 0.0
        verdict = "at the active site" if d_active < 6.0 else "elsewhere on the protein"
        results.append({"ligand": name, "best_affinity": best[1],
                        "closest_atom_to_active_site_A": round(d_active, 2),
                        "closest_atom_to_Mo_A": round(d_mo, 2),
                        "top3_centroid_spread_A": round(spread, 2),
                        "verdict": verdict,
                        "top3_centroids": [[round(c, 2) for c in cc] for cc in cents]})
        print(f"{name:12s} best {best[1]:>6} kcal/mol | closest atom to active site "
              f"{d_active:6.2f} A, to Mo {d_mo:6.2f} A | top-3 centroid spread {spread:5.2f} A"
              f"  -> {verdict}")

    if results:
        n_at = sum(1 for r in results if r["verdict"] == "at the active site")
        print(f"\n{n_at}/{len(results)} probe ligands place their best pose at the active site.")
        if n_at < len(results):
            print("The remaining probes prefer a different pocket, so the focused-box scores in the")
            print("main redo describe the active site under a constrained search, not the global")
            print("optimum of the receptor. Both should be reported.")
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump({"active_site_reference": [round(c, 3) for c in ACTIVE],
                   "mo": MO, "so4_s": SO4_S, "results": results}, f, indent=2)
    print(f"\n-> {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
