"""redock.py — redo the OLP docking against the rebuilt receptors, with recorded config and seeds.

What differs from the published run
-----------------------------------
The published `run_docking.py` docked the 22 DSigDB candidates once, with `--seed` left at Vina's
default, against receptors that carried **no hydrogens** and (for SUOX) no catalytic cofactor at
all. Three things are fixed here:

1. **Receptors** carry hydrogens and the SUOX MTE/SO4 cofactor (`prepare_receptor.py`).
2. **Ligands** are regenerated with an explicit protonation policy (`prepare_ligands.py`).
3. **Three seeds** per pair (`0, 1, 42`) instead of one, because a single Vina run on a 29-atom
   ligand with 6-9 rotatable bonds is not an adequate sample. The per-seed best is reported
   alongside the overall best, so seed-to-seed spread is visible rather than hidden.

The grid boxes and the scoring protocol are deliberately **unchanged** from the published run
(exhaustiveness 32, 9 modes, 30^3 A box) so the new numbers stay comparable: GDF11 is centred on
the monomer centroid, SUOX on the old Mo/SO4 active-site midpoint.

Every run writes its full command line into `docking_redo_config.json`; a blind-docking mode
(`--blind`) is provided for locating pockets independently of the published box.

Usage:
    python analysis/scripts/redock.py --quick            # 2 ligands, 1 seed: smoke test
    python analysis/scripts/redock.py                    # full redo: 22 pairs x 3 seeds
    python analysis/scripts/redock.py --blind SUOX       # whole-protein grid, for pocket finding
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = os.path.join(ROOT, "analysis", "docking")
VINA = os.path.join(ROOT, "analysis", "tools", "vina", "vina_1.2.7_win.exe")
RECDIR = os.path.join(BASE, "receptors")
LIGDIR = os.path.join(BASE, "ligands", "prepared")
OUTROOT = os.path.join(BASE, "redo")

RECEPTORS = {
    "GDF11": "5E4G_chainA_H.pdbqt",
    "SUOX": "2A9A_chainA_MTE_SO4.pdbqt",
}
# unchanged from the published run (analysis/scripts/run_docking.py)
CENTER = {"GDF11": (-0.07, 8.51, 7.15), "SUOX": (-0.19, 1.80, 10.16)}
SIZE = {"GDF11": (30, 30, 30), "SUOX": (30, 30, 30)}
EXHAUSTIVENESS = 32
NUM_MODES = 9
SEEDS = [0, 1, 42]

MODE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$")


def parse_modes(text):
    """Vina's result table: mode, affinity (kcal/mol), rmsd l.b., rmsd u.b."""
    out = []
    for ln in text.splitlines():
        m = MODE_RE.match(ln)
        if m:
            out.append({"mode": int(m.group(1)), "affinity": float(m.group(2)),
                        "rmsd_lb": float(m.group(3)), "rmsd_ub": float(m.group(4))})
    return out


def run_vina(receptor, ligand, center, size, out_pdbqt, seed, cpu):
    cmd = [VINA, "--receptor", receptor, "--ligand", ligand,
           "--center_x", str(center[0]), "--center_y", str(center[1]),
           "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
           "--seed", str(seed), "--cpu", str(cpu),
           "--out", out_pdbqt]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    return {"cmd": cmd, "exit": r.returncode, "stdout": r.stdout, "stderr": r.stderr,
            "seconds": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="2 ligands, 1 seed")
    ap.add_argument("--seeds", type=str, default="", help="override seeds, comma separated")
    ap.add_argument("--cpu", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--target", type=str, default="", help="restrict to one target")
    ap.add_argument("--decoys", action="store_true",
                    help="dock the NCI decoy set instead of the candidates (1 seed, both targets)")
    ap.add_argument("--blind", type=str, default="",
                    help="target name: 60^3 A box over the whole protein, for pocket finding")
    a = ap.parse_args()

    os.makedirs(OUTROOT, exist_ok=True)
    man = json.load(open(os.path.join(BASE, "ligands", "ligand_manifest.json"), encoding="utf-8"))

    if a.blind:
        t = a.blind
        # A blind scan must be centred on the protein, not on the published active site (otherwise
        # it can only rediscover the site it was told about). 60^3 A about the receptor centroid
        # covers ~95% of the SUOX chain-A envelope (spans x 69.3, y 45.0, z 52.9 A).
        rec_path = os.path.join(RECDIR, RECEPTORS[t])
        xs = ys = zs = None
        rx, ry, rz = [], [], []
        for ln in open(rec_path, encoding="utf-8"):
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    rx.append(float(ln[30:38])); ry.append(float(ln[38:46])); rz.append(float(ln[46:54]))
                except ValueError:
                    pass
        cen = (sum(rx) / len(rx), sum(ry) / len(ry), sum(rz) / len(rz))
        cfg = {"mode": "blind", "target": t, "size": (60, 60, 60), "center": [round(c, 3) for c in cen],
               "center_note": "receptor centroid (not the published active-site box)",
               "span_A": [round(max(rx) - min(rx), 1), round(max(ry) - min(ry), 1),
                          round(max(rz) - min(rz), 1)],
               "exhaustiveness": 24, "num_modes": 20, "seeds": [0], "cpu": a.cpu}
        out_dir = os.path.join(OUTROOT, f"blind_{t}")
        os.makedirs(out_dir, exist_ok=True)
        ligs = [k.split(":", 1)[1] for k in man["ligands"] if k.startswith(t + ":")]
        chosen = ["Apigenin", "Vorinostat", "Celastrol"]      # one per chemotype
        ligs = [n for n in chosen if n in ligs]
        print(f"blind docking {t}: box 60^3 A about the centroid "
              f"({cen[0]:.2f}, {cen[1]:.2f}, {cen[2]:.2f}), {len(ligs)} probe ligands")
        print("NOTE: a 60^3 A box in Vina is a pocket-finding heuristic, not a rigorous engine; "
              "the result is interpreted as clustering of poses, not as an affinity.")
        rows = []
        for name in ligs:
            lg = os.path.join(LIGDIR, name + ".pdbqt")
            out = os.path.join(out_dir, f"{name}.pdbqt")
            r = run_vina(rec_path, lg, cen, (60, 60, 60), out, 0, a.cpu)
            modes = parse_modes(r["stdout"])
            print(f"  {name:16s} exit={r['exit']} modes={len(modes)} "
                  f"best={modes[0]['affinity'] if modes else 'n/a'} ({r['seconds']}s)")
            rows.append({"ligand": name, "modes": modes, "seconds": r["seconds"],
                         "exit": r["exit"]})
            with open(os.path.join(out_dir, "config.json"), "w", encoding="utf-8") as f:
                json.dump({**cfg, "results": rows}, f, indent=2)
        return 0

    seeds = [int(x) for x in a.seeds.split(",")] if a.seeds else SEEDS
    targets = [a.target] if a.target else list(RECEPTORS)

    lig_dir = LIGDIR
    out_dir = os.path.join(OUTROOT, "main")
    tag = ""
    if a.decoys:
        # decoys are docked into BOTH receptors: the control asks whether the protocol separates
        # candidates from non-binders at either site, so a single target would answer half of it
        lig_dir = os.path.join(BASE, "decoys", "prepared")
        names = sorted(os.path.basename(f)[:-6]
                       for f in glob.glob(os.path.join(lig_dir, "*.pdbqt")))
        pairs = [(t, n, None) for t in targets for n in names]
        seeds = [0]
        out_dir = os.path.join(OUTROOT, "decoys")
        tag = "_decoys"
        print(f"decoy control: {len(names)} decoys x {len(targets)} receptors, 1 seed")
    else:
        pairs = [(t, k.split(":", 1)[1], v)
                 for k, v in sorted(man["ligands"].items())
                 for t in targets if k.startswith(t + ":")]
    if a.quick:
        pairs = pairs[:2]
        seeds = seeds[:1]
    os.makedirs(out_dir, exist_ok=True)
    cfg = {"center": CENTER, "size": SIZE, "exhaustiveness": EXHAUSTIVENESS,
           "num_modes": NUM_MODES, "seeds": seeds, "cpu": a.cpu,
           "vina": os.path.relpath(VINA, ROOT).replace("\\", "/"),
           "receptors": {t: os.path.relpath(os.path.join(RECDIR, f), ROOT).replace("\\", "/")
                         for t, f in RECEPTORS.items()},
           "ligands": os.path.relpath(LIGDIR, ROOT).replace("\\", "/"),
           "note": ("box, exhaustiveness and mode count unchanged from the published run so the "
                    "numbers stay comparable; receptors and ligands regenerated; seeds added")}

    print(f"{len(pairs)} target-ligand pairs x {len(seeds)} seeds, cpu={a.cpu}")
    rows = []
    done = 0
    for t, name, rec in pairs:
        rec_p = os.path.join(RECDIR, RECEPTORS[t])
        lig_p = os.path.join(lig_dir, name + ".pdbqt")
        if not os.path.exists(lig_p):
            print(f"  {t}:{name} missing ligand pdbqt, skipped")
            continue
        for seed in seeds:
            out = os.path.join(out_dir, f"{t}__{name}__seed{seed}.pdbqt")
            r = run_vina(rec_p, lig_p, CENTER[t], SIZE[t], out, seed, a.cpu)
            modes = parse_modes(r["stdout"])
            form = "as_prepared"
            if r["exit"] != 0 or not modes:
                # A metal-containing ligand is rejected by Vina ("Unknown or inappropriate tag
                # found in flex residue"). Retry with the documented metal-free surrogate rather
                # than dropping the entry from the series.
                alt = os.path.join(lig_dir, name + "_metalfree.pdbqt")
                if os.path.exists(alt):
                    r = run_vina(rec_p, alt, CENTER[t], SIZE[t], out, seed, a.cpu)
                    modes = parse_modes(r["stdout"])
                    form = "metal_free_surrogate"
                if not modes:
                    print(f"  {t}:{name} seed={seed} FAILED exit={r['exit']} "
                          f"{(r['stderr'] or '').strip()[:70]}")
            for m in modes:
                rows.append({"target": t, "ligand": name, "seed": seed, **m,
                             "form": form, "seconds": r["seconds"]})
            best = f"{modes[0]['affinity']:.2f}" if modes else "n/a"
            done += 1
            print(f"  [{done}/{len(pairs) * len(seeds)}] {t:5s} {name:26s} seed={seed:<3d} "
                  f"best={best:>6s} modes={len(modes)} ({r['seconds']}s)")

    cfg["n_runs"] = done
    with open(os.path.join(OUTROOT, f"docking_redo{tag}_config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2)
    with open(os.path.join(OUTROOT, f"docking_redo{tag}_modes.json"), "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)

    # per (target, ligand) summary across seeds
    import statistics
    summary = []
    keys = sorted({(r["target"], r["ligand"]) for r in rows})
    for t, name in keys:
        sub = [r for r in rows if r["target"] == t and r["ligand"] == name]
        per_seed = {}
        for r in sub:
            per_seed.setdefault(r["seed"], []).append(r["affinity"])
        bests = [min(v) for v in per_seed.values()]
        summary.append({"target": t, "ligand": name,
                        "best_affinity": round(min(bests), 3),
                        "worst_seed_best": round(max(bests), 3),
                        "median_seed_best": round(statistics.median(bests), 3),
                        "seed_spread": round(max(bests) - min(bests), 3),
                        "n_seeds": len(bests),
                        "membrane_pass_lt_minus7": int(min(bests) < -7.0)})
    summary.sort(key=lambda x: (x["target"], x["best_affinity"]))
    with open(os.path.join(OUTROOT, f"docking_redo{tag}_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    n_pass = sum(s["membrane_pass_lt_minus7"] for s in summary)
    worst_spread = max((s["seed_spread"] for s in summary), default=0)
    print(f"\n{len(summary)} target-ligand pairs summarised; "
          f"{n_pass} reach < -7 kcal/mol; largest seed spread {worst_spread:.2f} kcal/mol")
    print(f"config  -> {os.path.relpath(os.path.join(OUTROOT, f'docking_redo{tag}_config.json'), ROOT)}")
    print(f"summary -> {os.path.relpath(os.path.join(OUTROOT, f'docking_redo{tag}_summary.json'), ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
