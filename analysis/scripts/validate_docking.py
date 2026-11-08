"""validate_docking.py — positive control for the SUOX docking setup.

A docking protocol is only interpretable if it can recover a known binding geometry. 2A9A contains
the catalytic sulfate (SO4) sitting in the molybdenum site; here it is extracted, docked back into
a receptor from which it was removed, and the recovered position is compared with the crystal one.

Interpretation: S-atom displacement < 1.5 A means the box, receptor protonation and scoring
protocol are sound, so the ligand scores in the main redo are meaningful. A larger displacement
means the box or the receptor preparation is wrong and the main results should not be trusted.

Writes analysis/docking/redo/validation_so4.json.

Usage:
    python analysis/scripts/validate_docking.py [--cpu 8]
"""
import argparse
import json
import math
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE = os.path.join(ROOT, "analysis", "docking")
VINA = os.path.join(ROOT, "analysis", "tools", "vina", "vina_1.2.7_win.exe")
OBABEL = os.path.join(BASE, ".pixi", "envs", "default", "Library", "bin", "obabel.exe")
SRC = os.path.join(BASE, "structures", "2A9A.pdb")
WORK = os.path.join(BASE, "redo", "validation")
CHAIN = "A"


def xyz(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)

    # ---- 1. the crystal SO4, and its centroid -------------------------------------------------
    # 2A9A holds TWO sulfates in chain A. Only res 7503 sits at the catalytic centre (5.1 A from the
    # Mo); res 7504 is a surface ion ~10 A away. The published SUOX grid centre is exactly the
    # midpoint of Mo and 7503, so 7503 is the catalytic one and the only valid control ligand.
    all_h = [ln for ln in open(SRC, encoding="utf-8")
             if ln.startswith("HETATM") and ln[21] == CHAIN]
    mo = [xyz(ln) for ln in all_h if ln[17:20].strip() == "MO"]
    if not mo:
        print("no Mo in the crystal file - cannot identify the catalytic sulfate")
        return 1
    mo_xyz = mo[0]
    by_res = {}
    for ln in all_h:
        if ln[17:20].strip() == "SO4":
            by_res.setdefault(ln[22:26].strip(), []).append(ln)
    if not by_res:
        print("no SO4 found - cannot run the control")
        return 1

    def s_of(lines):
        for ln in lines:
            if ln[12:16].strip() == "S":
                return xyz(ln)
        return xyz(lines[0])

    ranked = sorted(by_res.items(), key=lambda kv: math.dist(s_of(kv[1]), mo_xyz))
    res_id, so4_lines = ranked[0]
    s_xyz_crystal = s_of(so4_lines)
    print(f"crystal Mo at ({mo_xyz[0]:.3f}, {mo_xyz[1]:.3f}, {mo_xyz[2]:.3f})")
    for rid, ls in ranked:
        print(f"  sulfate res {rid}: S at "
              f"({s_of(ls)[0]:.3f}, {s_of(ls)[1]:.3f}, {s_of(ls)[2]:.3f}), "
              f"{math.dist(s_of(ls), mo_xyz):.2f} A from Mo")
    print(f"-> using the catalytic sulfate res {res_id} as the control ligand")
    # the box is centred on the crystallographic sulfate itself: a self-docking recovery test
    cen = s_xyz_crystal

    so4_pdb = os.path.join(WORK, "so4.pdb")
    with open(so4_pdb, "w", encoding="utf-8") as f:
        f.writelines(so4_lines)
        f.write("END\n")

    # ---- 2. a receptor with MTE but WITHOUT SO4 ----------------------------------------------
    rec_noSo4_pdb = os.path.join(WORK, "2A9A_chainA_MTE_noSO4.pdb")
    kept = 0
    with open(rec_noSo4_pdb, "w", encoding="utf-8") as f:
        for ln in open(SRC, encoding="utf-8"):
            if ln.startswith("HETATM") and ln[21] == CHAIN:
                res = ln[17:20].strip()
                if res == "SO4":
                    continue                                    # the ligand must go
                if res == "MO":
                    continue                                    # Vina cannot parse Mo
                f.write(ln)                                     # HOH is dropped by obabel -h
                kept += 1
            elif ln.startswith("ATOM") and ln[21] == CHAIN:
                f.write(ln)
                kept += 1
        f.write("END\n")
    print(f"receptor rebuilt without SO4: {kept} HETATM/ATOM lines kept")

    rec_pdbqt = os.path.join(WORK, "2A9A_chainA_MTE_noSO4.pdbqt")
    r = subprocess.run([OBABEL, "-ipdb", rec_noSo4_pdb, "-opdbqt", "-xr", "-h", "-O", rec_pdbqt],
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(rec_pdbqt):
        print("receptor -> pdbqt failed:", (r.stderr or "").strip()[:200])
        return 1
    # the rigid-receptor flavour must carry no flex tags (the same fix as prepare_receptor.py)
    lines = [ln for ln in open(rec_pdbqt, encoding="utf-8")
             if not ln.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH",
                                   "TORSDOF", "BEGIN_RES", "END_RES"))]
    with open(rec_pdbqt, "w", encoding="utf-8", newline="\n") as f:
        f.writelines(lines)

    # ---- 3. dock the sulfate back --------------------------------------------------------------
    lig_pdbqt = os.path.join(WORK, "so4.pdbqt")
    r = subprocess.run([OBABEL, "-ipdb", so4_pdb, "-opdbqt", "-O", lig_pdbqt],
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(lig_pdbqt):
        print("SO4 -> pdbqt failed:", (r.stderr or "").strip()[:200])
        return 1

    out = os.path.join(WORK, "so4_redocked.pdbqt")
    box = 16
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cen[0]:.3f}", "--center_y", f"{cen[1]:.3f}", "--center_z", f"{cen[2]:.3f}",
           "--size_x", str(box), "--size_y", str(box), "--size_z", str(box),
           "--exhaustiveness", "32", "--num_modes", "9", "--seed", "0",
           "--cpu", str(a.cpu), "--out", out]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(out):
        print("vina failed:", (res.stderr or "").strip()[:200])
        return 1

    # ---- 4. displacement of the best mode's S atom from the crystal S ------------------------
    crystal_S = s_xyz_crystal
    best_aff, s_xyz, n_models = None, None, 0
    for ln in open(out, encoding="utf-8"):
        if ln.startswith("MODEL"):
            n_models += 1
            if n_models > 1:
                break
            continue
        if n_models == 1 and ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "S":
            s_xyz = xyz(ln)
    for ln in res.stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.\d+)", ln)
        if m:
            best_aff = float(m.group(1))
            break

    result = {"box_A": box, "crystal_centroid": [round(c, 3) for c in cen],
              "crystal_S": [round(c, 3) for c in crystal_S] if crystal_S else None,
              "docked_S": [round(c, 3) for c in s_xyz] if s_xyz else None,
              "best_affinity": best_aff, "n_modes": n_models}
    if crystal_S and s_xyz:
        d = math.dist(crystal_S, s_xyz)
        result["S_displacement_A"] = round(d, 3)
        result["pass"] = bool(d < 1.5)
        print(f"\ndocked sulfate best affinity {best_aff} kcal/mol; "
              f"S displaced {d:.2f} A from the crystal position")
        print("POSITIVE CONTROL:", "PASS" if d < 1.5 else "FAIL",
              "- box and receptor preparation recover the crystallographic site"
              if d < 1.5 else "- the protocol does not recover the known site")
    else:
        result["pass"] = False
        print("\ncould not compare S positions; control inconclusive")

    with open(os.path.join(BASE, "redo", "validation_so4.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"-> {os.path.relpath(os.path.join(BASE, 'redo', 'validation_so4.json'), ROOT)}")
    return 0 if result.get("pass") else 1


if __name__ == "__main__":
    sys.exit(main())
