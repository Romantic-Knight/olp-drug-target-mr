# Run AutoDock Vina docking for GDF11 (5E4G) and SUOX (2A9A) DSigDB candidates.
# Splits multi-model ligand PDBQT, runs vina per ligand, collects top affinity.
import os, subprocess, re, sys, time

BASE = r"D:\AI-Wisp-Science\Drug Target MR\analysis\docking"
VINA = os.path.join(BASE, r"..\tools\vina\vina_1.2.7_win.exe")
VINA = os.path.normpath(VINA)
RECEPTORS = {
    "GDF11": os.path.join(BASE, "receptors", "5E4G_receptor.pdbqt"),
    "SUOX": os.path.join(BASE, "receptors", "2A9A_receptor.pdbqt"),
}
CENTER = {
    "GDF11": (-0.07, 8.51, 7.15),   # protein centroid, monomer
    "SUOX": (-0.19, 1.80, 10.16),   # MO/SO4 active-site midpoint
}
SIZE = {"GDF11": (30, 30, 30), "SUOX": (30, 30, 30)}

def split_ligands(pdbqt_path, outdir):
    """split multi-MODEL pdbqt using official vina_split; names come from REMARK Name = X
    in each MODEL block (same order). Returns {name: path}"""
    split_exe = os.path.normpath(os.path.join(BASE, r"..\tools\vina\vina_split_1.2.7_win.exe"))
    # fresh dir
    if os.path.exists(outdir):
        import shutil; shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)
    subprocess.run([split_exe, "--input", pdbqt_path, "--ligand", outdir + os.sep],
                   capture_output=True, text=True)
    files = sorted(f for f in os.listdir(outdir) if f.endswith(".pdbqt"))
    # names in MODEL order
    txt = open(pdbqt_path).read()
    names = [m.strip() for m in re.findall(r"REMARK\s+Name\s*=\s*(.+)", txt)]
    out = {}
    for i, fn in enumerate(files):
        name = names[i] if i < len(names) else "lig_%d" % i
        safe = re.sub(r"[^\w.-]", "_", name)
        src = os.path.join(outdir, fn)
        dst = os.path.join(outdir, safe + ".pdbqt")
        os.replace(src, dst)
        out[safe] = dst
    return out

def run_vina(receptor, ligand, cx, cy, cz, sx, sy, sz, out_pdbqt, cpu=4, exhaust=32):
    cmd = [VINA, "--receptor", receptor, "--ligand", ligand,
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
           "--exhaustiveness", str(exhaust), "--num_modes", "9",
           "--out", out_pdbqt, "--cpu", str(cpu)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr

def parse_affinity(stdout):
    # vina writes "REMARK VINA RESULT:" lines in the output PDBQT, but stdout has a
    # "mode | affinity | rmsd l.b.| rmsd u.b." table; parse the table rows after the header
    best = None
    for line in stdout.splitlines():
        m = re.match(r"\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", line)
        if m and not line.lstrip().startswith(("mode", "-----")):
            a = float(m.group(1))
            if best is None or a < best:
                best = a
    return best

def main():
    targets = sys.argv[1].split(",") if len(sys.argv) > 1 else ["GDF11", "SUOX"]
    ligdir = os.path.join(BASE, "ligands")
    outdir = os.path.join(BASE, "docked")
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for t in targets:
        lf = os.path.join(ligdir, t + ".pdbqt")
        if not os.path.exists(lf):
            print("no ligand file", lf); continue
        split_dir = os.path.join(BASE, "ligands", t + "_split")
        os.makedirs(split_dir, exist_ok=True)
        ligs = split_ligands(lf, split_dir)
        rec = RECEPTORS[t]
        cx, cy, cz = CENTER[t]
        sx, sy, sz = SIZE[t]
        print("[%s] %d ligands; receptor %s" % (t, len(ligs), os.path.basename(rec)), flush=True)
        for name, lp in ligs.items():
            op = os.path.join(outdir, "%s_%s_out.pdbqt" % (t, name))
            t0 = time.time()
            rc, so, se = run_vina(rec, lp, cx, cy, cz, sx, sy, sz, op)
            aff = parse_affinity(so) if rc == 0 else None
            rows.append((t, name, aff, rc, so.splitlines()[-1] if so else se.strip().splitlines()[-1] if se else ""))
            print("  %-30s rc=%d best=%.2f  (%.0fs)" % (name, rc, aff if aff else -99, time.time()-t0), flush=True)
    # write results
    with open(os.path.join(BASE, "docking_results.tsv"), "w") as f:
        f.write("target\tligand\tbest_affinity_kcal_mol\trc\tnote\n")
        for t, n, a, rc, note in rows:
            f.write("%s\t%s\t%s\t%d\t%s\n" % (t, n, ("%.2f" % a) if a else "NA", rc, note.replace("\t", " ")))
    print("done; %d rows -> docking_results.tsv" % len(rows))

if __name__ == "__main__":
    main()
