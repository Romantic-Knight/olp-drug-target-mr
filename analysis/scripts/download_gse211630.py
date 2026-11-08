# -*- coding: utf-8 -*-
"""Download GSE211630 supplementary tar + series matrix into data/raw/geo/GSE211630/.

Downloads (https):
  - GSE211630_RAW.tar   (~110 MB, 6 per-sample count matrices)
  - GSE211630_series_matrix.txt.gz  (series metadata)
Also verifies tar integrity with tarfile + reports member listing.
Run inside the project root (run context cwd is project root).
"""
import os
import sys
import tarfile
import urllib.request

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211630"
OUT = os.path.join("data", "raw", "geo", "GSE211630")
os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(OUT, "matrix"), exist_ok=True)

UA = {"User-Agent": "wisp-science-OLP-MR"}


def download(url, dest, label):
    tmp = dest + ".part"
    if os.path.exists(dest):
        print(f"[skip] {label}: already at {dest}")
        return True
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as r:
        total = int(r.headers.get("Content-Length") or 0)
        done = 0
        with open(tmp, "wb") as f:
            while True:
                chunk = r.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
                done += len(chunk)
                pct = (100.0 * done / total) if total else 0
                if done % (8 * 1024 * 1024) < 1024 * 1024:
                    print(f"[{label}] {done/1e6:.1f}/{total/1e6:.1f} MB ({pct:.0f}%)",
                          flush=True)
    os.replace(tmp, dest)
    print(f"[done] {label}: {os.path.getsize(dest)} bytes")
    return True


# 1) series matrix metadata (tiny)
sm = os.path.join(OUT, "matrix", "GSE211630_series_matrix.txt.gz")
download(f"{BASE}/matrix/GSE211630_series_matrix.txt.gz", sm, "series_matrix")

# 2) supplementary raw tar (big)
tar = os.path.join(OUT, "GSE211630_RAW.tar")
download(f"{BASE}/suppl/GSE211630_RAW.tar", tar, "RAW.tar")

# 3) integrity + member listing
sz = os.path.getsize(tar)
print("tar bytes:", sz)
with tarfile.open(tar, "r") as tf:
    names = tf.getnames()
    bad = tf.members  # detect truncated: open(mode='r') raises on truncation
print("members:", len(names))
for n in names:
    info = [m for m in tarfile.open(tar).getmembers() if m.name == n][0]
    print(f"  {n}  {info.size} bytes")
print("ALL_OK")
