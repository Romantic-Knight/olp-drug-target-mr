# -*- coding: utf-8 -*-
"""Extract the 6 per-sample count matrices from GSE211630_RAW.tar.

The tar contains *_count_data-*.txt.gz (still gzip-compressed inside). We extract
them as-is into data/raw/geo/GSE211630/matrix/. tarfile mode 'r' with gz members.
"""
import os
import tarfile

TAR = os.path.join("data", "raw", "geo", "GSE211630", "GSE211630_RAW.tar")
OUT = os.path.join("data", "raw", "geo", "GSE211630", "matrix")
os.makedirs(OUT, exist_ok=True)

with tarfile.open(TAR, "r") as tf:
    for m in tf.getmembers():
        if not m.isfile():
            continue
        # member is stored as-is (gzip file inside tar)
        src = tf.extractfile(m)
        dest = os.path.join(OUT, os.path.basename(m.name))
        if os.path.exists(dest) and os.path.getsize(dest) == m.size:
            print(f"[skip] {os.path.basename(m.name)} already present")
            continue
        with open(dest, "wb") as f:
            f.write(src.read())
        print(f"[ok] {os.path.basename(m.name)} -> {m.size} bytes")
print("done")
