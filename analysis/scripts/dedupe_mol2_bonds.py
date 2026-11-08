#!/usr/bin/env python3
"""dedupe_mol2_bonds.py — remove duplicate bonds from an antechamber mol2.

Why this is needed
------------------
tleap refused the SUOX molybdenum cofactor:

    Atom .R<MTE 1>.A<C4 5> has force field coordination 4
    but only 3 bonded neighbors.

and the cause is in the mol2 rather than in the geometry: antechamber wrote the 5-6 bond **twice**.

    10     5     6 1
    11     5     6 1

tleap accumulates bonds, so the duplicate consumes two of C4's valence slots and the atom ends up
over-coordinated without having the neighbours to match - hence the message. Removing the repeat
fixes the valence.

Note on a previous failed attempt: this was first tried as an inline heredoc inside the SSH
command, and the heredoc silently did not execute (no output, no change to the file). Scripts that
must run are therefore written as files and uploaded, never pasted as heredocs.

Usage:
    python3 dedupe_mol2_bonds.py <file.mol2> [more.mol2 ...]
"""
import re
import sys

BOND_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$")


def main(paths):
    bad = 0
    for path in paths:
        lines = open(path, encoding="utf-8", errors="replace").read().split("\n")
        out, seen, in_bonds, removed = [], set(), False, []

        for ln in lines:
            if ln.startswith("@<TRIPOS>BOND"):
                in_bonds = True
                out.append(ln)
                continue
            if ln.startswith("@<TRIPOS>") and in_bonds:
                in_bonds = False

            if in_bonds:
                m = BOND_RE.match(ln)
                if m:
                    a, b = int(m.group(2)), int(m.group(3))
                    key = (min(a, b), max(a, b))
                    if key in seen:
                        removed.append((a, b, m.group(4)))
                        continue
                    seen.add(key)

            out.append(ln)

        if removed:
            open(path, "w", encoding="utf-8", newline="\n").write("\n".join(out))
            print(f"{path}: removed {len(removed)} duplicate bond(s): "
                  f"{', '.join(f'{a}-{b}(order {o})' for a, b, o in removed)}")
        else:
            print(f"{path}: no duplicate bonds found ({len(seen)} bonds)")

        # re-read and report, so the fix is verified rather than assumed
        lines2 = open(path, encoding="utf-8", errors="replace").read().split("\n")
        n_bonds = sum(1 for ln in lines2 if BOND_RE.match(ln))
        n_atoms = 0
        in_atoms = False
        for ln in lines2:
            if ln.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if ln.startswith("@<TRIPOS>") and in_atoms:
                in_atoms = False
            elif in_atoms and ln.strip():
                n_atoms += 1
        # valence sanity check for every heavy atom, using the bond list
        deg = {}
        for ln in lines2:
            m = BOND_RE.match(ln)
            if m:
                for i in (int(m.group(2)), int(m.group(3))):
                    deg[i] = deg.get(i, 0) + 1
        over = [i for i, d in deg.items() if d > 4]
        print(f"    verified: {n_atoms} atoms, {n_bonds} bonds"
              f"{', OVER-COORDINATED atoms: ' + str(over) if over else ', no atom has >4 bonds'}")
        if over:
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1:]))
