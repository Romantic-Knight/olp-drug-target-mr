#!/usr/bin/env bash
# 01_param_all.sh — force-field parameters for the ligands and the Mo-MPT site.
#
# Run inside the MD environment, from analysis/md/, with the complexes already built:
#     bash scripts/01_param_all.sh
#
# Produces, per complex, a ligand .mol2/.frcmod (GAFF2 + AM1-BCC) and — for SUOX — an MCPB.py bonded
# model for the Mo centre. Both are needed before tleap can build a topology.
#
# ⚠️ THE MO-MPT SITE IS THE HIGHEST-RISK STEP IN THIS PROJECT.
# MCPB.py is designed for a metal plus its coordinating residues. Here the "residue" that donates
# two of the five coordinating atoms is MTE, a 24-atom organometallic cofactor, and the other three
# donors are Cys185 SG and two water oxygens. Read section 0's output, check the generated model,
# and treat a first-attempt failure as expected rather than surprising. If MCPB.py cannot handle
# MTE, the fallback is a QM/MM or a published Mo-MPT parameter set — see README section "Known
# limitations", and note that any such change must be reported in the paper.
set -euo pipefail
cd "$(dirname "$0")/.."

COMPLEXES="${COMPLEXES:-GDF11_GW8510 SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
NET_CHARGE_FILE="complexes/net_charges.txt"

# ---------------------------------------------------------------- ligand net charges
# antechamber needs the total charge. It is taken from the docking ligand manifest, which recorded
# the protonation policy applied to each molecule, rather than guessed here.
if [ ! -f "$NET_CHARGE_FILE" ]; then
    echo "Building $NET_CHARGE_FILE from the docking ligand manifest..."
    python - <<'PY'
import json, os, re
root = os.path.dirname(os.path.abspath("__file__"))
man = json.load(open("../docking/ligands/ligand_manifest.json", encoding="utf-8"))
seen = {}
for key, v in man["ligands"].items():
    seen.setdefault(v["name"], (v["smiles"], v.get("net_charge", 0)))
os.makedirs("complexes", exist_ok=True)
with open("complexes/net_charges.txt", "w", encoding="utf-8") as f:
    f.write("# ligand\tcharge\tsmiles\n")
    for n, (smi, q) in sorted(seen.items()):
        f.write(f"{n}\t{q}\t{smi}\n")
print("wrote complexes/net_charges.txt")
PY
fi

charge_of() {  # charge_of <ligand_pdb_stem>
    local want="$1"
    awk -F'\t' -v w="$want" '$1==w {print $2; found=1} END{if(!found) print "NA"}' "$NET_CHARGE_FILE"
}

# ---------------------------------------------------------------- ligands (GAFF2 / AM1-BCC)
for cx in $COMPLEXES; do
    dir="complexes/$cx"
    [ -d "$dir" ] || { echo "skip $cx (no such directory)"; continue; }
    lig=$(awk '/^HETATM/ {print $4; exit}' "$dir/ligand.pdb")
    name=$(python -c "
import json;m=json.load(open('complexes/complex_manifest.json',encoding='utf-8'))
print(m['complexes']['$cx']['ligand'])")
    q=$(charge_of "$name")
    if [ "$q" = "NA" ]; then
        echo "!! $cx: no net charge recorded for ligand '$name'; add it to $NET_CHARGE_FILE"
        continue
    fi
    echo "=== $cx  ligand=$name  net charge=$q"
    (
        cd "$dir"
        # -c bcc: AM1-BCC charges via sqm (semi-empirical, free); -at gaff2: the GAFF2 atom types
        antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 \
                    -c bcc -s 2 -at gaff2 -nc "$q" -rn LIG \
            > antechamber.log 2>&1
        parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -s gaff2 > parmchk2.log 2>&1
        if [ -s ligand.mol2 ] && [ -s ligand.frcmod ]; then
            echo "    ligand.mol2 + ligand.frcmod written"
        else
            echo "    !! ligand parameterisation failed - see antechamber.log"
        fi
    )
done

# ---------------------------------------------------------------- Mo-MPT bonded model
for cx in $COMPLEXES; do
    case "$cx" in SUOX_*) ;; *) continue ;; esac
    dir="complexes/$cx"
    echo "=== $cx  Mo-MPT site (MCPB.py bonded model)"
    cat > "$dir/mcpb_$cx.in" <<EOF
# MCPB.py input for the SUOX molybdenum centre.
#
# Coordination shell read directly from 2A9A (see analysis/scripts/_probe_mo_coordination.py):
#   MTE 1501 S1'  2.38 A      Mo-molybdopterin dithiolene
#   MTE 1501 S2'  2.43 A      Mo-molybdopterin dithiolene
#   CYS 185  SG   2.41 A      cysteine ligand
#   HOH 7505 O    1.69 A      oxo ligand
#   HOH 7506 O    2.10 A      hydroxo / water
#
# VERIFY THIS FILE AGAINST THE MCPB.py MANUAL BEFORE RUNNING. The keyword layout below follows the
# MCPB.py documentation, but the tool is strict and version-dependent; a rejected input costs
# minutes to fix, an accepted-but-wrong one costs the whole run.
large_mol = receptor.pdb
small_mol = ligand.pdb
ion_ids = 2501
ion_mol2files = MO.mol2
CYS 185 SG
MTE 1501 S1'
MTE 1501 S2'
HOH 7505 O
HOH 7506 O
EOF
    echo "    wrote $dir/mcpb_$cx.in"
    echo "    steps to run in that directory, once the input is checked:"
    echo "      MCPB.py -i mcpb_$cx.in -s 1   # build the small model"
    echo "      MCPB.py -i mcpb_$cx.in -s 2   # write the QM input"
    echo "      sqm  -i ... -o ...            # or Gaussian/ORCA, for the QM charges"
    echo "      MCPB.py -i mcpb_$cx.in -s 3   # fit charges (RESP)"
    echo "      MCPB.py -i mcpb_$cx.in -s 4   # emit the tleap input"
done

echo
echo "parameterisation stage finished; inspect the logs before building the topology"
