# -*- coding: utf-8 -*-
"""Stage D (LIANA/CellChat cell-cell communication, GSE211630) table verification.

Independently re-derives the canonical numbers behind Fig D1 and checks the
final tables for internal consistency and cross-table agreement:

  * geo_gse211630_liana_signal_shift.csv        (merged differential table)
  * geo_gse211630_liana_headlines.csv           (40 presence headlines)
  * geo_gse211630_liana_incoming_focus.csv      (per receiver x source x class)
  * geo_gse211630_liana_receptor_module_corr.csv(392 receptor x module rho)
  * geo_gse211630_liana_figd1_summary.csv|_ligands.csv (canonical Fig D1 numbers)

Conventions locked for Stage D (identical to Fig D1):
  class: EOLP_lost = significant only in NEOLP; sustained = significant in both
         conditions (also both-present by construction); EOLP_gained =
         significant only in EOLP; none = significant in neither.
  focus receivers = {0,1,3,4,5,6,7,8,16}; T/NK focus = {3,4,5,6,7,8,16}.
Output: results/reports/stageD_table_verification.txt
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "tables"
OUT = ROOT / "results" / "reports" / "stageD_table_verification.txt"

RECEIVERS = [16, 6, 1, 3, 8, 4, 7, 5, 0]
TN_FOCUS = {3, 4, 5, 6, 7, 8, 16}
FOCUS = set(RECEIVERS)
CLASSES = ["EOLP_lost", "sustained", "EOLP_gained", "none"]

checks, fails = [], []


def ck(name, ok, detail=""):
    checks.append(name)
    if not ok:
        fails.append((name, detail))


def nearly(a, b, tol=1e-8):
    return bool(np.allclose(a, b, rtol=tol, atol=tol, equal_nan=True))


def nan_free(s):
    return bool(not pd.isna(s).any())


print("loading tables ...")
sh = pd.read_csv(TAB / "geo_gse211630_liana_signal_shift.csv", encoding="utf-8-sig")
hl = pd.read_csv(TAB / "geo_gse211630_liana_headlines.csv", encoding="utf-8-sig")
inf = pd.read_csv(TAB / "geo_gse211630_liana_incoming_focus.csv", encoding="utf-8-sig")
cr = pd.read_csv(TAB / "geo_gse211630_liana_receptor_module_corr.csv", encoding="utf-8-sig")
fg = pd.read_csv(TAB / "geo_gse211630_liana_figd1_summary.csv", encoding="utf-8-sig")
fl = pd.read_csv(TAB / "geo_gse211630_liana_figd1_summary_ligands.csv", encoding="utf-8-sig")

# ============================ 1. signal_shift (merged) ============================
ck("shift: n rows == 43,495", len(sh) == 43495, str(len(sh)))
req = ["source", "target", "ligand", "receptor", "source_family", "target_cluster",
       "receiver_focus", "present_N", "present_E", "prob_N", "prob_E", "p_N", "p_E",
       "sig_N", "sig_E", "sw_N", "sw_E", "rank_pct_N", "rank_pct_E",
       "delta_prob", "delta_sw", "class"]
ck("shift: required columns present", all(c in sh.columns for c in req))
ck("shift: no duplicated L-R rows",
   not sh.duplicated(subset=["source", "target", "ligand", "receptor"]).any())
ck("shift: prob in [0,1]", float(sh[["prob_N", "prob_E"]].min().min()) >= 0
   and float(sh[["prob_N", "prob_E"]].max().max()) <= 1)
ck("shift: sw finite & no NaN (specificity weight, unbounded)",
   nan_free(sh["sw_N"]) and nan_free(sh["sw_E"])
   and bool(np.isfinite(sh["sw_N"]).all()) and bool(np.isfinite(sh["sw_E"]).all()))
ck("shift: p in [0,1]", float(sh[["p_N", "p_E"]].min().min()) >= 0
   and float(sh[["p_N", "p_E"]].max().max()) <= 1)
ck("shift: delta_prob == prob_E - prob_N",
   nearly(sh["delta_prob"], sh["prob_E"] - sh["prob_N"]))
ck("shift: delta_sw == sw_E - sw_N",
   nearly(sh["delta_sw"], sh["sw_E"] - sh["sw_N"]))
n_foc = int(sh["receiver_focus"].sum())
ck("shift: receiver_focus rows == 20,396", n_foc == 20396, str(n_foc))
got_t = set(sh.loc[sh["receiver_focus"], "target_cluster"].unique())
ck("shift: focus targets == {0,1,3,4,5,6,7,8,16}", got_t == FOCUS, str(got_t))

cls_ok = {
    "sustained":   (True, True, True, True),
    "EOLP_lost":   (True, None, True, False),
    "EOLP_gained": (None, True, False, True),
    "none":        (None, None, False, False),
}
for cl, (pn, pe, sn, se) in cls_ok.items():
    sub = sh[sh["class"] == cl]
    rule = (
        (sub["sig_N"].eq(sn) if sn is not None else True) &
        (sub["sig_E"].eq(se) if se is not None else True) &
        (sub["present_N"].eq(pn) if pn is not None else True) &
        (sub["present_E"].eq(pe) if pe is not None else True)
    )
    ck(f"shift: class '{cl}' fully consistent with sig/present flags",
       int((~rule).sum()) == 0, f"viol={int((~rule).sum())}")
ck("shift: p_N NaN exactly when present_N False (whole table)",
   bool((sh["p_N"].isna() == (sh["present_N"] == False)).all()))
ck("shift: p_E NaN exactly when present_E False (whole table)",
   bool((sh["p_E"].isna() == (sh["present_E"] == False)).all()))
ck("shift: prob_N == 0 whenever absent in NEOLP",
   bool((sh.loc[sh["present_N"] == False, "prob_N"] == 0).all()))
ck("shift: prob_E == 0 whenever absent in EOLP",
   bool((sh.loc[sh["present_E"] == False, "prob_E"] == 0).all()))
ck("shift: class values restricted", set(sh["class"].unique()) <= set(CLASSES))
n_sus_tnk = int(((sh["class"] == "sustained") & sh["target_cluster"].isin(TN_FOCUS)).sum())
ck("shift: sustained T/NK pairs == 2,381", n_sus_tnk == 2381, str(n_sus_tnk))
n_sus_foc = int(((sh["class"] == "sustained") & sh["receiver_focus"]).sum())
ck("shift: sustained focus pairs == 4,033", n_sus_foc == 4033, str(n_sus_foc))

# ======================= 2. headlines (40 presence headlines) =====================
ck("headlines: n rows == 40", len(hl) == 40, str(len(hl)))
ck("headlines: only EOLP_gained / EOLP_lost classes",
   set(hl["class"].unique()) == {"EOLP_gained", "EOLP_lost"},
   str(hl["class"].value_counts().to_dict()))
ck("headlines: 20 gained / 20 lost",
   int((hl["class"] == "EOLP_gained").sum()) == 20
   and int((hl["class"] == "EOLP_lost").sum()) == 20)
ck("headlines: all rows receiver_focus", bool(hl["receiver_focus"].all()))
key = ["source", "target", "ligand", "receptor"]
hl2 = hl.merge(sh[key + ["class", "prob_N", "prob_E", "delta_prob", "delta_sw"]],
               on=key, how="left", suffixes=("", "_sh"))
ck("headlines: every row joins signal_shift", int(hl2["class_sh"].isna().sum()) == 0)
ck("headlines: class identical to signal_shift",
   bool((hl2["class"] == hl2["class_sh"]).all()))
ck("headlines: prob_E identical to signal_shift",
   nearly(hl2["prob_E"], hl2["prob_E_sh"]))
ck("headlines: delta_prob identical to signal_shift",
   nearly(hl2["delta_prob"], hl2["delta_prob_sh"]))
g_sig = hl[hl["class"] == "EOLP_gained"]
l_sig = hl[hl["class"] == "EOLP_lost"]
ck("headlines: gained rows sig_E & ~sig_N, lost rows sig_N & ~sig_E",
   bool(g_sig["sig_E"].all()) and bool((~g_sig["sig_N"]).all())
   and bool(l_sig["sig_N"].all()) and bool((~l_sig["sig_E"]).all()))

# ==================== 3. incoming_focus (receiver x source x class) ==============
ck("inf: n rows == 283", len(inf) == 283, str(len(inf)))
ck("inf: targets restricted to focus receivers",
   set(inf["target_cluster"].unique()) <= FOCUS)
ck("inf: classes restricted", set(inf["class"].unique()) <= set(CLASSES))
for t in RECEIVERS:  # cross-check with signal_shift per target
    n_sh = int((sh["receiver_focus"] & (sh["target_cluster"] == t)).sum())
    n_inf = int(inf.loc[inf["target_cluster"] == t, "n"].sum())
    ck(f"inf: c{t} total n == signal_shift focus rows ({n_sh})", n_sh == n_inf,
       f"inf={n_inf}")
for cl in ["EOLP_lost", "sustained", "EOLP_gained", "none"]:
    for t in RECEIVERS:
        v_sh = sh.loc[(sh["class"] == cl) & sh["receiver_focus"]
                      & (sh["target_cluster"] == t), "delta_prob"].sum()
        v_inf = inf.loc[(inf["class"] == cl) & (inf["target_cluster"] == t),
                        "sum_delta_prob"].sum()
        ck(f"inf: sum_delta_prob c{t} {cl} == signal_shift ({v_sh:.4f})",
           nearly(v_sh, v_inf), f"inf={v_inf:.4f}")
        w_sh = sh.loc[(sh["class"] == cl) & sh["receiver_focus"]
                      & (sh["target_cluster"] == t), "delta_sw"].sum()
        w_inf = inf.loc[(inf["class"] == cl) & (inf["target_cluster"] == t),
                        "sum_delta_sw"].sum()
        ck(f"inf: sum_delta_sw c{t} {cl} == signal_shift ({w_sh:.4f})",
           nearly(w_sh, w_inf), f"inf={w_inf:.4f}")
pres = inf.groupby(["target_cluster", "class"])["n"].sum().unstack(fill_value=0)
for t in RECEIVERS:
    nf = int(pres.loc[t, "EOLP_gained"] - pres.loc[t, "EOLP_lost"])
    ck(f"inf: net_flux c{t} arithmetic (= {nf:+d})",
       nf == int(fg.loc[fg["receiver"] == t, "net_flux"].iloc[0]))

# ==================== 4. receptor-module correlation =============================
ck("corr: n rows == 392", len(cr) == 392, str(len(cr)))
ck("corr: rho in [-1,1]", float(cr["rho"].min()) >= -1 and float(cr["rho"].max()) <= 1)
ck("corr: p in [0,1]", float(cr["p"].min()) >= 0 and float(cr["p"].max()) <= 1)
ck("corr: n_incoming_interactions integer >= 0",
   bool((cr["n_incoming_interactions"] >= 0).all())
   and bool(np.equal(cr["n_incoming_interactions"],
                     cr["n_incoming_interactions"].astype(int)).all()))
ck("corr: no NaN in rho/p", nan_free(cr["rho"]) and nan_free(cr["p"]))
ck("corr: no NaN in n_incoming_interactions", nan_free(cr["n_incoming_interactions"]))
# modules match Stage B module names
mod_exp = {"Cytotoxic", "Exhaustion", "Treg", "Naive", "Stress", "Activation",
           "Proliferation"}
ck("corr: module names within Stage B set", set(cr["module"].unique()) <= mod_exp,
   str(cr["module"].unique()))
per_rec = cr.groupby("receptor")["module"].size()
ck("corr: exactly 56 unique receptors x 7 modules (1 row per pair)",
   len(per_rec) == 56 and int((per_rec == 7).sum()) == 56
   and int(cr.duplicated(subset=["receptor", "module"]).sum()) == 0,
   f"n_rec={len(per_rec)} n_7={int((per_rec == 7).sum())}")

# ==================== 5. Fig D1 canonical summary CSVs ==========================
ck("figd1: 9 receiver rows in canonical order",
   list(fg["receiver"]) == RECEIVERS)
st = sh.loc[sh["class"] == "sustained"].groupby("target_cluster")["delta_prob"]
sus_focus = sh.loc[(sh["class"] == "sustained") & sh["receiver_focus"]]
for t in RECEIVERS:
    n_ = int((sus_focus["target_cluster"] == t).sum())
    s_ = float(st.get_group(t).sum())
    r_ = fg.loc[fg["receiver"] == t].iloc[0]
    ck(f"figd1: c{t} n_sustained {n_} == summary", n_ == int(r_["n_sustained_sig_pairs"]),
       f"summary={int(r_['n_sustained_sig_pairs'])}")
    ck(f"figd1: c{t} sum_delta_prob {s_:.4f} == summary",
       nearly(s_, float(r_["sum_delta_prob_sustained"])))
    ck(f"figd1: c{t} focus_TNK flag correct", bool(r_["focus_TNK"]) == (t in TN_FOCUS))
    ck(f"figd1: c{t} presence counts == incoming_focus",
       int(r_["pairs_EOLP_lost"]) == int(pres.loc[t, "EOLP_lost"])
       and int(r_["pairs_sustained"]) == int(pres.loc[t, "sustained"])
       and int(r_["pairs_EOLP_gained"]) == int(pres.loc[t, "EOLP_gained"]))
# ligand summary
lt = sh.loc[(sh["class"] == "sustained") & sh["target_cluster"].isin(TN_FOCUS)]
ltg = lt.groupby("ligand")["delta_prob"].agg(n_pairs="size", s="sum").query("n_pairs >= 14")
ck("figd1-lig: 10 displayed ligands (6 most negative reversed + 4 most positive, n>=14)",
   len(fl) == 10 and
   list(fl["ligand"]) == (list(ltg.sort_values("s").head(6).iloc[::-1].index)
                          + list(ltg.sort_values("s").tail(4).index)),
   str(list(fl["ligand"])))
for _, rw in fl.iterrows():
    row = ltg.loc[rw["ligand"]]
    ck(f"figd1-lig: {rw['ligand']} n {int(row['n_pairs'])} == summary",
       int(row["n_pairs"]) == int(rw["n_sustained_pairs_TNK"]))
    ck(f"figd1-lig: {rw['ligand']} sum {row['s']:.4f} == summary",
       nearly(float(row["s"]), float(rw["sum_delta_prob_sustained_TNK"])))
mhc1 = float(ltg.loc[["HLA-A", "HLA-C", "HLA-E"], "s"].sum())
ck("figd1-lig: MHC-I (HLA-A/C/E) sustained sum toward T/NK == -6.199",
   nearly(mhc1, -6.1985, tol=1e-3), f"{mhc1:.4f}")
# canonical headline numbers for narrative
ck("figd1: T/NK sustained total 2,381 reproduced",
   int(fg.loc[fg["focus_TNK"], "n_sustained_sig_pairs"].sum()) == 2381)
ck("figd1: all seven T/NK receivers negative sum_delta_prob",
   bool((fg.loc[fg["focus_TNK"], "sum_delta_prob_sustained"] < 0).all()))
pos_rec = list(fg.loc[fg["net_flux"] > 0, "receiver"])
neg_rec = list(fg.loc[fg["net_flux"] < 0, "receiver"])
ck("figd1: presence-flux net positive receivers == {0,1,3,6}",
   sorted(pos_rec) == [0, 1, 3, 6], str(pos_rec))
ck("figd1: presence-flux net negative receivers == {4,5,7,8,16}",
   sorted(neg_rec) == [4, 5, 7, 8, 16], str(neg_rec))
ck("figd1: total net flux across receivers == +298",
   int(fg["net_flux"].sum()) == 298, str(int(fg["net_flux"].sum())))

# ================================ report =========================================
OUT.parent.mkdir(parents=True, exist_ok=True)
lines = [
    "Stage D table verification (geo_gse211630 liana / cell-cell communication)\n",
    "scope : signal_shift | headlines | incoming_focus | receptor_module_corr | figd1 summary",
    "conventions : class from sig flags; focus receivers {0,1,3,4,5,6,7,8,16}; T/NK focus {3..8,16}",
    f"result : {len(checks)} checks, {len(checks) - len(fails)} PASS, {len(fails)} FAIL\n",
]
for name, _ in fails:
    lines.append(f"FAIL  {name}")
lines.append(f"\nall checks passed" if not fails else f"\n{len(fails)} failing checks")
OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"Stage D: {len(checks) - len(fails)}/{len(checks)} checks passed -> {OUT}")
if fails:
    for n, d in fails:
        print("  FAIL", n, "|", d[:120])
    sys.exit(1)
