# -*- coding: utf-8 -*-
"""Build results/reports/reference_pmid_list.md from the verified PMID lookup.

Reads results/tables/_ref_pmid_lookup.json (produced by _fetch_ref_pmids.py via
NCBI E-utilities esearch/esummary; all 34 entries re-resolved on 2026-09-13 after
the alphabetical insertions of Walker 2016, Machelak 2023, Yang 2025 and the four database
publications (Buniello 2025, Cannon 2024, Kelleher 2023, Kurki 2023); every PMID
title- or DOI-verified, then re-verified once more on 2026-09-13 after the
NLM renumbering) and emits a Markdown reference list in the manuscript's
first-citation order (1-34), each entry identified by its PubMed ID, plus a
plain-PMID block for batch import (EndNote, PubMed "find by ID", submission systems).

Author names use the real accented spellings to match manuscript_main_text.md.
"""
import json
import os

SRC = os.path.join("results", "tables", "_ref_pmid_lookup.json")
OUT = os.path.join("results", "reports", "reference_pmid_list.md")

# n -> short identifying tag exactly as in the manuscript list (authors. Journal. Year)
TAGS = {
    1: "Eisen D, Carrozzo M, Bagan Sebastian JV, Thongprasom K. Oral Dis. 2005",
    2: "McCartan BE, Healy CM. J Oral Pathol Med. 2008",
    3: "González-Moles MA, Ramos-García P, Warnakulasuriya S. Oral Dis. 2021",
    4: "Warnakulasuriya S, Johnson NW, van der Waal I. J Oral Pathol Med. 2007",
    5: "Akbari E, Epstein JB, Samim F. Head Neck. 2025",
    6: "Lodi G, Manfredi M, Mercadante V, Murphy R, Carrozzo M. Cochrane Database Syst Rev. 2020",
    7: "Ritari J, et al. HLA. 2026",
    8: "Sugerman PB, et al. Crit Rev Oral Biol Med. 2002",
    9: "Davey Smith G, Ebrahim S. Int J Epidemiol. 2003",
    10: "Schmidt AF, et al. Nat Commun. 2020",
    11: "Cao Y, et al. J Transl Med. 2023",
    12: "Qin Y, et al. Invest Ophthalmol Vis Sci. 2024",
    13: "Alayash Z, et al. J Clin Periodontol. 2024",
    14: "Liu L, et al. J Dent Res. 2024",
    15: "Georgakis MK, Gill D. Stroke. 2021",
    16: "Finan C, et al. Sci Transl Med. 2017",
    17: "Qing M, et al. eLife. 2023",
    18: "Skrivankova VW, et al. JAMA. 2021",
    19: "Kelleher KJ, et al. Nucleic Acids Res. 2023",
    20: "Cannon M, et al. Nucleic Acids Res. 2024",
    21: "Buniello A, et al. Nucleic Acids Res. 2025",
    22: "Võsa U, et al. Nat Genet. 2021",
    23: "Kurki MI, et al. Nature. 2023",
    24: "1000 Genomes Project Consortium (Auton A, et al.). Nature. 2015",
    25: "Sun BB, et al. Nature. 2018",
    26: "Yang C, et al. Nat Commun. 2025",
    27: "Eldjarn GH, et al. Nature. 2023",
    28: "Giambartolomei C, et al. PLoS Genet. 2014",
    29: "Wallace C. PLoS Genet. 2021",
    30: "Wang G, Sarkar A, Carbonetto P, Stephens M. J R Stat Soc B. 2020",
    31: "Hemani G, Tilling K, Davey Smith G. PLoS Genet. 2017",
    32: "Trott O, Olson AJ. J Comput Chem. 2010",
    33: "Walker RG, et al. Circ Res. 2016",
    34: "Machelak W, Szczepaniak A, Jacenik D, Zielińska M. Life Sci. 2023",
}


recs = json.load(open(SRC, encoding="utf-8"))
by_n = {r["n"]: r for r in recs}
assert list(by_n) == list(range(1, 35))
assert all(by_n[n].get("pmid") for n in by_n)

lines = []
lines.append("# References with PubMed IDs — OLP drug-target MR (FinnGen R13)")
lines.append("")
lines.append("> 主稿参考文献列表（`results/reports/manuscript_main_text.md`，编号 1–34，"
             "按正文首次引用顺序排列，NLM/ICMJE 编号制）"
             "的 PubMed ID 清单，供按 PMID 在 EndNote / PubMed（Find by ID）/ 投稿系统批量添加文献。")
lines.append("")
lines.append("**PMID 获取与校验**：全表于 2026-09-13 重新经 NCBI E-utilities（esearch 按标题精确检索 → "
             "esummary 核对）逐条解析（此前 26 条于 2026-09-06 首次解析，当时为字母序编号）；"
             "命中标准为 DOI 精确匹配或标题相似度 ≥ 0.84（全部命中，34/34）。"
             "2026-09-13 引文改为 NLM 编号制后，本表随参考文献表按正文首次引用顺序重排并重编号，"
             "PMID 于同日晚按重排后编号再次逐条 esummary 复核（34/34，编号与姓氏一一对应）。"
             "新增条目——19 Kelleher 2023（Pharos）、23 Kurki 2023"
             "（FinnGen）、34 Machelak 2023、33 Walker 2016、26 Yang 2025、以及 21 Buniello 2025"
             "（Open Targets Platform）、20 Cannon 2024（DGIdb 5.0）——另经 Crossref 与 Europe PMC 复核。"
             "未使用任何记忆值。"
             "校验台账：`results/tables/_ref_pmid_lookup.json`；"
             "完整核验证据：`results/reports/reference_verification_report.md`。")
lines.append("")
lines.append("| # | Reference (as in manuscript) | PMID |")
lines.append("|---|------------------------------|------|")
for n in range(1, 35):
    pmid = by_n[n]["pmid"]
    lines.append(f"| {n} | {TAGS[n]} | [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) |")
lines.append("")
lines.append("## Plain PMID list (batch import)")
lines.append("")
lines.append("```")
for n in range(1, 35):
    lines.append(str(by_n[n]["pmid"]))
lines.append("```")
lines.append("")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(lines))
print("written:", OUT, "| entries:", sum(1 for l in lines if l.startswith("| ") and not l.startswith("| #")))
