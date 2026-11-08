# -*- coding: utf-8 -*-
"""Emit the audit record for the NLM citation renumbering (2026-09-13).

Prints a file/MD5 table and the old->new citation-number mapping, both read from disk so the
audit note cannot drift from the artefacts. Idempotent: prints only, writes nothing.
"""
import hashlib
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FILES = [
    ("主稿", "results/reports/manuscript_main_text.md"),
    ("图注（主图）", "figures/figure_legends_main.md"),
    ("图清单（主图）", "figures/figure_manifest.md"),
    ("核验报告", "results/reports/reference_verification_report.md"),
    ("编号对照表", "results/tables/_citation_number_map.json"),
    ("PMID 台账", "results/tables/_ref_pmid_lookup.json"),
    ("PMID 清单", "results/reports/reference_pmid_list.md"),
    ("RIS", "results/reports/reference_pmid_list.ris"),
    ("转换脚本", "analysis/scripts/convert_citations_to_nlm.py"),
    ("QA：引文", "analysis/scripts/qa_citations_nlm.py"),
    ("QA：DOCX 引文", "analysis/scripts/qa_docx_citations.py"),
    ("PMID 抓取", "analysis/scripts/_fetch_ref_pmids.py"),
    ("PMID 清单构建", "analysis/scripts/build_reference_pmid_list.py"),
    ("RIS 构建", "analysis/scripts/build_reference_ris.py"),
    ("docx：打印版", "results/reports/OLP_MR_manuscript_print.docx"),
    ("docx：审查合集", "results/reports/OLP_MR_review_compilation.docx"),
    ("docx：图表审查包", "results/reports/OLP_MR_figures_tables_review.docx"),
]

print("| 文件 | 字节 | MD5(前16) |")
print("|---|---:|---|")
for label, rel in FILES:
    path = os.path.join(ROOT, rel)
    if not os.path.exists(path):
        print(f"| {label} `{rel}` | MISSING | - |")
        continue
    b = open(path, "rb").read()
    print(f"| {label} `{rel}` | {len(b):,} | `{hashlib.md5(b).hexdigest()[:16].upper()}` |")

mp = json.load(open(os.path.join(ROOT, "results/tables/_citation_number_map.json"),
                    encoding="utf-8-sig"))
o2n = {int(k): v for k, v in mp["old_number_to_new"].items()}
s2n = mp["surname_to_number"]
new_to_surname = {n: s for s, n in s2n.items()}
print("\n| 旧（字母序） | 新（首现序） | 姓氏 |")
print("|---:|---:|---|")
for old in sorted(o2n):
    new = o2n[old]
    print(f"| {old} | {new} | {new_to_surname[new]} |")
