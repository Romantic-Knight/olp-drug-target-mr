"""Build two Word deliverables for the OLP drug-target MR manuscript.

1) OLP_MR_manuscript_print.docx        - clean printable main text (author print
                                          revision copy: title/authors/abstract/
                                          body/references + 1-page author checklist).
2) OLP_MR_figures_tables_review.docx   - figure & table review pack: main figures
                                          1-7 and supplementary figures 1-10 with
                                          their official legends and embedded PNGs,
                                          plus supplementary tables 1-10 rendered
                                          as Word tables (full for the small- and
                                          medium-sized tables; column structure +
                                          first 8 rows preview + full CSV path for
                                          the two large tables S2 and S5), a
                                          numerical verification-notes appendix and
                                          abbreviations.

Scope note (2026-09-23, after the supplementary-table trim): Part III of the review pack renders
the five on-request table legends, in first-citation order, together with the tables they describe
(Tables 1 and 3 - the 5,040-row screen and the 756-row credible-set list - as a preview plus CSV
path; Tables 2, 4 and 5 in full). The ten former tables that were retired because a figure already
displays their content are listed afterwards as numberless source data. Supplementary tables are
NOT part of the submission (author decision: supplementary material is figures-only, tables
available from the corresponding author on request), so this pack doubles as the on-request
delivery document; the STROBE-MR checklist is the one tabular item that accompanies the submission
as a separate file.

Uses the rendering helpers of build_review_docx.py (same fonts/layout language).

NOTE (2026-09-13): numbering systems in the manuscript are the CURRENT ones chosen by
the authors: in-text citations are NLM/ICMJE bracketed numerals in first-citation
order and the reference list is numbered 1-34 in that same order, with
[verified ...] annotations, supplementary figures S1-S10 (the STROBE-MR checklist
carries no "Supplementary Table" number because it is the only tabular item
submitted). No renumbering is performed here.
"""
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_review_docx as brv  # noqa: E402  (helpers: add_heading, add_par, ...)

from docx import Document  # noqa: E402
from docx.shared import Pt, Inches, RGBColor  # noqa: E402
from docx.enum.section import WD_SECTION, WD_ORIENT  # noqa: E402
from docx.oxml.ns import qn  # noqa: E402

# ----------------------------------------------------------------------------- paths
ROOT = brv.ROOT
REPORTS = brv.REPORTS
MAIN_MD = os.path.join(REPORTS, "manuscript_main_text.md")
MAIN_LEGEND = os.path.join(ROOT, "figures", "main_figures_v2", "figure_legends_main.md")
SUPP_LEGEND = os.path.join(ROOT, "figures", "supplementary_v2", "figure_legends_supplementary.md")
# The 2026-09-30 restructure (6 main + 3 supplementary) lives in the _v2 folders; the original
# folders keep the previous 7 + 10 set untouched.
MAIN_FIG_DIR = os.path.join(ROOT, "figures", "main_figures_v2")
SUPP_FIG_DIR = os.path.join(ROOT, "figures", "supplementary_v2")
SUPP_TAB_DIR = os.path.join(ROOT, "results", "tables", "supplementary")
OUT_PRINT = os.path.join(REPORTS, "OLP_MR_manuscript_print.docx")
OUT_REVIEW = os.path.join(REPORTS, "OLP_MR_figures_tables_review.docx")

# Final title (round-4 scRNA restructure, 2026-09-13); kept identical to the H1 of
# results/reports/manuscript_main_text.md and to results/reports/manuscript_title_page.md.
TITLE = ("From germline risk to lesional tissue: druggable-genome Mendelian randomization "
         "and single-cell profiling nominate GDF11 inhibition for oral lichen planus")

AUTHORS_HTML = ("Yahui Lu<sup>1</sup>, Xiteng Yin<sup>1</sup>, Jianchuan Ran<sup>1</sup>, "
                 "Wenguang Xu<sup>1</sup>, Yufeng Wang<sup>1</sup>, Wei Han<sup>1,*</sup>")
AFFIL_HTML = ("<sup>1</sup> Nanjing Stomatological Hospital, Affiliated Hospital of Medical School, "
              "Institute of Stomatology, Nanjing University, Nanjing 210008, China.")
CORRESP_HTML = ("*Correspondence:* Wei Han, Nanjing Stomatological Hospital, Affiliated Hospital of "
                "Medical School, Institute of Stomatology, Nanjing University, Nanjing 210008, China. "
                "Email: doctorhanwei@nju.edu.cn")

MAIN_FIG_MAP = {
    "Figure 1.": "fig01_design_screen.png",
    "Figure 2.": "fig02_replication_colocalisation.png",
    "Figure 3.": "fig03_pleiotropy_docking.png",
    "Figure 4.": "fig04_docking_md.png",
    "Figure 5.": "fig05_atlas_composition.png",
    "Figure 6.": "fig06_candidate_expression.png",
    "Figure 7.": "fig07_target_mapping.png",
}
SUPP_FIG_MAP = {
    "Supplementary Figure 1.": "suppfig01_regional_plots.png",
    "Supplementary Figure 2.": "suppfig02_dataset_qc.png",
    "Supplementary Figure 3.": "suppfig03_detection_composition.png",
    "Supplementary Figure 4.": "suppfig04_candidate_landscape.png",
}
TABLE_FILES = {
    1: "suppl_table_01_full_mr_screen_5040.csv",
    2: "suppl_table_02_wide_hits_mr_replication_coloc.csv",
    3: "suppl_table_03_susie_credible_sets.csv",
    4: "suppl_table_04_phewas_gw_hits.csv",
    5: "suppl_table_05_pqtl_assessment.csv",
}
STROBE_FILE = "strobe_mr_checklist.csv"
# Retired on 2026-09-23 (former Tables 4/7/8/9 plus the five GSE211630 single-cell summary
# tables): their content is already displayed in a figure, so they carry no table number.
# They are listed in Part III as source data so the pack still inventories every deliverable.
RETIRED_SOURCES = [
    ("results/tables/supplementary/_retired/coloc_abf_susie_44runs.csv",
     "Former Table 4 (44 runs) - the same PP.H0-H4 posteriors are plotted in Supplementary Figure 3"),
    ("results/tables/supplementary/_retired/pleiotropy_verdict_5x9.csv",
     "Former Table 7 (45 runs) - Figure 5D prints each abf PP.H4 cell; Figure 5C gives the SuSiE verdicts"),
    ("results/tables/supplementary/_retired/docking_results.csv",
     "Former Table 8 (22 ligands) - Figure 3B labels name and affinity of every passing ligand"),
    ("results/tables/supplementary/_retired/steiger_direction.csv",
     "Former Table 9 (50 runs) - Figure 5B covers the OLP WIDE endpoint for all five instruments"),
    ("results/tables/geo_gse211630_erosion_candidate_celltype.csv + _erosion_per_sample.csv",
     "Former Table 11 - candidate-gene expression by cell type (Figure 6, Supplementary Figure 8)"),
    ("results/tables/geo_gse211630_composition_{per_sample,shift}.csv",
     "Former Table 12 - composition shift (Supplementary Figure 9A/B)"),
    ("results/tables/geo_gse211630_tnk_module_{per_sample,shift}.csv",
     "Former Table 13 - T/NK module scores (Supplementary Figure 9C)"),
    ("results/tables/geo_gse211630_mr_backmap_{gene,subtype}.csv",
     "Former Table 14 - tissue-level reverse mapping and Oaxaca terms (Figure 7)"),
    ("results/tables/geo_gse211630_liana_*.csv",
     "Former Table 15 - LIANA cell-cell communication (Main Figure 6D-F)"),
]
BIG_TABLES = {1, 3}
PREVIEW_ROWS = 8
GRAY = RGBColor(0x66, 0x66, 0x66)

# ----------------------------------------------------------------------------- document scaffolding


def new_doc(header_text, landscape=False):
    doc = Document()
    sec = doc.sections[0]
    if landscape:
        sec.orientation = WD_ORIENT.LANDSCAPE
        sec.page_width = brv.PAGE_H
        sec.page_height = brv.PAGE_W
        sec.left_margin = sec.right_margin = Inches(0.5)
        sec.top_margin = sec.bottom_margin = Inches(0.6)
    else:
        sec.page_width = brv.PAGE_W
        sec.page_height = brv.PAGE_H
        sec.left_margin = sec.right_margin = brv.MARGIN
        sec.top_margin = sec.bottom_margin = Inches(0.79)
    st = doc.styles['Normal']
    st.font.name = brv.BODY
    st.font.size = Pt(11)
    st.element.rPr.rFonts.set(qn('w:eastAsia'), brv.CJK)
    brv.add_header_footer(doc, header_text)
    return doc


def note_par(doc, text, size=9.5):
    brv.add_par(doc, text, size=size, italic=True, color=GRAY)


def csv_rows(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.reader(f))


def clean_cell(v):
    v = (v or "").strip()
    if v.lower() == "nan":
        return ""
    return v


def csv_column_blocks(n_cols, key_cols=1, per_block=8):
    """[(start, end)] column slices for readable rendering of a wide CSV.

    Word gives every column of a table the same width, so a 29-column table rendered as one
    block collapses each cell into one character per line. A landscape A4 page fits roughly
    eight usable columns at this font size, so wider tables are emitted as horizontal blocks
    that repeat the key column(s) (e.g. `gene`): every row is preserved and no cell is
    squeezed. Returns [] only for an empty CSV.
    """
    if n_cols <= 0:
        return []
    span = max(1, per_block - key_cols)
    return [(k, min(k + span, n_cols)) for k in range(key_cols, n_cols, span)] or \
           [(min(key_cols, n_cols), n_cols)]


def add_csv_table_blocks(doc, path, max_rows=None, key_cols=1, per_block=8,
                         cell_size=7.5, hdr_size=8.0):
    """CSV -> one or more Word tables split with `csv_column_blocks`. Returns table count."""
    rows = csv_rows(path)
    if not rows:
        brv.add_par(doc, "[empty CSV]", italic=True, color=GRAY)
        return 0
    hdr = [clean_cell(h) for h in rows[0]]
    body = [r for r in rows[1:] if any(clean_cell(c) for c in r)]
    if max_rows is not None:
        body = body[:max_rows]
    blocks = csv_column_blocks(len(hdr), key_cols=key_cols, per_block=per_block)
    for start, end in blocks:
        idx = list(range(0, min(key_cols, len(hdr)))) + list(range(start, end))
        if len(blocks) > 1:
            note_par(doc, f"columns {start + 1}-{end} of {len(hdr)} "
                          f"(key column repeated): " +
                          ", ".join(hdr[i] for i in range(start, end)))
        tbl = doc.add_table(rows=1, cols=len(idx))
        tbl.style = 'Table Grid'
        tbl.autofit = True
        for j, i in enumerate(idx):
            cell = tbl.rows[0].cells[j]
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            r = p.add_run(hdr[i])
            brv.set_run_font(r, brv.BODY, hdr_size, bold=True)
            brv._shade(cell, "E8E8E8")
        for row in body:
            cells = tbl.add_row().cells
            for j, i in enumerate(idx):
                v = clean_cell(row[i]) if i < len(row) else ""
                cell = cells[j]
                cell.text = ""
                p = cell.paragraphs[0]
                p.paragraph_format.space_after = Pt(0)
                r = p.add_run(v)
                brv.set_run_font(r, brv.BODY, cell_size)
        brv.add_par(doc, "", size=2)
    return len(blocks)


def parse_legend_blocks(path):
    """[(heading, [(kind,payload), ...])] from a legend markdown (h2 blocks)."""
    blocks = []
    cur = None
    for kind, payload in brv.iter_md(path):
        if kind == "h2":
            if cur:
                blocks.append(cur)
            cur = [payload, []]
        elif kind in ("h1", "h3", "hr", "blank"):
            continue
        else:
            if cur is not None:
                cur[1].append((kind, payload))
    if cur:
        blocks.append(cur)
    return blocks


def render_block_paras(doc, paras, size=10.5):
    for kind, payload in paras:
        if kind == "text":
            brv.add_par(doc, payload, size=size)
        elif kind == "bullet":
            brv.add_bullet(doc, payload, size=size)
        elif kind == "num":
            brv.add_numbered(doc, payload, size=size)
        else:
            brv.add_par(doc, payload, size=size)


def section_to_landscape(doc):
    sec = doc.add_section(WD_SECTION.NEW_PAGE)
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width = brv.PAGE_H
    sec.page_height = brv.PAGE_W
    sec.left_margin = sec.right_margin = Inches(0.5)
    sec.top_margin = sec.bottom_margin = Inches(0.6)
    return sec


def section_to_portrait(doc):
    sec = doc.add_section(WD_SECTION.NEW_PAGE)
    sec.orientation = WD_ORIENT.PORTRAIT
    sec.page_width = brv.PAGE_W
    sec.page_height = brv.PAGE_H
    sec.left_margin = sec.right_margin = brv.MARGIN
    sec.top_margin = sec.bottom_margin = Inches(0.79)
    return sec


def collect_checks(main_md_text):
    out = []
    for m in re.finditer(r'\[[Cc]heck[^\]]*\]|\[To be completed[^\]]*\]', main_md_text):
        s = max(0, m.start() - 95)
        ctx = re.sub(r"\s+", " ", main_md_text[s:m.end() + 25]).strip()
        out.append((m.group(0), ctx))
    seen, uniq = set(), []
    for g, ctx in out:
        if ctx in seen:
            continue
        seen.add(ctx)
        uniq.append((g, ctx))
    return uniq


# =============================================================================
# Deliverable 1 - printable main text
# =============================================================================
def build_print():
    doc = new_doc("")   # print copy of the manuscript: no running header

    # --- title block -----------------------------------------------------
    brv.add_heading(doc, TITLE, 1)
    brv.add_par(doc, AUTHORS_HTML, size=12, align=brv.WD_ALIGN_PARAGRAPH.CENTER)
    brv.add_par(doc, AFFIL_HTML, size=10.5, align=brv.WD_ALIGN_PARAGRAPH.CENTER)
    brv.add_par(doc, CORRESP_HTML, size=10, italic=True,
                align=brv.WD_ALIGN_PARAGRAPH.CENTER, color=GRAY)
    brv.add_hr(doc)

    # --- main body, HTML comments stripped --------------------------------
    raw = open(MAIN_MD, encoding="utf-8").read()
    clean = re.sub(r"<!--.*?-->", "", raw, flags=re.S)
    clean = re.sub(r"\n{3,}", "\n\n", clean)
    tmp = os.path.join(REPORTS, "_tmp_main_clean.md")
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        f.write(clean)
    try:
        brv.render_md(doc, tmp, stop_heading="References", skip_first_h1=True)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    # --- references (1-34), stop before the trailing draft footnote -------
    brv.add_heading(doc, "References", 2)
    for kind, payload in brv.iter_md(MAIN_MD):
        if kind == "h2" and payload.startswith("References"):
            break
    started = False
    for kind, payload in brv.iter_md(MAIN_MD):
        if kind == "h2" and payload.startswith("References"):
            started = True
            continue
        if not started:
            continue
        if payload.startswith("*End of draft") or (kind == "hr"):
            break
        if kind in ("blank",):
            continue
        if kind == "num":
            brv.add_numbered(doc, payload, size=10.5)
        elif kind == "text":
            brv.add_par(doc, payload, size=10.5, italic=True)
        else:
            brv.add_par(doc, payload, size=10.5)
    # trailing note that the references carry verification annotations
    note_par(doc, ("Note: each entry carries a [verified ...] annotation "
                   "(Crossref/PubMed evidence ledger in "
                   "results/reports/reference_verification_report.md); remove before submission "
                   "if the target journal forbids such annotations."))

    # --- author checklist page --------------------------------------------
    brv.add_page_break(doc)
    brv.add_heading(doc, "Author checklist for this print copy (not for journal)", 1)
    brv.add_par(doc, ("本打印版 = 当前主稿 working copy。文内引用为 NLM/ICMJE 方括号编号制"
                      "（[1]、[3,4]、[11-15]，按正文首次引用顺序），文末参考文献为同一顺序的 "
                      "编号列表 1–34（条目内保留 [verified …] 核验注记；编号 2026-09-13 由 "
                      "author-year 转换而来，故非字母序）；补充材料按 2026-09-13 裁决"
                      "只含补充图（S1–S10），补充表不随稿提交、正文表格引用已全部删除、"
                      "Data availability 改为可向通讯作者索取；随稿的表格项仅 STROBE-MR 报告清单一项。"
                      "2026-09-23 进一步把「按需索取」补充表由 15 张精简为 5 张并按正文首现顺序"
                      "重编号（1 全量筛选 → 2 22 命中 MR/复制/共定位 → 3 SuSiE 可信集 → "
                      "4 PheWAS 全基因组显著 → 5 pQTL 评估），被删的 10 张因内容已被图完整呈现。"
                      "Discussion 已按意见改为纯段落（无 4.1–4.5 小标题），Limitations 已压缩为 4 段流畅行文、无小标题；参考文献表已补齐至 34 条（含 Pharos/DGIdb/Open Targets/FinnGen 四个数据库论文），正文 [check] 占位已清零。"),
                size=10)
    brv.add_heading(doc, "P0 待办（正文 §5 Declarations 与标题页均有占位）", 3)
    for t in [
        "Funding —— 资助来源或 “This research received no specific grant…” 一句",
        "Competing interests —— 无则写 “The authors declare no competing interests.”",
        "Ethics —— 二次分析公开汇总数据豁免表述（STROBE-MR 20d）",
        "CRediT contribution statement —— 逐人勾选（title page md）",
        "目标期刊选定 —— 决定格式/摘要长度/图规格",
    ]:
        brv.add_bullet(doc, t, size=10)
    brv.add_heading(doc, "正文仍含 [check] / [To be completed] 占位处", 3)
    checks = collect_checks(raw)
    for tag, ctx in checks:
        brv.add_bullet(doc, f"`{tag}` — …{ctx}…", size=9.5)
    # Author-facing navigation block: paths are for the author's own workspace and are NOT part of
    # the submitted manuscript (the submission build drops every note of this kind; see
    # build_manuscript_submission_docx.py). Kept in step with the 2026-09-30 restructure, whose
    # figures live in the _v2 folders.
    brv.add_heading(doc, "关键文件（作者自查用；不随投稿）", 3)
    for t in [
        "主稿: results/reports/manuscript_main_text.md",
        "标题页/作者: results/reports/manuscript_title_page.md",
        "主图图注: figures/main_figures_v2/figure_legends_main.md",
        "补充图注: figures/supplementary_v2/figure_legends_supplementary.md",
        "主图/补充图: figures/main_figures_v2/、figures/supplementary_v2/",
        "补充表 CSV: results/tables/supplementary/",
        "补充表目录 README: results/tables/supplementary/README.md",
        "参考文献核验报告: results/reports/reference_verification_report.md",
        "STROBE-MR 审计: results/reports/manuscript_audit_strobe_mr.md",
        "图表重组映射: results/reports/figure_restructure_plan3_mapping.md",
    ]:
        brv.add_bullet(doc, t, size=9.5)
    brv.add_hr(doc)
    brv.add_par(doc, ("※ 图表请另见审查版 Word：OLP_MR_figures_tables_review.docx "
                      "（主图 1–6、补充图 1–4 及各表完整/预览渲染）。"), size=10, italic=True, color=GRAY)

    doc.save(OUT_PRINT)
    return doc


# =============================================================================
# Deliverable 2 - figures & tables review pack
# =============================================================================
def build_review():
    doc = new_doc("OLP drug-target MR (FinnGen R13) - figures & tables review pack")

    # --- cover -------------------------------------------------------------
    brv.add_heading(doc, "Figures and Tables Review Pack", 1)
    brv.add_par(doc, ("OLP drug-target Mendelian randomization (FinnGen R13) — 图表集中审查版。"
                      "编号体系与主稿一致：主图 Figure 1–6（2026-09-30 重组）；补充图 Supplementary Figure 1–4。"
                      "按作者 2026-09-13 裁决，投稿的补充材料只含补充图，补充表及表注"
                      "（Part III）为按需索取内容，不随稿提交；随稿的表格项仅 STROBE-MR 报告清单一项。"
                      "2026-09-23 进一步把按需表由 15 张精简为 5 张并按正文首现顺序重编号。"
                      "图注正文与投稿用版本逐字一致（figures/main_figures_v2/figure_legends_main.md、"
                      "figures/supplementary_v2/figure_legends_supplementary.md）。"), size=10.5)
    brv.add_heading(doc, "Contents", 3)
    for t in [
        "Part I — Main Figures 1–6（正式图注 + 嵌入 PNG）",
        "Part II — Supplementary Figures 1–4（图注 + PNG）",
        "Part III — Supplementary Tables 1–5（不随稿，按需索取；横向页；"
        "Table 1、3 仅列结构 + 前 8 行预览 + 完整 CSV 路径，Table 2、4、5 全部行完整渲染，"
        "宽表按 8 列分块并重复 `gene` 列；末尾附无编号的 STROBE-MR 清单与源数据清单）",
        "Part IV — 数值核验注（中文，不随稿，供作者复核）",
        "Part V — Abbreviations（主图注 / 补充图注两套）",
    ]:
        brv.add_bullet(doc, t, size=10.5)
    note_par(doc, ("PNG 均为可发表分辨率（figures/main_figures_v2/、figures/supplementary_v2/，"
                   "同目录附 PDF 矢量版）。大表完整数据以 CSV 提供，不进 Word。"))
    brv.add_page_break(doc)

    # --- Part I main figures ------------------------------------------------
    brv.add_heading(doc, "Part I — Main Figures 1–7 (legends + images)", 1)
    main_blocks = parse_legend_blocks(MAIN_LEGEND)
    main_figs = [b for b in main_blocks if b[0].startswith("Figure ")]
    abbr_main = [b for b in main_blocks if b[0].startswith("Abbreviations")]
    for k, (head, paras) in enumerate(main_figs):
        key = head.split(".")[0] + "."
        brv.add_heading(doc, head, 2)
        render_block_paras(doc, paras, size=10.5)
        fname = MAIN_FIG_MAP.get(key)
        brv.add_figure(doc, os.path.join(MAIN_FIG_DIR, fname), width_in=6.5)
        note_par(doc, f"file: figures/main_figures_v2/{fname}")
        if k < len(main_figs) - 1:
            brv.add_page_break(doc)

    # --- Part II supplementary figures --------------------------------------
    brv.add_page_break(doc)
    brv.add_heading(doc, "Part II — Supplementary Figures 1–4 (legends + images)", 1)
    supp_blocks = parse_legend_blocks(SUPP_LEGEND)
    supp_figs = [b for b in supp_blocks if b[0].startswith("Supplementary Figure ")]
    # Table legends live in their own file since the 2026-09-13 decision that the
    # submitted supplementary material is figures-only (tables on request).
    SUPP_TAB_LEGEND = os.path.join(ROOT, "results", "tables", "supplementary",
                                   "supplementary_table_legends_on_request.md")
    supp_tabs = [b for b in parse_legend_blocks(SUPP_TAB_LEGEND)
                 if b[0].startswith("Supplementary Table ")]
    abbr_supp = [b for b in supp_blocks if b[0].startswith("Abbreviations")]
    for k, (head, paras) in enumerate(supp_figs):
        key = head.split(".")[0] + "."
        brv.add_heading(doc, head, 2)
        render_block_paras(doc, paras, size=10.5)
        fname = SUPP_FIG_MAP.get(key)
        brv.add_figure(doc, os.path.join(SUPP_FIG_DIR, fname), width_in=6.3)
        note_par(doc, f"file: figures/supplementary_v2/{fname}")
        if k < len(supp_figs) - 1:
            brv.add_page_break(doc)

    # --- Part III supplementary tables (landscape) --------------------------
    section_to_landscape(doc)
    brv.add_heading(doc, "Part III — Supplementary Tables 1–5 (on-request set; NOT submitted)", 1)
    brv.add_par(doc, ("横向页面排版。按作者 2026-09-13 裁决，投稿的补充材料只含补充图，"
                      "补充表不随稿提交、正文引用已全部删除；本部分保留完整表格与表注，"
                      "供审稿人索取时直接交付。"
                      "2026-09-23 进一步精简：15 张减为 5 张，并按正文首现顺序重排"
                      "（§3.1 全量筛选 → §3.1–3.2 22 命中 MR/复制/共定位 → §3.2 SuSiE 可信集 "
                      "→ §3.3 PheWAS 全基因组显著 → §3.4 pQTL 评估）；"
                      "被删的 10 张因内容已被图完整呈现，改列为无编号源数据（见本部分末尾）。"
                      "Table 1（5,040 行）与 Table 3（SuSiE 756 行）体积过大，"
                      "此处仅给出列结构、前 8 行预览与完整 CSV 路径；其余表全部行完整渲染。"
                      "列数超过 8 的宽表按 8 列横向分块渲染并重复首列（`gene`），"
                      "以避免单元格被挤成单字符换行；每个分块标注其列范围。"
                      "表格内容与 results/tables/supplementary/*.csv 完全一致。"), size=10)
    for k, (head, paras) in enumerate(supp_tabs):
        m = re.match(r"Supplementary Table (\d+)\.", head)
        num = int(m.group(1)) if m else 0
        brv.add_heading(doc, head, 2)
        render_block_paras(doc, paras, size=10)
        fpath = os.path.join(SUPP_TAB_DIR, TABLE_FILES[num])
        nrows = max(0, len(csv_rows(fpath)) - 1)
        if num in BIG_TABLES:
            note_par(doc, (f"Full table file: results/tables/supplementary/{TABLE_FILES[num]} "
                           f"({nrows:,} data rows) — not embedded in full because of size. "
                           f"Preview below shows the column structure and the first "
                           f"{PREVIEW_ROWS} rows."))
            add_csv_table_blocks(doc, fpath, max_rows=PREVIEW_ROWS)
        else:
            note_par(doc, (f"file: results/tables/supplementary/{TABLE_FILES[num]} "
                           f"({nrows:,} data rows) — all rows rendered below."))
            add_csv_table_blocks(doc, fpath)
        if k < len(supp_tabs) - 1:
            brv.add_page_break(doc)

    # --- the one submitted tabular item (numberless) -------------------------
    brv.add_page_break(doc)
    brv.add_heading(doc, "STROBE-MR reporting checklist (submitted separately; no "
                         "Supplementary Table number)", 2)
    strobe_path = os.path.join(SUPP_TAB_DIR, STROBE_FILE)
    n_strobe = max(0, len(csv_rows(strobe_path)) - 1)
    note_par(doc, (f"file: results/tables/supplementary/{STROBE_FILE} ({n_strobe} items) "
                   f"— the only tabular item that accompanies the manuscript as a separate "
                   f"file; it carries no table number because it is a reporting-compliance "
                   f"item rather than data."))
    add_csv_table_blocks(doc, strobe_path)
    brv.add_par(doc, "", size=2)

    # --- numberless source data (retired tables + single-cell summary tables) ---
    brv.add_page_break(doc)
    brv.add_heading(doc, "Source data without a supplementary-table number", 2)
    brv.add_par(doc, ("以下文件不再是编号补充表（内容已在图中呈现），但作为图件源数据保留，"
                      "随「按需索取」材料一并提供。"), size=10)
    for fn, desc in RETIRED_SOURCES:
        note_par(doc, f"{fn} — {desc}")
    brv.add_par(doc, "", size=2)

    # --- Part IV numerical verification notes (portrait) --------------------
    section_to_portrait(doc)
    brv.add_heading(doc, "Part IV — Numerical verification notes (author review, not for journal)", 1)
    found_any = False
    for path, blocks in ((MAIN_LEGEND, main_blocks), (SUPP_LEGEND, supp_blocks)):
        for head, paras in blocks:
            if head.startswith("数值核验注"):
                found_any = True
                src = ("figures/main_figures_v2/figure_legends_main.md"
                       if path == MAIN_LEGEND
                       else "figures/supplementary_v2/figure_legends_supplementary.md")
                brv.add_heading(doc, f"{head} — source: {src}", 3)
                render_block_paras(doc, paras, size=10)
    if not found_any:
        brv.add_par(doc, "(no verification notes found)", italic=True, color=GRAY)

    # --- Part V abbreviations ------------------------------------------------
    brv.add_heading(doc, "Part V — Abbreviations", 1)
    for head, paras in abbr_main:
        brv.add_heading(doc, "Main figure legends — " + head, 3)
        render_block_paras(doc, paras, size=10)
    for head, paras in abbr_supp:
        brv.add_heading(doc, "Supplementary legends — " + head, 3)
        render_block_paras(doc, paras, size=10)

    doc.save(OUT_REVIEW)
    # embedded figure count is derived from the legend sources, not hard-coded, so the
    # QA assertion below cannot go stale when plates are added or removed
    return doc, len(main_figs) + len(supp_figs)


def qa(doc, path, expect_images=None):
    n_imgs = len(doc.inline_shapes)
    n_tabs = len(doc.tables)
    n_pars = len(doc.paragraphs)
    size_kb = os.path.getsize(path) // 1024
    print(f"[QA] {os.path.basename(path)}: paragraphs={n_pars} tables={n_tabs} "
          f"inline_images={n_imgs} size={size_kb} KB")
    if expect_images is not None:
        assert n_imgs == expect_images, f"expected {expect_images} images, got {n_imgs}"
    return n_pars, n_tabs, n_imgs


if __name__ == "__main__":
    d1 = build_print()
    d2, n_fig_imgs = build_review()
    qa(d1, OUT_PRINT, expect_images=0)
    qa(d2, OUT_REVIEW, expect_images=n_fig_imgs)
    print("OK ->", OUT_PRINT)
    print("OK ->", OUT_REVIEW)
