"""Build the bilingual supervisor-review Word deliverable for the OLP drug-target MR study.

Output: results/reports/OLP_MR_supervisor_review_bilingual.docx

Structure
  1. English title (read from the manuscript H1) + Chinese title/authors/scope note
  2. Abstract, English original (read from the manuscript) + Chinese translation
  3. Research analysis workflow (Chinese, step by step) - the part to be reviewed
  4. Part I  - Figures and Tables: main figures 1-7, each as an embedded PNG followed by
               the verbatim English legend and, underneath it, the Chinese legend
  5. Part II - Supplementary figures 1-4 (same bilingual block) and supplementary tables
               S1-S15 (bilingual legend + table body; landscape section)

Bilingual legend sources (English is never re-typed here):
  figures/figure_legends_main.md
  figures/supplementary/figure_legends_supplementary.md
  results/tables/supplementary/supplementary_table_legends_on_request.md
Chinese legends were authored in results/reports/supervisor_review/0{3,4,5}_*.md and are matched
to the English blocks by figure/table number (asserted at build time, so numbering drift fails).

Reuses the rendering helpers of build_review_docx.py / build_print_review_docx.py.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_review_docx as brv          # noqa: E402  (add_heading, add_par, add_figure, ...)
import build_print_review_docx as bpr    # noqa: E402  (table/legend helpers, figure maps)

ROOT = brv.ROOT
REPORTS = os.path.join(ROOT, "results", "reports")
AUTH = os.path.join(REPORTS, "supervisor_review")
MAIN_MD = os.path.join(REPORTS, "manuscript_main_text.md")

FRONT = os.path.join(AUTH, "00_title_note.md")
ABS_ZH = os.path.join(AUTH, "01_abstract_zh.md")
WORKFLOW = os.path.join(AUTH, "02_analysis_workflow_zh.md")
ZH_MAIN = os.path.join(AUTH, "03_legends_zh_main.md")
ZH_SUPPFIG = os.path.join(AUTH, "04_legends_zh_suppfig.md")
ZH_TABLES = os.path.join(AUTH, "05_legends_zh_tables.md")

EN_MAIN = os.path.join(ROOT, "figures", "main_figures_v2", "figure_legends_main.md")
EN_SUPPFIG = os.path.join(ROOT, "figures", "supplementary_v2",
                           "figure_legends_supplementary.md")
EN_TABLES = os.path.join(ROOT, "results", "tables", "supplementary",
                         "supplementary_table_legends_on_request.md")
OUT = os.path.join(REPORTS, "OLP_MR_supervisor_review_bilingual.docx")

MAIN_FIG_DIR = os.path.join(ROOT, "figures", "main_figures_v2")
SUPP_FIG_DIR = os.path.join(ROOT, "figures", "supplementary_v2")
SUPP_TAB_DIR = os.path.join(ROOT, "results", "tables", "supplementary")

# Numbered on-request tables (2026-09-23 set): Tables 1 and 3 are previewed because of size.
# The STROBE-MR checklist is numberless; the ten retired tables are listed as source data.
BIG_TABLES = bpr.BIG_TABLES                    # {1: 5,040 rows, 3: 756 rows}
RETIRED_SOURCES = bpr.RETIRED_SOURCES
STROBE_FILE = bpr.STROBE_FILE
PREVIEW_ROWS = 8
LABEL_EN = 9.0
GRAY = bpr.GRAY
KEY_RE = re.compile(r"^(Figure|Supplementary Figure|Supplementary Table)\s+(\d+)\.")


# --------------------------------------------------------------------------- sources
def manuscript_title(path):
    """First '# ' heading of the manuscript draft."""
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("# "):
            return ln[2:].strip()
    raise SystemExit(f"no H1 title in {path}")


def manuscript_abstract(path):
    """Abstract paragraphs of the manuscript draft (heading -> next h2/hr)."""
    out, inside = [], False
    for raw in open(path, encoding="utf-8"):
        s = raw.strip()
        if s == "## Abstract":
            inside = True
            continue
        if inside and (s.startswith("## ") or s.startswith("---")):
            break
        if inside and s and not s.startswith("<!--"):
            out.append(s)
    if len(out) < 3:
        raise SystemExit(f"abstract extraction returned {len(out)} paragraphs")
    return out


def numbered_blocks(path, kind):
    """[(number, heading, paras)] for one block kind, keyed by the number in the heading."""
    out = []
    for head, paras in bpr.parse_legend_blocks(path):
        m = KEY_RE.match(head)
        if m and m.group(1) == kind:
            out.append((int(m.group(2)), head, paras))
    return out


def by_number(blocks, label):
    d = {}
    for num, head, paras in blocks:
        if num in d:
            raise SystemExit(f"duplicate {label} number {num} in {blocks}")
        d[num] = (head, paras)
    return d


def check_same_numbers(en, zh, label):
    if set(en) != set(zh):
        raise SystemExit(f"{label} numbering mismatch: EN-only {sorted(set(en) - set(zh))}, "
                         f"ZH-only {sorted(set(zh) - set(en))}")
    return sorted(en)


# --------------------------------------------------------------------------- blocks
def bilingual(doc, en_head, en_paras, zh_head, zh_paras, img_path, img_width=6.5,
              size=10.5, label_en="English legend (verbatim from the submission file)",
              label_zh="中文图注（与上方英文图注逐项对应）"):
    brv.add_heading(doc, en_head, 2)
    if img_path:
        brv.add_figure(doc, img_path, width_in=img_width)
        bpr.note_par(doc, f"file: {os.path.relpath(img_path, ROOT).replace(os.sep, '/')}")
    brv.add_par(doc, label_en, size=LABEL_EN, italic=True, color=GRAY)
    bpr.render_block_paras(doc, en_paras, size=size)
    brv.add_par(doc, label_zh, size=LABEL_EN, italic=True, color=GRAY)
    brv.add_par(doc, f"**{zh_head}**", size=size)
    bpr.render_block_paras(doc, zh_paras, size=size)


def add_table_body(doc, num):
    """CSV body (full, or a preview for the two oversized tables)."""
    fname = bpr.TABLE_FILES[num]
    fpath = os.path.join(SUPP_TAB_DIR, fname)
    nrows = max(0, len(bpr.csv_rows(fpath)) - 1)
    if num in BIG_TABLES:
        bpr.note_par(doc, (f"Full table file: results/tables/supplementary/{fname} "
                           f"({nrows:,} data rows) — not embedded in full because of size. "
                           f"Preview below shows the column structure and the first "
                           f"{PREVIEW_ROWS} rows."))
        return bpr.add_csv_table_blocks(doc, fpath, max_rows=PREVIEW_ROWS)
    bpr.note_par(doc, (f"file: results/tables/supplementary/{fname} ({nrows:,} data rows) "
                       f"— all rows rendered below."))
    return bpr.add_csv_table_blocks(doc, fpath)


# --------------------------------------------------------------------------- build
def build():
    title = manuscript_title(MAIN_MD)
    if title != bpr.TITLE:
        raise SystemExit("manuscript H1 and build_print_review_docx.TITLE have diverged")

    en_main = by_number(numbered_blocks(EN_MAIN, "Figure"), "EN main figure")
    en_supp = by_number(numbered_blocks(EN_SUPPFIG, "Supplementary Figure"), "EN supp figure")
    en_tab = by_number(numbered_blocks(EN_TABLES, "Supplementary Table"), "EN table")
    zh_main = by_number(numbered_blocks(ZH_MAIN, "Figure"), "ZH main figure")
    zh_supp = by_number(numbered_blocks(ZH_SUPPFIG, "Supplementary Figure"), "ZH supp figure")
    zh_tab = by_number(numbered_blocks(ZH_TABLES, "Supplementary Table"), "ZH table")
    check_same_numbers(en_main, zh_main, "main figure")
    check_same_numbers(en_supp, zh_supp, "supplementary figure")
    check_same_numbers(en_tab, zh_tab, "supplementary table")
    if len(en_main) != 7 or len(en_supp) != 4 or len(en_tab) != 5:
        raise SystemExit(f"unexpected block counts: {len(en_main)}/{len(en_supp)}/{len(en_tab)}")

    doc = bpr.new_doc("OLP drug-target MR (FinnGen R13) — supervisor review, bilingual "
                      "figures & tables")

    # --- title + scope note -------------------------------------------------
    brv.add_heading(doc, title, 1)
    brv.render_md(doc, FRONT)

    brv.add_heading(doc, "Contents", 2)
    for t in [
        "Abstract（英文原文 + 中文摘要）",
        "研究分析思路（Analysis workflow，分 8 步，含统计口径与待审点）",
        "Part I — Figures and Tables：主图 Figure 1–7（图像 + 英文图注 + 中文图注）",
        "Part II — Supplementary Figures and Tables：补充图 Supplementary Figure 1–10、"
        "补充表 Supplementary Table 1–5（双语表注 + 表体；横向页）+ 无编号的 STROBE-MR 清单",
    ]:
        brv.add_bullet(doc, t, size=10.5)
    bpr.note_par(doc, "图与表均为「上方英文 Legend（与投稿版逐字一致）+ 下方中文对照」；"
                      "英文图注/表注直接读取三个定稿文件，未作改写。")
    brv.add_page_break(doc)

    # --- abstract -----------------------------------------------------------
    brv.add_heading(doc, "Abstract (English — as in the manuscript)", 2)
    for para in manuscript_abstract(MAIN_MD):
        brv.add_par(doc, para, size=10.5)
    brv.add_par(doc, "", size=4)
    brv.render_md(doc, ABS_ZH)
    brv.add_page_break(doc)

    # --- analysis workflow --------------------------------------------------
    brv.render_md(doc, WORKFLOW)
    brv.add_page_break(doc)

    # --- Part I: main figures ----------------------------------------------
    brv.add_heading(doc, "Part I — Figures and Tables", 1)
    brv.add_par(doc, ("本部分给出全部**主图 Figure 1–7**，每张图先嵌入图像，再给出英文图注"
                      "（与投稿版逐字一致）与中文图注。"), size=10.5)
    brv.add_par(doc, ("**关于表格**：按作者 2026-09-13 裁决，正文不含编号表格（原 31 处 "
                      "`Supplementary Table N` 引用已全部删除）；投稿的补充材料只含补充图。"
                      "按作者 2026-09-23 裁决，「按需索取」补充表由 15 张精简为 **5 张**，"
                      "并按正文首现顺序重编号，统一放在 **Part II**；被删的 10 张因内容已被图"
                      "完整呈现，改列为无编号源数据。唯一随稿的表格项是 STROBE-MR 报告清单"
                      "（无编号、单独文件提交）。"
                      "图上每一个数字的底表对应关系为："), size=10.5)
    for t in [
        "Figure 1 ← 可成药基因与 eQTLGen 交集表（druggable_x_eqtlgen.csv 等；无编号表）",
        "Figure 2 ← **Table 1**（5,040 行全量 MR 结果）、**Table 2**（22 个命中）；"
        "敏感性分析见 Supplementary Figure 2",
        "Figure 3 ← **Table 2**（候选证据；对接见 `docking_results.csv`，无编号）",
        "Figure 4 ← **Table 3**（SuSiE 可信集；44 次共定位后验见 Supplementary Figure 3）",
        "Figure 5 ← **Table 4**（121 个全基因组显著关联）、**Table 5**（pQTL 评估，无图）；"
        "多效性判定与 Steiger 见 Figure 5C/D、5B",
        "Figure 6 ← `geo_gse211630_erosion_*`（95 个基因×亚型行 + 逐样本表；无编号）",
        "Figure 7 ← `geo_gse211630_mr_backmap_*`（逐基因分解 + 逐亚型 Oaxaca 项；无编号）",
        "Supplementary Figure 9 ← `geo_gse211630_composition_*`、`geo_gse211630_tnk_module_*`；"
        "Supplementary Figure 10 ← `geo_gse211630_liana_*`（均无编号）",
    ]:
        brv.add_bullet(doc, t, size=10)

    for k, num in enumerate(sorted(en_main)):
        brv.add_page_break(doc)
        en_head, en_paras = en_main[num]
        zh_head, zh_paras = zh_main[num]
        fname = bpr.MAIN_FIG_MAP[f"Figure {num}."]
        bilingual(doc, en_head, en_paras, zh_head, zh_paras,
                  os.path.join(MAIN_FIG_DIR, fname))

    # --- Part II: supplementary figures ------------------------------------
    brv.add_page_break(doc)
    brv.add_heading(doc, "Part II — Supplementary Figures and Tables", 1)
    brv.add_par(doc, ("**补充图 Supplementary Figure 1–10**（顺序与正文行文一致；单细胞图集中在 "
                      "S5–S10）与**补充表 Supplementary Table 1–5**（不随稿、按需索取；"
                      "按正文首现顺序编号：1 全量筛选 → 2 22 命中 MR/复制/共定位 → 3 SuSiE 可信集 "
                      "→ 4 PheWAS 全基因组显著 → 5 pQTL 评估）。Table 1（5,040 行）与 Table 3"
                      "（756 行）体积过大，仅给列结构与前 8 行预览；末尾另附无编号的 STROBE-MR "
                      "清单与无编号源数据清单。"), size=10)

    for k, num in enumerate(sorted(en_supp)):
        brv.add_page_break(doc)
        en_head, en_paras = en_supp[num]
        zh_head, zh_paras = zh_supp[num]
        fname = bpr.SUPP_FIG_MAP[f"Supplementary Figure {num}."]
        bilingual(doc, en_head, en_paras, zh_head, zh_paras,
                  os.path.join(SUPP_FIG_DIR, fname), img_width=6.4)

    # --- Part II continued: supplementary tables (landscape) ---------------
    bpr.section_to_landscape(doc)
    brv.add_heading(doc, "Supplementary Tables 1–5 (on-request set; not submitted)", 1)
    brv.add_par(doc, ("横向页面。表注同样采用「英文在上、中文在下」。表格内容与 "
                      "results/tables/supplementary/ 下的 CSV 完全一致。"), size=10)

    n_csv = 0
    for k, num in enumerate(sorted(en_tab)):
        if k:
            brv.add_page_break(doc)
        en_head, en_paras = en_tab[num]
        zh_head, zh_paras = zh_tab[num]
        brv.add_heading(doc, en_head, 2)
        brv.add_par(doc, "English legend (verbatim from the submission file)",
                    size=LABEL_EN, italic=True, color=GRAY)
        bpr.render_block_paras(doc, en_paras, size=10)
        brv.add_par(doc, "中文表注（与上方英文表注逐项对应）", size=LABEL_EN, italic=True,
                    color=GRAY)
        brv.add_par(doc, f"**{zh_head}**", size=10)
        bpr.render_block_paras(doc, zh_paras, size=10)
        brv.add_par(doc, "", size=4)
        n_csv += add_table_body(doc, num)

    # --- the one submitted tabular item (numberless) ------------------------
    brv.add_page_break(doc)
    brv.add_heading(doc, "STROBE-MR reporting checklist (submitted separately; no table number)",
                    2)
    strobe_path = os.path.join(SUPP_TAB_DIR, STROBE_FILE)
    n_strobe = max(0, len(bpr.csv_rows(strobe_path)) - 1)
    bpr.note_par(doc, (f"file: results/tables/supplementary/{STROBE_FILE} ({n_strobe} items) "
                       f"— 申报合规文件、非数据；因是唯一随稿的表格项，稿件中不占编号。"))
    n_csv += bpr.add_csv_table_blocks(doc, strobe_path)

    # --- numberless source data (retired tables + single-cell summary tables) ---
    brv.add_page_break(doc)
    brv.add_heading(doc, "Source data without a supplementary-table number", 2)
    brv.add_par(doc, ("以下 10 张原补充表不再是编号表（内容已在图中完整呈现），作为图件源数据保留，"
                      "随「按需索取」材料一并提供。"), size=10)
    for fn, desc in RETIRED_SOURCES:
        bpr.note_par(doc, f"{fn} — {desc}")

    doc.save(OUT)
    return doc, n_csv, len(en_main) + len(en_supp)


def qa(doc, n_csv, n_figs):
    n_imgs = len(doc.inline_shapes)
    n_tabs = len(doc.tables)
    n_pars = len(doc.paragraphs)
    size_kb = os.path.getsize(OUT) // 1024
    print(f"[QA] {os.path.basename(OUT)}: paragraphs={n_pars} tables={n_tabs} "
          f"inline_images={n_imgs} size={size_kb} KB (csv tables rendered={n_csv})")
    assert n_imgs == n_figs, f"expected {n_figs} embedded figures, got {n_imgs}"
    # +1 = the pipe table inside 02_analysis_workflow_zh.md
    assert n_tabs == n_csv + 1, f"expected {n_csv + 1} tables, got {n_tabs}"
    for need in ("研究分析思路", "中文图注", "中文表注", "Supplementary Table 5",
                 "STROBE-MR reporting checklist", "Source data without a supplementary-table "
                 "number"):
        assert need in "\n".join(p.text for p in doc.paragraphs), f"missing section: {need}"


if __name__ == "__main__":
    d, n_csv, n_figs = build()
    qa(d, n_csv, n_figs)
