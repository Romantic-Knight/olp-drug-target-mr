"""Build the clean, submission-oriented figures & tables Word file.

Output: results/reports/OLP_MR_submission_figures_tables.docx

Contents (nothing else - no Chinese, no author notes, no file paths, no QA material):
  1. manuscript title
  2. Main Figures 1-7          - embedded PNG + the official English legend
  3. Supplementary Figures 1-4 - embedded PNG + the official English legend
  4. Supplementary Tables 1-5   - official English legend + embedded table body
  5. STROBE-MR reporting checklist (the one tabular item submitted, numberless)

Stripped from the legend sources because they are internal provenance rather than legend text:
  `**Source tables:** ...`, `**File:** ...`, `**中文对照：** ...`, the `results/tables/...` path
  lines in the table legends, the `Abbreviations` blocks and the author-only numerical
  verification notes.

Tables 1 (5,040 rows) and 3 (756 rows) exceed EMBED_MAX_ROWS; a Word file containing them in full
would run to hundreds of pages, so those two are embedded by reference with a one-line English
note (standard practice: oversized tables travel as separate files). Everything else is embedded
in full. Wide tables are split into readable column blocks by the shared renderer.

Reuses the rendering helpers of build_review_docx.py / build_print_review_docx.py.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_review_docx as brv          # noqa: E402
import build_print_review_docx as bpr    # noqa: E402

ROOT = brv.ROOT
REPORTS = os.path.join(ROOT, "results", "reports")
MAIN_MD = os.path.join(REPORTS, "manuscript_main_text.md")
EN_MAIN = os.path.join(ROOT, "figures", "main_figures_v2", "figure_legends_main.md")
EN_SUPPFIG = os.path.join(ROOT, "figures", "supplementary_v2",
                           "figure_legends_supplementary.md")
EN_TABLES = os.path.join(ROOT, "results", "tables", "supplementary",
                         "supplementary_table_legends_on_request.md")
MAIN_FIG_DIR = os.path.join(ROOT, "figures", "main_figures_v2")
SUPP_FIG_DIR = os.path.join(ROOT, "figures", "supplementary_v2")
SUPP_TAB_DIR = os.path.join(ROOT, "results", "tables", "supplementary")
OUT = os.path.join(REPORTS, "OLP_MR_submission_figures_tables.docx")

EMBED_MAX_ROWS = 300        # above this a table is embedded by reference, not in full
HEADER = "OLP drug-target Mendelian randomization (FinnGen R13) — figures and tables"
KEY_RE = re.compile(r"^(Figure|Supplementary Figure|Supplementary Table)\s+(\d+)\.")
# legend lines that are repo provenance (or the Chinese parallel text), not legend text
STRIP_PREFIXES = ("**Source tables:**", "**File:**", "**中文对照", "`results/tables/")


def manuscript_title(path):
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("# "):
            return ln[2:].strip()
    raise SystemExit(f"no H1 title in {path}")


def official_paras(paras):
    """Drop the internal provenance / Chinese lines so only legend text survives."""
    return [(k, p) for k, p in paras if not p.strip().startswith(STRIP_PREFIXES)]


def numbered_blocks(path, kind):
    out = []
    for head, paras in bpr.parse_legend_blocks(path):
        m = KEY_RE.match(head)
        if m and m.group(1) == kind:
            out.append((int(m.group(2)), head, official_paras(paras)))
    return sorted(out)


def figure_section(doc, blocks, fig_dir, fig_map, width_in):
    for num, head, paras in blocks:
        brv.add_heading(doc, head, 2)
        brv.add_figure(doc, os.path.join(fig_dir, fig_map[f"{head.split('.')[0]}."]),
                       width_in=width_in)
        bpr.render_block_paras(doc, paras, size=10.5)
        if num != blocks[-1][0]:
            brv.add_page_break(doc)


def table_section(doc, blocks):
    """Legend + embedded body (or a reference note for oversized tables). Returns #tables."""
    n = 0
    for num, head, paras in blocks:
        brv.add_heading(doc, head, 2)
        bpr.render_block_paras(doc, paras, size=10)
        fname = bpr.TABLE_FILES[num]
        fpath = os.path.join(SUPP_TAB_DIR, fname)
        nrows = max(0, len(bpr.csv_rows(fpath)) - 1)
        if nrows > EMBED_MAX_ROWS:
            brv.add_par(doc, f"The complete table ({nrows:,} data rows) is provided as the "
                             f"separate file {fname}.", size=10, italic=True)
        else:
            n += bpr.add_csv_table_blocks(doc, fpath)
        if num != blocks[-1][0]:
            brv.add_page_break(doc)
    return n


def build():
    title = manuscript_title(MAIN_MD)
    if title != bpr.TITLE:
        raise SystemExit("manuscript H1 and build_print_review_docx.TITLE have diverged")
    main_figs = numbered_blocks(EN_MAIN, "Figure")
    supp_figs = numbered_blocks(EN_SUPPFIG, "Supplementary Figure")
    tables = numbered_blocks(EN_TABLES, "Supplementary Table")
    if len(main_figs) != 7 or len(supp_figs) != 4 or len(tables) != 5:
        raise SystemExit(f"unexpected blocks: {len(main_figs)}/{len(supp_figs)}/{len(tables)}")
    if [n for n, _, _ in main_figs] != list(range(1, 8)) \
            or [n for n, _, _ in supp_figs] != list(range(1, 5)) \
            or [n for n, _, _ in tables] != list(range(1, 6)):
        raise SystemExit("figure/table numbering is not contiguous")

    doc = bpr.new_doc(HEADER)
    brv.add_heading(doc, title, 1)

    brv.add_heading(doc, "Main Figures", 1)
    figure_section(doc, main_figs, MAIN_FIG_DIR, bpr.MAIN_FIG_MAP, 6.5)
    brv.add_page_break(doc)

    brv.add_heading(doc, "Supplementary Figures", 1)
    figure_section(doc, supp_figs, SUPP_FIG_DIR, bpr.SUPP_FIG_MAP, 6.4)

    bpr.section_to_landscape(doc)
    brv.add_heading(doc, "Supplementary Tables", 1)
    n_tab = table_section(doc, tables)

    brv.add_page_break(doc)
    brv.add_heading(doc, "STROBE-MR reporting checklist", 2)
    strobe = os.path.join(SUPP_TAB_DIR, bpr.STROBE_FILE)
    n_tab += bpr.add_csv_table_blocks(doc, strobe)

    doc.save(OUT)
    return doc, n_tab, len(main_figs) + len(supp_figs)


def qa(doc, n_tab, n_figs):
    n_imgs = len(doc.inline_shapes)
    n_tabs = len(doc.tables)
    txt = "\n".join(p.text for p in doc.paragraphs)
    print(f"[QA] {os.path.basename(OUT)}: paragraphs={len(doc.paragraphs)} tables={n_tabs} "
          f"inline_images={n_imgs} size={os.path.getsize(OUT) // 1024} KB (csv tables={n_tab})")
    assert n_imgs == n_figs, f"expected {n_figs} embedded figures, got {n_imgs}"
    assert n_tabs == n_tab, f"expected {n_tab} tables, got {n_tabs}"
    for need in ("Main Figures", "Supplementary Figures", "Supplementary Tables",
                 "STROBE-MR reporting checklist"):
        assert need in txt, f"missing section: {need}"
    # nothing internal may survive
    for bad in ("Source tables", "中文", "数值核验注", "results/tables/", "Abbreviations",
                "English legend", "on-request", "not submitted", "_retired"):
        assert bad not in txt, f"internal content leaked into the submission file: {bad!r}"
    print("[QA] no internal/Chinese/provenance content present")


if __name__ == "__main__":
    d, n_tab, n_figs = build()
    qa(d, n_tab, n_figs)
