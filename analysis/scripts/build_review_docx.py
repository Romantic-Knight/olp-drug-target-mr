# -*- coding: utf-8 -*-
"""
Build a single printable review .docx for the OLP drug-target MR manuscript.

Contents (in order):
  A. Title page / author block (manuscript_title_page.md)
  B. Main text: Abstract .. References (manuscript_main_text.md, up to '# References')
  C. Main figures 1-7: legend (figure_legends_main.md) + embedded PNG
  D. Supplementary material:
       S-Figs 1-10: legend + embedded PNG
       STROBE-MR checklist rendered as a compact table (30 rows; submitted separately)
  E. Author review appendix (not for journal):
       verification notes (Chinese + English numeric checks), [check] to-do list
       S-Table S1-S15 legends (on request only - by author decision 2026-09-13 the
       supplementary material submitted with the manuscript is figures-only)

Companion data tables (CSV) are deliberately not embedded and not submitted; full paths
are listed so the tables can be sent on request.
"""
import csv, os, re, sys
from docx import Document
from docx.shared import Pt, Inches, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = r"D:\AI-Wisp-Science\Drug Target MR"
REPORTS = os.path.join(ROOT, "results", "reports")
MAIN_FIG = os.path.join(ROOT, "figures", "main_figures")
SUPP_FIG = os.path.join(ROOT, "figures", "supplementary")
SUPP_TAB = os.path.join(ROOT, "results", "tables", "supplementary")
OUT = os.path.join(REPORTS, "OLP_MR_review_compilation.docx")

PAGE_W = Inches(8.27)   # A4 portrait
PAGE_H = Inches(11.69)
MARGIN = Inches(0.79)
CONTENT_W_IN = 8.27 - 2 * 0.79  # ~6.69 inches (float, for width computation)
CONTENT_W = Inches(CONTENT_W_IN)

CHECK_COLOR = RGBColor(0xC0, 0x39, 0x2B)  # red-ish for [check] flags

# ----------------------------------------------------------------------------
# Fonts: body = Times New Roman (msyh reserved for Chinese), fallback handled
# by Word automatically. Set eastasia font for CJK runs.
# ----------------------------------------------------------------------------
BODY = "Times New Roman"
CJK = "Microsoft YaHei"
MONO = "Consolas"

# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def set_run_font(run, name=BODY, size=11, bold=None, italic=None, color=None,
                 eastasia=CJK):
    run.font.name = name
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic
    if color is not None:
        run.font.color.rgb = color
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn('w:rFonts'))
    if rfonts is None:
        rfonts = OxmlElement('w:rFonts')
        rpr.append(rfonts)
    rfonts.set(qn('w:ascii'), name)
    rfonts.set(qn('w:hAnsi'), name)
    rfonts.set(qn('w:eastAsia'), eastasia)


def add_par(doc, text="", size=11, bold=False, italic=False, align=None,
            space_before=0, space_after=6, color=None, style=None,
            font=BODY, keep_with_next=False, line=None):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.space_before = Pt(space_before)
    pf.space_after = Pt(space_after)
    if align is not None:
        pf.alignment = align
    if line is not None:
        pf.line_spacing = line
    if keep_with_next:
        pf.keep_with_next = True
    if style:
        p.style = doc.styles[style]
    if text:
        _add_runs(p, text, size=size, bold=bold, italic=italic, color=color,
                  font=font)
    return p


_INLINE_RE = re.compile(r'(<sup>.*?</sup>|<sub>.*?</sub>|`[^`]*`|\*\*[^*]+\*\*|\*[^*]+\*)')


def _add_runs(p, text, size=11, bold=False, italic=False, color=None, font=BODY):
    """Recursive minimal md inline: <sup>/<sub>, **bold**, *italic*, `code`."""
    pos = 0
    for m in _INLINE_RE.finditer(text):
        if m.start() > pos:
            _emit(p, text[pos:m.start()], size, bold, italic, color, font)
        tok = m.group(0)
        if tok.startswith('<sup>'):
            _emit(p, tok[6:-7], size, bold, italic, color, font, superscript=True)
        elif tok.startswith('<sub>'):
            _emit(p, tok[6:-7], size, bold, italic, color, font, subscript=True)
        elif tok.startswith('`'):
            _emit(p, tok[1:-1], size, bold, italic, color, MONO)
        else:
            inner = tok[2:-2] if tok.startswith('**') else tok[1:-1]
            _add_runs(p, inner, size,
                      bold=True if tok.startswith('**') else bold,
                      italic=True if tok.startswith('*') and not tok.startswith('**') else italic,
                      color=color, font=font)
        pos = m.end()
    if pos < len(text):
        _emit(p, text[pos:], size, bold, italic, color, font)


def _emit(p, text, size, bold, italic, color, font, superscript=False, subscript=False):
    if not text:
        return
    # split out [check] ... flags to colour them
    parts = re.split(r'(\[[Cc]heck[^\]]*\])', text)
    for part in parts:
        if not part:
            continue
        r = p.add_run(part)
        if re.fullmatch(r'\[[Cc]heck[^\]]*\]', part):
            set_run_font(r, font, size, True, italic, CHECK_COLOR)
        else:
            set_run_font(r, font, size, bold, italic, color)
        if superscript:
            r.font.superscript = True
        if subscript:
            r.font.subscript = True


def add_heading(doc, text, level=1):
    """Custom heading using body font (Word's default headings are Calibri blue)."""
    sizes = {1: 16, 2: 13.5, 3: 12, 4: 11}
    size = sizes.get(level, 11)
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.space_before = Pt(14 if level <= 2 else 10)
    pf.space_after = Pt(6)
    pf.keep_with_next = True
    _add_runs(p, text, size=size, bold=True)
    # bottom border for level-1/2
    if level <= 2:
        ppr = p._element.get_or_add_pPr()
        pbdr = OxmlElement('w:pBdr')
        bottom = OxmlElement('w:bottom')
        bottom.set(qn('w:val'), 'single')
        bottom.set(qn('w:sz'), '6')
        bottom.set(qn('w:space'), '2')
        bottom.set(qn('w:color'), '444444')
        pbdr.append(bottom)
        ppr.append(pbdr)
    return p


def add_bullet(doc, text, size=11, level=0):
    p = doc.add_paragraph(style='List Bullet' if level == 0 else 'List Bullet 2')
    pf = p.paragraph_format
    pf.space_after = Pt(2)
    pf.left_indent = Inches(0.25 + 0.25 * level)
    _add_runs(p, text, size=size)
    return p


def add_numbered(doc, text, size=11):
    # literal numbering: keep the "1. " prefix in text, avoid Word's global
    # List Number continuation across sections
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.left_indent = Inches(0.25)
    p.paragraph_format.first_line_indent = Inches(-0.25)
    _add_runs(p, text, size=size)
    return p


def add_hr(doc):
    p = doc.add_paragraph()
    ppr = p._element.get_or_add_pPr()
    pbdr = OxmlElement('w:pBdr')
    bottom = OxmlElement('w:bottom')
    bottom.set(qn('w:val'), 'single')
    bottom.set(qn('w:sz'), '4')
    bottom.set(qn('w:space'), '1')
    bottom.set(qn('w:color'), '999999')
    pbdr.append(bottom)
    ppr.append(pbdr)
    p.paragraph_format.space_after = Pt(4)
    return p


def add_figure(doc, path, width_in=None):
    """Embed image, centered. path must exist."""
    from docx.shared import Inches as _In
    if not os.path.exists(path):
        add_par(doc, f"[image missing: {path}]", italic=True, color=CHECK_COLOR)
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(2)
    r = p.add_run()
    if width_in is None:
        width_in = CONTENT_W_IN
    r.add_picture(path, width=_In(width_in))


def add_page_break(doc):
    from docx.enum.text import WD_BREAK
    p = doc.add_paragraph()
    p.add_run().add_break(WD_BREAK.PAGE)


def add_header_footer(doc, title_text):
    sec = doc.sections[0]
    header = sec.header
    hp = header.paragraphs[0]
    hp.text = ""
    r = hp.add_run(title_text)
    set_run_font(r, BODY, 8.5, italic=True, color=RGBColor(0x66, 0x66, 0x66))
    hp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    # footer page numbers
    footer = sec.footer
    fp = footer.paragraphs[0]
    fp.text = ""
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r1 = fp.add_run("Page ")
    set_run_font(r1, BODY, 9)
    fld = OxmlElement('w:fldSimple')
    fld.set(qn('w:instr'), 'PAGE')
    rEl = OxmlElement('w:r')
    rpr = OxmlElement('w:rPr')
    sz = OxmlElement('w:sz'); sz.set(qn('w:val'), '18')
    rpr.append(sz)
    rEl.append(rpr)
    t = OxmlElement('w:t'); t.text = "1"
    rEl.append(t)
    fld.append(rEl)
    fp._element.append(fld)


# ----------------------------------------------------------------------------
# Table rendering for STROBE checklist
# ----------------------------------------------------------------------------
def add_strobe_table(doc, csv_path):
    with open(csv_path, encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    hdr = rows[0]
    body = rows[1:]
    tbl = doc.add_table(rows=1, cols=len(hdr))
    tbl.style = 'Table Grid'
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    widths = [Inches(1.05), Inches(0.35), Inches(1.9), Inches(0.75), Inches(2.64)]
    for i, h in enumerate(hdr):
        cell = tbl.rows[0].cells[i]
        cell.text = ""
        p = cell.paragraphs[0]
        r = p.add_run(h)
        set_run_font(r, BODY, 8.5, bold=True)
        p.paragraph_format.space_after = Pt(0)
        _shade(cell, "E8E8E8")
    # merge domain column? keep simple: repeat domain, small caps color
    for row in body:
        cells = tbl.add_row().cells
        for i, val in enumerate(row):
            if i >= len(hdr):
                break
            cell = cells[i]
            cell.text = ""
            p = cell.paragraphs[0]
            r = p.add_run(val)
            set_run_font(r, BODY, 8)
            p.paragraph_format.space_after = Pt(0)
    # fix widths per cell
    for row in tbl.rows:
        for i, cell in enumerate(row.cells):
            if i < len(widths):
                cell.width = widths[i]
    return tbl


def _shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:color'), 'auto')
    shd.set(qn('w:fill'), hexcolor)
    tcPr.append(shd)


# ----------------------------------------------------------------------------
# CSV path catalog (shown in supp section)
# ----------------------------------------------------------------------------
CSV_CATALOG = [
    ("Supplementary Table 1", "suppl_table_01_full_mr_screen_5040.csv",
     "Full screen: all 5,040 druggable genes, WIDE/INCLAVO MR + multiplicity columns"),
    ("Supplementary Table 2", "suppl_table_02_wide_hits_mr_replication_coloc.csv",
     "22 WIDE Bonferroni hits: IV annotation + F, dual-cohort MR, replication, coloc PP.H4, screen tier"),
    ("Supplementary Table 3", "suppl_table_03_susie_credible_sets.csv",
     "SuSiE credible sets (756): eqtl_full 284 + eqtl 262 + WIDE 158 + INCLAVO 52"),
    ("Supplementary Table 4", "suppl_table_04_phewas_gw_hits.csv",
     "PheWAS genome-wide hits (121) for the five candidate instruments"),
    ("Supplementary Table 5", "suppl_table_05_pqtl_assessment.csv",
     "Cross-ancestry pQTL assessment, GDF11/SUOX (INTERVAL / Knight ADRC / deCODE)"),
]
STROBE_CSV = "strobe_mr_checklist.csv"
# Numberless source data for the figures and for the retired tables (2026-09-23 trim): the
# single-cell summary tables and four curated tables that lost their "Supplementary Table"
# number because a figure already displays their content. Kept so the pack still inventories
# every file a reviewer may ask for.
SUMMARY_SOURCES = [
    ("results/tables/", "geo_gse211630_erosion_candidate_celltype.csv + _erosion_per_sample.csv",
     "Candidate-gene expression by cell type (Figure 6)"),
    ("results/tables/", "geo_gse211630_composition_{per_sample,shift}.csv",
     "Cell-type composition (Supplementary Figure 9A/B)"),
    ("results/tables/", "geo_gse211630_tnk_module_{per_sample,shift}.csv",
     "T/NK module scores (Supplementary Figure 9C)"),
    ("results/tables/", "geo_gse211630_mr_backmap_{gene,subtype}.csv",
     "Tissue-level reverse mapping and Oaxaca terms (Figure 7)"),
    ("results/tables/", "geo_gse211630_liana_*.csv",
     "LIANA cell-cell communication (Supplementary Figure 10)"),
    ("results/tables/supplementary/_retired/",
     "coloc_abf_susie_44runs.csv, pleiotropy_verdict_5x9.csv, docking_results.csv, "
     "steiger_direction.csv",
     "Former Tables 4, 7, 8, 9 - content displayed in Supplementary Figure 3, Figures 5C/D, "
     "3B and 5B respectively"),
]


# ----------------------------------------------------------------------------
# md section splitter: read a .md file, yield (kind, payload)
#   kind: 'title'/'h1'/'h2'/'h3'/'h4'/'text'/'bullet'/'hr'/'blank'
# We treat a leading '# ' line as the document title.
# ----------------------------------------------------------------------------
def iter_md(path, stop_heading=None):
    with open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")
    in_code = False
    for ln in lines:
        s = ln.rstrip()
        if s.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        if not s.strip():
            yield ("blank", "")
            continue
        if re.match(r'^---+$', s.strip()):
            yield ("hr", "")
            continue
        if s.lstrip().startswith(">"):
            yield ("text", re.sub(r'^\s*>\s?', '', s))
            continue
        m = re.match(r'^(#{1,6})\s+(.*)$', s)
        if m:
            lvl = len(m.group(1))
            text = m.group(2).strip()
            if stop_heading and text.startswith(stop_heading):
                return
            yield (f"h{lvl}", text)
            continue
        # bullets: - or *
        if re.match(r'^\s*[-*]\s+', s):
            yield ("bullet", re.sub(r'^\s*[-*]\s+', '', s))
            continue
        # numbered list "1. " style -> keep literal marker in payload
        if re.match(r'^\s*\d+\.\s+', s):
            yield ("num", s.strip())
            continue
        # image syntax ![](path) => placeholder (handled by caller normally)
        if s.startswith("![") :
            yield ("text", s)
            continue
        yield ("text", s)


def render_md(doc, path, base_heading_level=1, stop_heading=None,
              skip_heading=None, include_hr=True, skip_first_h1=False):
    """Render an md file into doc. skips lines whose heading text starts with
    one of skip_heading. Simple pipe tables are rendered as Word tables."""
    first_h1_seen = False
    table_buf = []

    def flush_table():
        nonlocal table_buf
        if table_buf:
            _render_pipe_table(doc, table_buf)
            table_buf = []

    for kind, payload in iter_md(path, stop_heading=stop_heading):
        if skip_heading and kind.startswith("h") and any(
                payload.startswith(sh) for sh in skip_heading):
            continue
        if kind == "blank":
            flush_table()
            continue
        if kind == "hr":
            flush_table()
            if include_hr:
                add_hr(doc)
            continue
        if kind.startswith("h"):
            flush_table()
            lvl = int(kind[1])
            if skip_first_h1 and not first_h1_seen and lvl == 1:
                first_h1_seen = True
                continue
            add_heading(doc, payload, level=min(lvl, 4))
            continue
        if payload.strip().startswith("|"):
            table_buf.append(payload)
            continue
        flush_table()
        if kind == "bullet":
            add_bullet(doc, payload)
        elif kind == "num":
            add_numbered(doc, payload)
        else:
            add_par(doc, payload, size=11)
    flush_table()


def _render_pipe_table(doc, rows):
    """rows: raw md pipe lines -> Word table (skip --- separator rows)."""
    parsed = []
    for line in rows:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if all(re.fullmatch(r':?-{2,}:?', c or "-") for c in cells):
            continue  # separator
        parsed.append(cells)
    if not parsed:
        return
    ncols = max(len(r) for r in parsed)
    tbl = doc.add_table(rows=0, cols=ncols)
    tbl.style = 'Table Grid'
    tbl.alignment = WD_TABLE_ALIGNMENT.LEFT
    for ri, cells in enumerate(parsed):
        row = tbl.add_row().cells
        for ci in range(ncols):
            val = cells[ci] if ci < len(cells) else ""
            cell = row[ci]
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            _add_runs(p, val, size=10, bold=(ri == 0))
            if ri == 0:
                _shade(cell, "E8E8E8")
    add_par(doc, "", size=4, space_after=2)
    return tbl


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def build():
    doc = Document()
    # page geometry
    sec = doc.sections[0]
    sec.page_width = PAGE_W
    sec.page_height = PAGE_H
    sec.left_margin = sec.right_margin = MARGIN
    sec.top_margin = sec.bottom_margin = Inches(0.79)
    # default style
    st = doc.styles['Normal']
    st.font.name = BODY
    st.font.size = Pt(11)
    st.element.rPr.rFonts.set(qn('w:eastAsia'), CJK)

    add_header_footer(doc, "OLP drug-target MR (FinnGen R13) - review compilation (draft, 2026-09)")

    # ==================================================================
    # PART A - title page block
    # ==================================================================
    add_heading(doc, "From germline risk to lesional tissue: druggable-genome Mendelian randomization and single-cell profiling nominate GDF11 inhibition for oral lichen planus", 1)
    render_md(doc, os.path.join(REPORTS, "manuscript_title_page.md"),
              skip_heading=["# Manuscript Title Page"])
    add_par(doc, "Note: CRediT / Funding / COI / Ethics / target journal remain P0-pending (marked in-text).", size=9, italic=True, color=RGBColor(0x88, 0x88, 0x88))
    add_page_break(doc)

    # ==================================================================
    # PART B - main text through References
    # ==================================================================
    main_md = os.path.join(REPORTS, "manuscript_main_text.md")
    render_md(doc, main_md, stop_heading="References", skip_first_h1=True)
    # references section rendered separately (numbered plain lines already
    # come through as 'num')
    # manual: find References heading and list
    ref_started = False
    ref_lines = []
    for kind, payload in iter_md(main_md):
        if kind == "h2" and payload.startswith("References"):
            ref_started = True
            continue
        if ref_started:
            ref_lines.append((kind, payload))
    if ref_lines:
        add_heading(doc, "References", 2)
        for kind, payload in ref_lines:
            if kind in ("blank",):
                continue
            if kind == "num":
                add_numbered(doc, payload)
            elif kind == "text":
                add_par(doc, payload, size=10.5, italic=True)
            elif kind == "hr":
                continue
            else:
                add_par(doc, payload, size=10.5)
    add_page_break(doc)

    # ==================================================================
    # PART C - main figures 1-7 (legend + figure)
    # ==================================================================
    main_legend_md = os.path.join(ROOT, "figures", "figure_legends_main.md")
    fig_map = {
        "Figure 1.": "fig01_design_sources.png",
        "Figure 2.": "fig02_mr_screens.png",
        "Figure 3.": "fig03_candidates_docking.png",
        "Figure 4.": "fig04_colocalisation.png",
        "Figure 5.": "fig05_pleiotropy_susie.png",
        "Figure 6.": "fig06_target_expression.png",
        "Figure 7.": "fig07_backmap_liana.png",
    }
    # Segment legend md into (heading, [paras]) blocks
    blocks = []
    cur_head = None
    cur_paras = []
    for kind, payload in iter_md(main_legend_md):
        if kind == "h2":
            if cur_head:
                blocks.append((cur_head, cur_paras))
            cur_head = payload
            cur_paras = []
        elif kind == "h3" or kind == "h1" or kind == "hr" or kind == "blank":
            continue
        else:
            cur_paras.append((kind, payload))
    if cur_head:
        blocks.append((cur_head, cur_paras))
    add_heading(doc, "Main Figures 1-7 (with legends)", 1)
    for head, paras in blocks:
        if head.startswith("Abbreviations"):
            add_heading(doc, head, 2)
            for k, pl in paras:
                add_par(doc, pl, size=10)
            continue
        add_heading(doc, head, 2)
        for k, pl in paras:
            add_par(doc, pl, size=10.5)
        # find image
        key = head.split(".")[0] + "."
        fname = fig_map.get(key)
        if fname:
            add_figure(doc, os.path.join(MAIN_FIG, fname))
    add_page_break(doc)

    # ==================================================================
    # PART D - supplementary material
    # ==================================================================
    add_heading(doc, "Supplementary Material", 1)
    add_par(doc, "Scope: 10 supplementary figures (S1-S10). By author decision (2026-09-13) the "
                 "submitted supplementary material contains figures only: supplementary tables "
                 "are not submitted, are available from the corresponding author on request, and "
                 "all table citations have been removed from the main text. By a further author "
                 "decision (2026-09-23) the on-request set was trimmed from 15 tables to 5 and "
                 "renumbered into first-citation order; ten former tables were retired because a "
                 "figure already displays their content. The STROBE-MR checklist is a "
                 "reporting-compliance item rather than data, so it accompanies the submission "
                 "as a separate file and is rendered in full below.",
            size=10, italic=True)

    supp_legend_md = os.path.join(SUPP_FIG, "figure_legends_supplementary.md")
    supp_fig_map = {
        "Supplementary Figure 1.": "supp_fig01_22hits_forest.png",
        "Supplementary Figure 2.": "supp_fig02_fdr_sensitivity.png",
        "Supplementary Figure 3.": "supp_fig03_coloc_full_posteriors.png",
        "Supplementary Figure 4.": "supp_fig04_regional_plots.png",
        "Supplementary Figure 5.": "supp_fig05_dataset_qc.png",
        "Supplementary Figure 6.": "supp_fig06_olp_umap.png",
        "Supplementary Figure 7.": "supp_fig07_marker_dotplot.png",
        "Supplementary Figure 8.": "supp_fig08_candidate_expression.png",
        "Supplementary Figure 9.": "supp_fig09_composition_tnk.png",
        "Supplementary Figure 10.": "supp_fig10_liana_communication.png",
    }
    supp_blocks = []
    cur_head = None
    cur_paras = []
    for kind, payload in iter_md(supp_legend_md):
        if kind == "h2":
            if cur_head:
                supp_blocks.append((cur_head, cur_paras))
            cur_head = payload
            cur_paras = []
        elif kind == "h3":
            # the Chinese numerical verification note heading -> appendix; skip
            continue
        elif kind in ("h1", "hr", "blank"):
            continue
        else:
            cur_paras.append((kind, payload))
    if cur_head:
        supp_blocks.append((cur_head, cur_paras))

    for head, paras in supp_blocks:
        if head.startswith("Abbreviations"):
            add_heading(doc, head, 2)
            for k, pl in paras:
                add_par(doc, pl, size=10)
            continue
        # Supplementary Table heading
        is_fig = head.startswith("Supplementary Figure")
        add_heading(doc, head, 2)
        for k, pl in paras:
            add_par(doc, pl, size=10.5)
        # match csv or figure
        if not is_fig:
            # table: add path line
            num = re.match(r"Supplementary Table (\d+)\.", head)
            if num:
                idx = int(num.group(1))
                fn = CSV_CATALOG[idx-1][1] if idx <= len(CSV_CATALOG) else None
                if fn:
                    add_par(doc, f"Full data: results/tables/supplementary/{fn} "
                                f"(not embedded in this compilation).",
                            size=9, italic=True, color=RGBColor(0x44, 0x44, 0x99))
        else:
            key = head.split(".")[0] + "."
            fname = supp_fig_map.get(key)
            if fname:
                add_figure(doc, os.path.join(SUPP_FIG, fname))

    # ---- STROBE-MR checklist rendered in full (submitted as a separate file) ----
    add_heading(doc, "STROBE-MR reporting checklist (submitted separately; no "
                     "Supplementary Table number)", 2)
    strobe_csv = os.path.join(SUPP_TAB, STROBE_CSV)
    add_strobe_table(doc, strobe_csv)

    # ---- Appendix ----
    add_page_break(doc)
    add_heading(doc, "Appendix (for authors only - not for submission)", 1)

    # A. pending [check] todo list harvested from main text + legends
    add_heading(doc, "A. Open [check] items found in this draft", 2)
    add_par(doc, "Extracted verbatim from the manuscript and legend files. "
                 "Resolve each before submission.", size=10, italic=True)
    found = {}
    for path in [main_md, os.path.join(ROOT, "figures", "figure_legends_main.md"),
                 supp_legend_md]:
        with open(path, encoding="utf-8") as f:
            content = f.read()
        for m in re.finditer(r'\[[Cc]heck[^\]]*\]', content):
            snippet = content[max(0, m.start()-90):m.start()].strip().replace("\n", " ")
            if len(snippet) > 90:
                snippet = snippet[-90:]
            found.setdefault((os.path.relpath(path, ROOT), snippet), m.group(0))
    if found:
        for (rel, ctx), flag in sorted(found.items()):
            add_bullet(doc, f"{flag} — \"…{ctx}…\" (file: {rel})", size=10)
    else:
        add_bullet(doc, "No [check] flags remain.", size=10)

    add_hr(doc)
    add_heading(doc, "B. Numerical verification notes (figure/table legends, author review)", 2)
    add_par(doc, "Chinese verification notes from the two legend files (originally marked "
                 "'not for journal'). Kept here for your print-review.", size=10, italic=True)
    for path, skip in [(os.path.join(ROOT, "figures", "figure_legends_main.md"), []),
                       (supp_legend_md, [])]:
        capture = False
        for kind, payload in iter_md(path):
            if kind == "h3" and payload.startswith("数值核验注"):
                capture = True
                continue
            if capture:
                if kind == "h2":
                    break
                if kind == "text" or kind == "bullet" or kind == "num":
                    add_bullet(doc, payload, size=10) if kind == "bullet" else add_par(doc, payload, size=10)

    # C. data/file inventory
    add_hr(doc)
    add_heading(doc, "C. File inventory (for traceability)", 2)
    inv = [
        ("Manuscript main text", "results/reports/manuscript_main_text.md"),
        ("Title page / author block", "results/reports/manuscript_title_page.md"),
        ("Reference verification report", "results/reports/reference_verification_report.md"),
        ("Main figure legends (final)", "figures/figure_legends_main.md"),
        ("Supplementary legends (final)", "figures/supplementary_v2/figure_legends_supplementary.md"),
        ("Main figures (PNG 300dpi + PDF)", "figures/main_figures_v2/"),
        ("Supplementary figures (PNG + PDF)", "figures/supplementary_v2/"),
        ("Supplementary tables (CSV, 10 in supplementary/ + 5 single-cell in tables/; "
         "on request only - not submitted)", "results/tables/"),
        ("Supplementary table legends (on request)",
         "results/tables/supplementary/supplementary_table_legends_on_request.md"),
    ]
    for name, p in inv:
        add_bullet(doc, f"{name}: {p}", size=10)

    # D. on-request supplementary table legends (NOT part of the submission)
    add_hr(doc)
    add_heading(doc, "D. Supplementary table legends (on request - not part of the submission)", 2)
    add_par(doc, "Retained so the paste-ready legends travel with the tables when a reviewer "
                 "requests them. Source file: "
                 "results/tables/supplementary/supplementary_table_legends_on_request.md",
            size=10, italic=True)
    for kind, payload in iter_md(os.path.join(
            SUPP_TAB, "supplementary_table_legends_on_request.md")):
        if kind == "h1":
            continue
        if kind.startswith("h"):
            add_heading(doc, payload, level=3)
            continue
        if kind in ("hr", "blank"):
            continue
        if kind == "bullet":
            add_bullet(doc, payload, size=10)
        else:
            add_par(doc, payload, size=10)

    add_par(doc, "CSV paths for the on-request set (files are delivered with the tables "
                 "when a reviewer asks for them):", size=10, italic=True)
    for label, fn, desc in CSV_CATALOG:
        add_bullet(doc, f"{label}: results/tables/supplementary/{fn} — {desc}", size=9)
    add_bullet(doc, f"STROBE-MR checklist (submitted separately, no table number): "
                    f"results/tables/supplementary/{STROBE_CSV} — 30-item reporting checklist "
                    f"rendered in full above", size=9)
    add_par(doc, "Source data carrying no supplementary-table number (their content is "
                 "displayed in figures; released with the on-request set):",
            size=10, italic=True)
    for folder, fn, desc in SUMMARY_SOURCES:
        add_bullet(doc, f"{fn} — {desc}", size=9)

    doc.save(OUT)
    print("SAVED:", OUT)


if __name__ == "__main__":
    build()
