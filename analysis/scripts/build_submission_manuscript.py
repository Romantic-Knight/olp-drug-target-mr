"""Build the clean, submission-ready manuscript Word file.

Output: results/reports/OLP_MR_submission_manuscript.docx

Contents (English only, nothing else):
  1. Title, authors, affiliation, correspondence
  2. Abstract + keywords
  3. Main text, Sections 1-5 (Introduction ... Declarations and availability)
  4. References, numbered 1-34 (internal verification annotations removed)
  5. Main figure legends, Figure 1-7
  6. Supplementary figure legends, Supplementary Figure 1-10

Deliberately NOT included, and why:
  * figures themselves            - delivered as OLP_MR_submission_figures_tables.docx
  * supplementary tables           - the manuscript cites no numbered table (author decision:
                                     the on-request set is delivered from the figures & tables
                                     file), so none belong in this file
  * STROBE-MR checklist            - submitted as its own separate file
  * author checklist, [check] flags, numerical-verification notes, Chinese parallel text,
    repository paths               - internal working material

Stripped from the manuscript source because a journal must not see it: the two working HTML
comments, the reference-list verification preamble, the per-entry `[verified ...]` annotations,
the trailing internal file note, and the repository paths inside the Data/Code availability and
author-contribution statements. No scientific statement, number or citation was altered; the four
author-completion fields (funding, competing interests, ethics, contributions) are kept because
they are genuinely pending author input.
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
TITLE_MD = os.path.join(REPORTS, "manuscript_title_page.md")
EN_MAIN = os.path.join(ROOT, "figures", "figure_legends_main.md")
EN_SUPPFIG = os.path.join(ROOT, "figures", "supplementary",
                          "figure_legends_supplementary.md")
OUT = os.path.join(REPORTS, "OLP_MR_submission_manuscript.docx")

HEADER = "OLP drug-target Mendelian randomization (FinnGen R13) — manuscript"
VERIFIED_RE = re.compile(r"\s*\[[Vv]erified[^\]]*\]")
KEY_RE = re.compile(r"^(Figure|Supplementary Figure)\s+(\d+)\.")
STRIP_PARA_LEAD = ("*Draft reference list.", "*End of draft main text.")
# internal repository paths inside the availability / contribution statements
PATH_FIXES = [
    ("**Code availability.** Analysis scripts (R 4.4.2 and Python) are available in "
     "`analysis/scripts/` of the project repository; docking inputs/results in "
     "`analysis/docking/`. A reproducibility package with environment locking and a data "
     "dictionary is in preparation (see `results/reports/manuscript_audit_strobe_mr.md`).",
     "**Code availability.** Analysis scripts (R 4.4.2 and Python) are available from the "
     "corresponding author on reasonable request, as are the docking inputs and results. A "
     "reproducibility package with environment locking and a data dictionary is in preparation."),
    ("**Author contributions.** [To be completed by the authors; author/affiliation block in "
     "`results/reports/manuscript_title_page.md`.]",
     "**Author contributions.** [To be completed by the authors.]"),
    ("Analysis scripts are available in `analysis/scripts/` (see Data and code availability).",
     "Analysis scripts are available from the corresponding author on reasonable request "
     "(see Data and code availability)."),
]


def read_lines(path):
    return open(path, encoding="utf-8").read().split("\n")


def clean_manuscript():
    """Transformed manuscript lines, with every strip/fix asserted."""
    lines = read_lines(MAIN_MD)
    out, dropped, fixed, stripped_refs = [], [], [], 0
    for ln in lines:
        s = ln.rstrip()
        if s.lstrip().startswith("<!--"):
            dropped.append(s[:60])
            continue
        if any(s.strip().startswith(p) for p in STRIP_PARA_LEAD):
            dropped.append(s[:60])
            continue
        for old, new in PATH_FIXES:
            if s.strip() == old:
                s = new
                fixed.append(new[:48])
        if re.match(r"^\d+\.\s+\S", s):          # reference entry
            s2 = VERIFIED_RE.sub("", s)
            if s2 != s:
                stripped_refs += 1
            s = s2
        out.append(s)
    assert len(dropped) == 4, \
        f"expected 2 working comments + reference preamble + trailing note, got {dropped}"
    assert len(fixed) == 2, f"expected 2 path fixes, got {fixed}"
    assert stripped_refs == 34, f"expected 34 reference annotations stripped, got {stripped_refs}"
    text = "\n".join(out)
    assert "<!--" not in text and "[verified" not in text.lower()
    assert "results/reports/" not in text and "analysis/scripts" not in text
    assert "## Abstract" in text and "## References" in text
    return out


def title_page_lines():
    """Title, authors, affiliation and correspondence only; working notes dropped."""
    lines = read_lines(TITLE_MD)
    keep, section, fields = [], None, {"Title": [], "Authors": [], "Correspondence": []}
    for ln in lines:
        s = ln.rstrip()
        if s.startswith("## "):
            section = s[3:].strip()
            continue
        if section in fields:
            fields[section].append(s)
    assert fields["Title"] and fields["Authors"] and fields["Correspondence"]
    title = next(l[2:].strip() for l in fields["Title"] if l.startswith("> "))
    authors = [l for l in fields["Authors"] if l.strip() and not l.startswith("#")]
    corr = [l[2:].strip() if l.startswith("<sup>*</sup> ") else l.strip()
            for l in fields["Correspondence"] if l.strip() and not l.startswith("#")]
    for blk in (authors, corr):
        assert not any(re.search(r"[\u4e00-\u9fff]", x) for x in blk), "Chinese leaked"
    return title, authors, corr


def split_sections(lines):
    """(title_line, abstract_lines, body_lines, reference_lines)."""
    idx = {}
    for i, ln in enumerate(lines):
        if ln.startswith("# "):
            idx["title"] = i
        elif ln.startswith("## Abstract"):
            idx["abstract"] = i
        elif ln.startswith("## 1. Introduction"):
            idx["body"] = i
        elif ln.startswith("## References"):
            idx["refs"] = i
    for k in ("title", "abstract", "body", "refs"):
        assert k in idx, f"missing section {k} in the manuscript"
    return (lines[idx["title"]], lines[idx["abstract"]:idx["body"]],
            lines[idx["body"]:idx["refs"]], lines[idx["refs"]:])


def render_lines(doc, lines, size=11, skip_h1=False):
    for ln in lines:
        s = ln.rstrip()
        if not s.strip() or s.lstrip().startswith("<!--"):
            continue
        if re.match(r"^---+$", s.strip()):
            brv.add_hr(doc)
            continue
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            lvl = len(m.group(1))
            if lvl == 1 and skip_h1:
                continue
            brv.add_heading(doc, m.group(2).strip(), level=min(lvl, 4))
            continue
        if re.match(r"^\s*[-*]\s+", s):
            brv.add_bullet(doc, re.sub(r"^\s*[-*]\s+", "", s), size=size)
            continue
        if re.match(r"^\s*\d+\.\s+", s):
            brv.add_numbered(doc, s.strip(), size=size)
            continue
        if re.match(r"^\s*\|", s):            # pipe table row (not used by the body)
            brv.add_par(doc, s.strip(), size=size)
            continue
        if s.strip().startswith("> "):
            brv.add_par(doc, s.strip()[2:], size=size)
            continue
        if s.strip().startswith("*") and s.strip().endswith("*") and len(s.strip()) > 2:
            brv.add_par(doc, s.strip(), size=size)
            continue
        brv.add_par(doc, s.strip(), size=size)


def legend_blocks(path, kind):
    out = []
    for head, paras in bpr.parse_legend_blocks(path):
        m = KEY_RE.match(head)
        if m and m.group(1) == kind:
            keep = [(k, p) for k, p in paras
                    if not p.strip().startswith(("**Source tables:**", "**File:**",
                                                 "**中文对照", "`results/tables/"))]
            out.append((int(m.group(2)), head, keep))
    return sorted(out)


def build():
    lines = clean_manuscript()
    title_src, authors, corr = title_page_lines()
    title, abstract, body, refs = split_sections(lines)
    assert title_src == title[2:].strip() or title_src in title, \
        "title page title and manuscript H1 have diverged"
    main_legends = legend_blocks(EN_MAIN, "Figure")
    supp_legends = legend_blocks(EN_SUPPFIG, "Supplementary Figure")
    assert len(main_legends) == 7 and len(supp_legends) == 10

    doc = bpr.new_doc(HEADER)
    brv.add_heading(doc, title[2:].strip(), 1)
    for a in authors:
        brv.add_par(doc, a, size=11)
    brv.add_par(doc, "", size=6)
    for c in corr:
        brv.add_par(doc, c, size=10.5)
    brv.add_page_break(doc)

    render_lines(doc, abstract, size=10.5)
    brv.add_hr(doc)
    render_lines(doc, body, size=10.5)
    brv.add_page_break(doc)

    render_lines(doc, refs, size=10)
    brv.add_page_break(doc)

    brv.add_heading(doc, "Figure Legends", 1)
    for num, head, paras in main_legends:
        brv.add_heading(doc, head, 2)
        bpr.render_block_paras(doc, paras, size=10.5)
        if num != main_legends[-1][0]:
            brv.add_page_break(doc)
    brv.add_page_break(doc)

    brv.add_heading(doc, "Supplementary Figure Legends", 1)
    for num, head, paras in supp_legends:
        brv.add_heading(doc, head, 2)
        bpr.render_block_paras(doc, paras, size=10.5)
        if num != supp_legends[-1][0]:
            brv.add_page_break(doc)

    doc.save(OUT)
    return doc, len(main_legends) + len(supp_legends)


def qa(doc, n_legends):
    txt = "\n".join(p.text for p in doc.paragraphs)
    print(f"[QA] {os.path.basename(OUT)}: paragraphs={len(doc.paragraphs)} "
          f"tables={len(doc.tables)} size={os.path.getsize(OUT) // 1024} KB "
          f"(legends={n_legends})")
    assert len(doc.tables) == 0, "this manuscript file should contain no embedded tables"
    for need in ("Abstract", "1. Introduction", "## 2. Methods".replace("## ", ""),
                 "References", "Figure Legends", "Supplementary Figure Legends"):
        assert need in txt, f"missing section: {need}"
    for bad in ("[verified", "<!--", "Source tables", "中文", "数值核验注", "results/",
                "analysis/scripts", "on-request", "not submitted"):
        assert bad not in txt, f"internal content leaked: {bad!r}"
    n_refs = len(re.findall(r"(?m)^\d+\.\s\S", txt))
    print(f"[QA] numbered reference entries: {n_refs}; no internal/Chinese content present")
    assert n_refs == 34, f"expected 34 references, found {n_refs}"


if __name__ == "__main__":
    d, n = build()
    qa(d, n)
