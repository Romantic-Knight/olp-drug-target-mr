"""Build the latest submission-ready manuscript Word file.

Output: results/reports/OLP_MR_manuscript_submission.docx

Content: title page (title / authors / affiliation / correspondence) -> Abstract -> main text
(Sections 1-5, i.e. Introduction, Methods, Results, Discussion, Declarations and availability)
-> References 1-34 in NLM first-citation order. Nothing else: no Chinese, no author checklist, no
internal provenance.

Stripped because these are internal review artefacts, not manuscript content:
  - HTML comments (`<!-- ... -->`, e.g. the title-page pointer and the authors note)
  - the `[verified ...]` annotations on the 34 reference entries
  - the italic draft note above the reference list ("*Draft reference list. The original 18
    entries were verified on ...*")
  - the trailing "*End of draft main text. Supplementary legends: ...*" line

Kept verbatim: all body text, including the four author-completion fields in Section 5
(Funding / Competing interests / Ethics / Author contributions), which are genuine pending author
input and are reported to the user rather than invented here. The figure and supplementary-figure
legends are NOT duplicated in this file: they travel with the figures and tables in
`OLP_MR_submission_figures_tables.docx`.

Also reports, without changing the text, the body's figure cross-references and any internal
repository path that would be meaningless to a reader.
"""
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_review_docx as brv          # noqa: E402
import build_print_review_docx as bpr    # noqa: E402

ROOT = brv.ROOT
REPORTS = os.path.join(ROOT, "results", "reports")
MAIN_MD = bpr.MAIN_MD
OUT = os.path.join(REPORTS, "OLP_MR_manuscript_submission.docx")
HEADER = "OLP drug-target Mendelian randomization (FinnGen R13) — manuscript"

VERIFIED_RE = re.compile(r"\s*\[verified.*?\]", flags=re.S)
INTERNAL_PATHS = re.compile(r"`?(?:analysis/scripts|analysis/docking|results/reports)/[^`\s)]*`?")


def strip_verified(text):
    return re.sub(r"\s{2,}", " ", VERIFIED_RE.sub("", text)).strip()


def sections(text):
    """(abstract paragraphs, body lines, reference entries, dropped artifacts)."""
    body_text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    body_text = re.sub(r"\n{3,}", "\n\n", body_text)
    head, refs = body_text.split("## References", 1)

    abstract = head.split("## Abstract", 1)[1].split("## 1. Introduction", 1)[0]
    abstract = [ln.strip() for ln in abstract.strip().split("\n")
                if ln.strip() and not re.match(r"^-{3,}$", ln.strip())]
    body = "## 1. Introduction" + head.split("## 1. Introduction", 1)[1]

    entries = [strip_verified(ln) for ln in refs.split("\n") if re.match(r"^\d+\. ", ln)]
    dropped = [ln.strip() for ln in refs.split("\n")
               if ln.strip().startswith("*Draft reference list") or ln.strip().startswith("*End of")]
    return abstract, body, entries, dropped


def build():
    raw = open(MAIN_MD, encoding="utf-8").read()
    abstract, body, entries, dropped = sections(raw)
    title = next(ln[2:].strip() for ln in raw.split("\n") if ln.startswith("# "))
    if title != bpr.TITLE:
        raise SystemExit("manuscript H1 and build_print_review_docx.TITLE have diverged")

    doc = bpr.new_doc(HEADER)
    # ---- title page --------------------------------------------------------
    brv.add_heading(doc, title, 1)
    brv.add_par(doc, bpr.AUTHORS_HTML, size=12, align=brv.WD_ALIGN_PARAGRAPH.CENTER)
    brv.add_par(doc, bpr.AFFIL_HTML, size=10.5, align=brv.WD_ALIGN_PARAGRAPH.CENTER)
    brv.add_par(doc, bpr.CORRESP_HTML, size=10, italic=True,
                align=brv.WD_ALIGN_PARAGRAPH.CENTER, color=bpr.GRAY)
    brv.add_hr(doc)

    # ---- abstract ----------------------------------------------------------
    brv.add_heading(doc, "Abstract", 2)
    for para in abstract:
        brv.add_par(doc, para, size=10.5)

    # ---- body (Sections 1-5) ----------------------------------------------
    tmp = os.path.join(REPORTS, "_tmp_manuscript_body.md")
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        f.write(body)
    try:
        brv.render_md(doc, tmp, stop_heading="References")
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)

    # ---- references 1-34, annotations stripped -----------------------------
    brv.add_page_break(doc)
    brv.add_heading(doc, "References", 2)
    for i, entry in enumerate(entries, 1):
        if not entry.startswith(f"{i}. "):
            raise SystemExit(f"reference numbering out of order at entry {i}: {entry[:40]!r}")
        brv.add_numbered(doc, entry, size=10)

    doc.save(OUT)
    return doc, entries, dropped


def report(doc):
    text = "\n".join(p.text for p in doc.paragraphs)
    body = text.split("References", 1)[0]
    figs = sorted({int(n) for n in re.findall(r"(?<!Supplementary )\bFigure (\d+)", body)})
    sfigs = sorted({int(n) for n in re.findall(r"Supplementary Figure (\d+)", body)})
    print(f"[QA] {os.path.basename(OUT)}: paragraphs={len(doc.paragraphs)} "
          f"tables={len(doc.tables)} images={len(doc.inline_shapes)} "
          f"size={os.path.getsize(OUT) // 1024} KB")
    for need in ("Abstract", "1. Introduction", "2. Methods", "3. Results", "4. Discussion",
                 "5. Declarations and availability", "References"):
        assert need in text, f"missing section: {need}"
    for bad in ("[verified", "<!--", "中文", "数值核验", "Author checklist", "End of draft"):
        assert bad not in text, f"internal artefact leaked into the manuscript: {bad!r}"
    print("[QA] internal artefacts absent (no [verified] notes, comments or Chinese)")
    print(f"[QA] main figures cited : {figs}  -> {'OK' if figs == list(range(1, 8)) else 'CHECK'}")
    n_supp = len(glob.glob(os.path.join(ROOT, "figures", "supplementary_v2", "suppfig*_*.png")))
    miss = [n for n in range(1, n_supp + 1) if n not in sfigs]
    print(f"[QA] supp figures cited : {sfigs}")
    if miss:
        print(f"[QA] !! supplementary figures NOT cited in the text: {miss} "
              f"(journals normally require every supplementary item to be cited)")
    pending = re.findall(r"\[To be completed[^\]]*\]", text)
    print(f"[QA] author-completion fields still pending: {len(pending)}")
    paths = sorted(set(INTERNAL_PATHS.findall(text)))
    if paths:
        print(f"[QA] !! internal repository paths still in the text: {paths}")


if __name__ == "__main__":
    d, ents, drops = build()
    print(f"reference entries: {len(ents)} | dropped draft artefacts: {len(drops)}")
    report(d)
