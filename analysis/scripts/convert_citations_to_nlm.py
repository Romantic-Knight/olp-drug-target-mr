#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert the OLP drug-target MR manuscript from author-year citations to NLM/ICMJE
numbered citations in first-citation order, and re-order the reference list to match.

Style produced (user decision, 2026-09-13):
  * in-text markers are bracketed numerals: ``[1]``, ``[3,4]``, ``[11-15]``
  * the reference list is ordered by first citation in the text and renumbered 1..N
    (strict NLM/ICMJE; an alphabetical list with numeral markers is *not* NLM)
  * multiple references in one bracket are ascending; >= 3 consecutive numbers collapse
    to a hyphen range
  * ``[verified ...]`` provenance annotations inside the reference list are preserved
  * bracketed citations inside a parenthesis that carries other content keep the
    surrounding text: ``(n = 3,301; Sun et al. 2018)`` -> ``(n = 3,301 [26])``
  * a parenthesis that carried nothing but citations is replaced by the bracket alone,
    because the bracket already delimits the citation: ``(Cao et al. 2023)`` -> ``[11]``
    (never ``([11])``)
  * no values are invented: entry metadata is read from the manuscript, and the
    transform refuses to run unless every in-text form maps 1:1 onto a list entry

Only two hand-written exceptions exist (both listed in SPECIAL_GROUPS below and
asserted, so a silent miss is impossible):
  * the FinnGen group ``(Kurki et al. 2023; https://.../access_results)`` 鈥?the citation
    is lifted in front of the URL parenthesis so the URL keeps its own brackets
  * the reference-list preamble loses its now-brittle ``entry N`` numbers; the numbers
    are replaced by the first-citation-order statement

Idempotent guard: refuses to run when no author-year citation forms are found.
Call ``--dry-run`` to print the mapping and the plan without writing anything.

Downstream files that must be re-ordered by hand after this script runs
(they carry their own numbering): ``analysis/scripts/_fetch_ref_pmids.py`` REFS,
``analysis/scripts/build_reference_pmid_list.py`` TAGS.

Usage:
  python analysis/scripts/convert_citations_to_nlm.py --dry-run
  python analysis/scripts/convert_citations_to_nlm.py
"""
import argparse
import ast
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPORTS = os.path.join(ROOT, "results", "reports")
MAIN_MD = os.path.join(REPORTS, "manuscript_main_text.md")
MAP_JSON = os.path.join(ROOT, "results", "tables", "_citation_number_map.json")
# downstream files that carry their own copy of the reference numbering
FETCH_PY = os.path.join(ROOT, "analysis", "scripts", "_fetch_ref_pmids.py")
PMIDLIST_PY = os.path.join(ROOT, "analysis", "scripts", "build_reference_pmid_list.py")
# ... reordered from the archived pre-NLM copies, so re-running cannot double-map the numbers
ARCHIVE = os.path.join(REPORTS, "_pre_renumber_archive")
FETCH_PY_PRE = os.path.join(ARCHIVE, "_fetch_ref_pmids_pre_nlm_20260913.py")
PMIDLIST_PY_PRE = os.path.join(ARCHIVE, "build_reference_pmid_list_pre_nlm_20260913.py")

# Files rewritten with the same in-text rule. The first carries the reference list;
# the others are legend/manifest copies of the same figure caption text.
TEXT_FILES = [
    MAIN_MD,
    os.path.join(ROOT, "figures", "figure_legends_main.md"),
    os.path.join(ROOT, "figures", "figure_manifest.md"),
]

# In-text citation form: "Surname et al. 2005", "Surname & Surname 2008", "Surname 2021",
# plus the numbered consortium name. The negative lookbehind rejects dates such as
# "26 Aug 2026" / "1 Sep 2026" in the availability statements.
CITE_RE = re.compile(
    r"(?<!\d )(?:"
    r"1000 Genomes Project Consortium 2015"
    r"|[A-Z][\w\u00C0-\u024F.\-]*(?: [A-Z][\w\u00C0-\u024F.\-]*)?"
    r"(?:\s+et al\.|\s+(?:and|&)\s+[A-Z][\w\u00C0-\u024F.\-]*)?,?\s*(?:19|20)\d\d"
    r")"
)
# month names / captions that superficially look like "Surname Year" (access dates, figure refs)
STOP_TOKEN = {"Aug", "Sep", "Sept", "September", "Figure", "Supplementary", "Table", "Methods"}

# group text (exactly as it appears, without the outer parentheses) -> full replacement
SPECIAL_GROUPS = {
    "Kurki et al. 2023; https://www.finngen.fi/en/access_results":
        "{n} (https://www.finngen.fi/en/access_results)",
}

PREAMBLE_MARK = "*Draft reference list."
PREAMBLE_ORDER_NOTE = (
    " The list is ordered by first citation in the text and numbered in that order "
    "(NLM/ICMJE numbered style), converted from author-year citations on 2026-09-13 by "
    "`analysis/scripts/convert_citations_to_nlm.py`; in-text markers are bracketed "
    "numerals, e.g. [1], [3,4], [11-15]. The numbering below is therefore not "
    "alphabetical and changes whenever a citation is inserted."
)


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def write(path, text):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def token_surname(tok):
    """'Gonz谩lez-Moles et al. 2021' -> 'Gonz谩lez-Moles'; 'Davey Smith & Ebrahim 2003' -> 'Davey Smith'."""
    body = CITE_RE.match(tok).group(0)
    body = re.sub(r"\s+(?:et al\.|and|&)\s*.*$", "", body)
    body = re.sub(r",?\s*(?:19|20)\d\d$", "", body)
    return body.strip()


def token_year(tok):
    return re.search(r"(?:19|20)\d\d$", CITE_RE.match(tok).group(0)).group(0)


def match_entry(surname, year, entries):
    """Locate the reference-list entry that starts with this author surname, and check its year."""
    hits = [(n, t) for _, n, t in entries
            if re.match(re.escape(surname) + r"(?:[,;]|\s+(?:[A-Z]\.|[A-Z]{1,3}\b))", t)]
    if len(hits) != 1:
        return None
    n, text = hits[0]
    return (n, text) if year in text else None


def parse_ref_section(md):
    """Split the manuscript into (body, lines, preamble index, [(line index, ref number, entry text)], tail index)."""
    idx = md.index("## References")
    body = md[:idx]
    rest = md[idx:]
    lines = rest.split("\n")
    pre_i = next(k for k, l in enumerate(lines) if l.startswith(PREAMBLE_MARK))
    tail_i = next((k for k, l in enumerate(lines) if l.strip() == "---" and k > pre_i), len(lines))
    entries = []
    for k in range(pre_i, tail_i):
        m = re.match(r"^(\d+)\. (.+)$", lines[k])
        if m:
            entries.append((k, int(m.group(1)), m.group(2)))
    return body, lines, pre_i, entries, tail_i


def is_citation(tok):
    """Reject access dates, figure cross-references and GWAS accessions that look like 'Surname Year'."""
    s = token_surname(tok)
    if s in STOP_TOKEN:
        return False
    return not re.search(r"\d", s) or s == "1000 Genomes Project Consortium"


def scan(body):
    """Find every citation in document order.

    Returns (occ, spans, groups) with
      occ    = [(start, end, token, parenthetical group index or None)]
      spans  = [(inner start, inner end)] per parenthetical group
      groups = {group index: [occ indices]}
    """
    spans = [(m.start(1), m.end(1)) for m in re.finditer(r"\(([^()]*)\)", body)]
    occ, groups = [], {}
    for m in CITE_RE.finditer(body):
        tok = m.group(0)
        if not is_citation(tok):
            continue
        gi = next((i for i, (a, b) in enumerate(spans) if a <= m.start() < b), None)
        occ.append((m.start(), m.end(), tok, gi))
        if gi is not None:
            groups.setdefault(gi, []).append(len(occ) - 1)
    return occ, spans, groups


def order_of(occ):
    """First-citation order as a list of author surnames."""
    order = []
    for _, _, tok, _ in occ:
        s = token_surname(tok)
        if s not in order:
            order.append(s)
    return order


def fmt(nums):
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 2:
            out.append(f"{nums[i]}-{nums[j]}")
        else:
            out.extend(str(n) for n in nums[i:j + 1])
        i = j + 1
    return "[" + ",".join(out) + "]"


def rewrite_group(inner, toks, num_of, style):
    """Rewrite one parenthesis (inner text, without the outer parentheses).

    toks are (relative start, relative end, token) inside ``inner``.
    Returns None when the parenthesis carried nothing but citations: the caller then
    drops the parentheses and keeps the bracket alone.
    """
    left = inner.replace("; ", ";")
    for t in toks:
        left = left.replace(t[2], "")
    if not left.replace(";", "").strip():
        return None
    out, last = [], 0
    for s, e, tok in toks:
        out.append(inner[last:s])
        out.append(marker([num_of[token_surname(tok)]], style))
        last = e
    out.append(inner[last:])
    txt = "".join(out)
    txt = re.sub(r";\s*(\[[\d,\-]+\]|\d+)", r" \1", txt)    # drop separator before a marker
    return re.sub(r"\s+\)$", ")", txt)


SUP = str.maketrans("0123456789,-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079,\u207b")


def marker(nums, style):
    m = fmt(nums)
    if style == "superscript":
        return m[1:-1].translate(SUP)
    return m


def narrative(tok, n, style):
    """'Eldjarn et al. 2023' -> 'Eldjarn et al. [27]' (NLM narrative citation keeps the author)."""
    author = re.sub(r",?\s*(?:19|20)\d\d$", "", CITE_RE.match(tok).group(0)).strip()
    return f"{author} {marker([n], style)}"


def reorder_downstream(old_to_new, surname_to_number, dry_run):
    """Renumber/ reorder the two downstream files that carry their own reference numbering.

    ``_fetch_ref_pmids.py`` REFS and ``build_reference_pmid_list.py`` TAGS are keyed by
    reference number, so an unreordered copy would attach every PMID and every short tag to
    the wrong entry. Only the keys and the order change: labels, titles, DOIs, journals,
    years and recorded PMIDs are moved verbatim and asserted to survive.
    """
    new_to_surname = {n: s for s, n in surname_to_number.items()}
    assert sorted(old_to_new) == list(range(1, 35)), sorted(old_to_new)
    assert sorted(old_to_new.values()) == list(range(1, 35))
    written = []

    # --- REFS in _fetch_ref_pmids.py: (n, label, title, doi, journal, year, known_pmid) ---
    src = read(FETCH_PY_PRE)
    m = re.search(r"^REFS = \[\n(.*?)^\]$", src, re.S | re.M)
    assert m, f"REFS literal not found in {os.path.relpath(FETCH_PY_PRE, ROOT)}"
    refs = ast.literal_eval("[\n" + m.group(1) + "\n]")
    assert len(refs) == 34 and all(len(r) == 7 for r in refs)
    assert [r[0] for r in refs] == list(range(1, 35)), "source REFS is not in old-number order"
    rows = {}
    for row in refs:
        n = row[0]
        new = old_to_new[n]
        assert new not in rows, f"duplicate new number {new}"
        rows[new] = (new,) + row[1:]          # the tuple's own number must become the new one
    for n, row in rows.items():
        assert row[0] == n, f"REFS row {n} still carries number {row[0]}"
    new_block = "REFS = [\n" + "".join(
        "    (" + ", ".join(json.dumps(v, ensure_ascii=False) if isinstance(v, str) else str(v)
                           for v in rows[n]) + "),\n" for n in sorted(rows)
    ) + "]\n"
    assert sorted(r[2] for r in rows.values()) == sorted(r[2] for r in refs), "REFS titles changed"
    assert sorted(r[6] for r in rows.values()) == sorted(r[6] for r in refs), "REFS PMIDs changed"
    src2 = src[:m.start()] + new_block + src[m.end() + 1:]
    written.append((FETCH_PY, src, src2))

    # --- TAGS in build_reference_pmid_list.py: {n: "Surname A, et al. Journal. Year"} ---
    src = read(PMIDLIST_PY_PRE)
    m = re.search(r"^TAGS = \{\n(.*?)^\}$", src, re.S | re.M)
    assert m, f"TAGS literal not found in {os.path.relpath(PMIDLIST_PY_PRE, ROOT)}"
    tags = ast.literal_eval("{" + m.group(1) + "}")
    assert len(tags) == 34 and sorted(tags) == list(range(1, 35))
    new_tags = {old_to_new[n]: t for n, t in tags.items()}
    for n, tag in sorted(new_tags.items()):
        # the tag at a new number must name the surname the manuscript cites at that number
        if not tag.startswith(new_to_surname[n]):
            sys.exit(f"TAGS[{n}] = {tag!r} does not start with {new_to_surname[n]!r}")
    new_block = "TAGS = {\n" + "".join(
        f"    {n}: {json.dumps(new_tags[n], ensure_ascii=False)},\n" for n in sorted(new_tags)
    ) + "}\n"
    src2 = src[:m.start()] + new_block + src[m.end() + 1:]
    written.append((PMIDLIST_PY, src, src2))

    for path, old, new in written:
        print(f"{os.path.relpath(path, ROOT)}: {len(old)} -> {len(new)} chars, "
              f"{'unchanged' if old == new else 'renumbered'}")
        if not dry_run:
            write(path, new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--style", choices=["bracket", "superscript"], default="bracket")
    ap.add_argument("--map-in", help="reuse a previously written _citation_number_map.json")
    args = ap.parse_args()

    md = read(MAIN_MD)
    body, lines, pre_i, entries, tail_i = parse_ref_section(md)
    assert md == body + "\n".join(lines), "reference-section split is not lossless"

    occ, spans, groups = scan(body)
    if not occ:
        sys.exit("no author-year citations found - already converted?")

    years = {}
    for _, _, tok, _ in occ:
        years.setdefault(token_surname(tok), token_year(tok))
    order = order_of(occ)
    matched = {s: match_entry(s, years[s], entries) for s in order}
    unmatched = [s for s, v in matched.items() if v is None]
    if unmatched:
        sys.exit(f"in-text forms with no unique reference entry: {unmatched}")
    ref_surnames = {s: v[0] for s, v in matched.items()}
    uncited = [n for _, n, _ in entries if n not in set(ref_surnames.values())]
    if uncited:
        sys.exit(f"reference entries never cited: {uncited}")

    if args.map_in:
        saved = json.load(open(args.map_in, encoding="utf-8"))
        num_of = {s: saved["surname_to_number"][s] for s in order}
        assert sorted(num_of.values()) == list(range(1, len(order) + 1))
    else:
        num_of = {s: i + 1 for i, s in enumerate(order)}

    print(f"in-text citation occurrences: {len(occ)} | parenthetical groups: {len(groups)}")
    print(f"distinct cited references: {len(order)} | reference-list entries: {len(entries)}")
    print("\nnew# | old# | in-text form")
    for s, new in sorted(num_of.items(), key=lambda kv: kv[1]):
        print(f"{new:4d} | {ref_surnames[s]:4d} | {s}")
    assert len(order) == len(entries)

    # --- in-text rewrite -----------------------------------------------------
    edits = []
    for gi, idxs in groups.items():
        a, b = spans[gi]
        inner = body[a:b]
        special = SPECIAL_GROUPS.get(inner)
        nums = sorted({num_of[token_surname(occ[i][2])] for i in idxs})
        if special:                                          # replaces the whole parenthesis
            edits.append((a - 1, b + 1, special.format(n=marker(nums, args.style))))
            continue
        toks = [(occ[i][0] - a, occ[i][1] - a, occ[i][2]) for i in idxs]
        rep = rewrite_group(inner, toks, num_of, args.style)
        sa, ea = a, b
        if rep is None:                                      # citations only -> drop the parenthesis
            rep = marker(nums, args.style)
            # ``a`` is the first character *inside* the parenthesis, so the character
            # preceding the removed '(' (and its ')' ) is at a - 2
            if a < 2 or body[a - 2] not in " \t\n":
                rep = " " + rep
            sa, ea = a - 1, b + 1
        edits.append((sa, ea, rep))
    for s, e, tok, gi in occ:
        if gi is None:                                       # narrative citation
            edits.append((s, e, narrative(tok, num_of[token_surname(tok)], args.style)))
    edits.sort()
    out, last = [], 0
    for s, e, rep in edits:
        assert s >= last, f"overlapping citation edits at {s}"
        out.append(body[last:s])
        out.append(rep)
        last = e
    out.append(body[last:])
    body2 = "".join(out)
    for key in SPECIAL_GROUPS:
        assert key in body, f"special-case group no longer present: {key!r}"
    residual = [t.group(0) for t in CITE_RE.finditer(body2) if is_citation(t.group(0))]
    if residual:
        sys.exit(f"residual author-year citations after rewrite: {residual}")

    # --- reference list reorder ---------------------------------------------
    by_old = {n: text for _, n, text in entries}
    new_entries = []
    for s, new in sorted(num_of.items(), key=lambda kv: kv[1]):
        assert ref_surnames[s] in by_old, s
        new_entries.append(f"{new}. {by_old[ref_surnames[s]]}")

    preamble = lines[pre_i]
    preamble = re.sub(r"\(entry \d+, ", "(", preamble)
    preamble = re.sub(r" \(entry \d+\)", "", preamble)
    if "ordered by first citation" not in preamble:
        note = PREAMBLE_ORDER_NOTE
        preamble = (preamble.rstrip()[:-1] + note + "*" if preamble.rstrip().endswith("*")
                    else preamble.rstrip() + note)
    assert not re.search(r"\bentry \d+", preamble), "entry numbers left in preamble"

    new_lines = list(lines)
    new_lines[pre_i] = preamble
    for k, _, _ in entries:
        new_lines[k] = ""
    for off, txt in enumerate(new_entries):
        new_lines[entries[0][0] + off] = txt
    assert len(new_entries) == len(entries)
    assert entries[0][0] + len(new_entries) - 1 == entries[-1][0], "reference list would not fit its block"
    out_md = body2 + "\n".join(new_lines)

    # --- other files carrying the same in-text forms ------------------------
    others = {}
    for path in TEXT_FILES[1:]:
        t = read(path)
        # legend text lists resources by name, so the author name is kept: "Finan et al. [16]"
        toks = [(m.start(), m.end(), m.group(0)) for m in CITE_RE.finditer(t) if is_citation(m.group(0))]
        for s, e, tok in reversed(toks):
            t = t[:s] + narrative(tok, num_of[token_surname(tok)], args.style) + t[e:]
        left = [tok for _, _, tok in toks if token_surname(tok) not in num_of]
        if left:
            sys.exit(f"{path}: citation form absent from the manuscript map: {left}")
        others[path] = (t, len(toks))

    old_to_new = {ref_surnames[s]: n for s, n in num_of.items()}

    if args.dry_run:
        print("\n--- DRY RUN, nothing written ---")
        print("citation edits (before -> after):")
        for s, e, rep in edits:
            print(f"  {body[s:e]}\n    -> {rep}")
        for path, (t, n) in others.items():
            print(f"{os.path.relpath(path, ROOT)}: {n} in-text marker(s)")
        print(f"{os.path.relpath(MAIN_MD, ROOT)}: body markers {len(occ)}, list entries {len(new_entries)}")
        reorder_downstream(old_to_new, num_of, dry_run=True)
        return

    write(MAIN_MD, out_md)
    for path, (t, n) in others.items():
        write(path, t)
    os.makedirs(os.path.dirname(MAP_JSON), exist_ok=True)
    json.dump({
        "converted": "2026-09-13",
        "style": args.style,
        "order": "first citation in text (NLM/ICMJE)",
        "surname_to_number": num_of,
        "old_number_to_new": {str(k): v for k, v in old_to_new.items()},
        "in_text_forms": {s: s for s in num_of},
        "rules": {
            "multiple": "ascending inside one bracket; >=3 consecutive collapse to a-b",
            "mixed_parenthesis": "citation replaced in place, separator before marker dropped",
            "citations_only_parenthesis": "parentheses dropped, bracket alone: (Cao et al. 2023) -> [11]",
            "exceptions": sorted(SPECIAL_GROUPS),
        },
    }, open(MAP_JSON, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"\nwritten: manuscript + {len(others)} legend/manifest file(s)")
    print("written:", os.path.relpath(MAP_JSON, ROOT))
    reorder_downstream(old_to_new, num_of, dry_run=False)
    print("NOTE: re-run analysis/scripts/_fetch_ref_pmids.py and then "
          "analysis/scripts/build_reference_pmid_list.py so the PMID ledger follows the new order.")


if __name__ == "__main__":
    main()
