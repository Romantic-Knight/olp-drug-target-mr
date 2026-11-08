# -*- coding: utf-8 -*-
"""Build an EndNote-importable RIS file for the 34 manuscript references.

1. Reads the verified PMIDs from results/tables/_ref_pmid_lookup.json (order 1-26).
2. efetch (single request) full PubMed XML for all 26 PMIDs.
3. Converts each PubmedArticle -> RIS (JOUR) with FULL author list, article title,
   full journal title, volume/issue/pages, year, DOI, PMID (AN).
4. Writes results/reports/reference_pmid_list.ris (UTF-8 with BOM for EndNote).

No fields are taken from memory; everything comes from the efetch XML.
"""
import json
import os
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SRC = os.path.join("results", "tables", "_ref_pmid_lookup.json")
OUT = os.path.join("results", "reports", "reference_pmid_list.ris")


def fetch_xml(pmids):
    url = (EUTILS + "/efetch.fcgi?db=pubmed&retmode=xml&id=" +
           ",".join(pmids) + "&tool=wisp-science-OLP-MR")
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8")


def strip_ns(root):
    for el in root.iter():
        if el.tag.startswith("{"):
            el.tag = el.tag.split("}", 1)[1]
    return root


def txt(el):
    if el is None:
        return ""
    return "".join(el.itertext()).strip()


def child(el, path):
    if el is None:
        return None
    return el.find(path)


def doid_of(rec, art):
    """DOI from ELocationID or PubmedData ArticleIdList (both under PubmedArticle)."""
    out = ""
    el = art.find("ELocationID")
    if el is not None and el.get("EIdType") == "doi":
        out = txt(el)
    if not out:
        pdata = rec.find("PubmedData")
        if pdata is not None:
            for aid in pdata.findall("ArticleIdList/ArticleId"):
                if aid.get("IdType") == "doi":
                    out = txt(aid)
                    break
    return out.lower()


def split_pages(pg):
    pg = re.sub(r"\s+", "", pg or "")
    if not pg:
        return None, None
    parts = pg.split("-")
    if len(parts) == 2 and parts[1]:
        return parts[0], parts[1]
    return parts[0], None


def parse_articles(xml_text):
    root = strip_ns(ET.fromstring(xml_text))
    arts = {}
    for rec in root.findall("PubmedArticle"):
        pmid = txt(child(rec, "MedlineCitation/PMID"))
        if not pmid:
            continue
        med = rec.find("MedlineCitation")
        art = med.find("Article") if med is not None else None
        if art is None:
            continue
        issue = child(art, "Journal/JournalIssue")
        pubdate = child(issue, "PubDate") if issue is not None else None
        year = txt(child(pubdate, "Year"))
        if not year:
            md = txt(child(pubdate, "MedlineDate"))
            year = md[:4]
        authors = []
        alist = art.find("AuthorList")
        if alist is not None:
            for a in alist.findall("Author"):
                ln = txt(child(a, "LastName"))
                coll = txt(child(a, "CollectiveName"))
                if ln:
                    fn = txt(child(a, "ForeName")) or txt(child(a, "Initials"))
                    authors.append(f"{ln}, {fn}" if fn else ln)
                elif coll:
                    authors.append(coll)
        pages = txt(child(art, "Pagination/MedlinePgn")) or txt(child(art, "Pagination"))
        sp, ep = split_pages(pages)
        rec2 = {
            "pmid": pmid,
            "authors": authors,
            "title": txt(child(art, "ArticleTitle")),
            "journal": txt(child(art, "Journal/Title")) or txt(child(art, "Journal/ISOAbbreviation")),
            "volume": txt(child(issue, "Volume")) if issue is not None else "",
            "issue": txt(child(issue, "Issue")) if issue is not None else "",
            "year": year,
            "sp": sp,
            "ep": ep,
            "doi": doid_of(rec, art),
        }
        arts[pmid] = rec2
    return arts


def ris_lines(r):
    lines = ["TY  - JOUR"]
    for au in r["authors"]:
        lines.append(f"AU  - {au}")
    lines.append(f"TI  - {r['title'].rstrip('.')}")
    if r["journal"]:
        lines.append(f"JO  - {r['journal']}")
    if r["volume"]:
        lines.append(f"VL  - {r['volume']}")
    if r["issue"]:
        lines.append(f"IS  - {r['issue']}")
    if r["sp"]:
        lines.append(f"SP  - {r['sp']}")
        if r["ep"] and r["ep"] != r["sp"]:
            lines.append(f"EP  - {r['ep']}")
    if r["year"]:
        lines.append(f"PY  - {r['year']}")
    if r["doi"]:
        lines.append(f"DO  - {r['doi']}")
    lines.append(f"AN  - {r['pmid']}")
    lines.append("ER  - ")
    return lines


def main():
    recs = json.load(open(SRC, encoding="utf-8"))
    pmids = [str(r["pmid"]) for r in recs]
    # count follows the reference list (34 since 2026-09-13); records keep the
    # manuscript's first-citation order, which is the order EndNote will import
    assert len(pmids) == 34 and all(pmids), "lookup json must hold the 34 verified PMIDs"

    xml_text = fetch_xml(pmids)
    arts = parse_articles(xml_text)
    missing = [p for p in pmids if p not in arts]
    if missing:
        raise SystemExit(f"efetch missing PMIDs: {missing}")

    blocks = []
    for r in recs:
        art = arts[str(r["pmid"])]
        lines = ris_lines(art)
        if not art["authors"]:
            print(f"[warn] PMID {r['pmid']}: no authors parsed")
        blocks.append("\n".join(lines))

    with open(OUT, "w", encoding="utf-8-sig", newline="\n") as f:
        f.write("\n\n".join(blocks) + "\n")
    print("written:", OUT)
    print("records:", sum(b.startswith("TY  - ") for b in blocks))


if __name__ == "__main__":
    main()
