from datetime import datetime, timedelta
from Bio import Entrez
import sqlite3
from pathlib import Path
import requests
import re
from html import unescape
from xml.etree import ElementTree as ET
import html as html_lib

# Metti una tua email vera (consigliato da NCBI)
Entrez.email = "tuamail@example.com"

# Query: stretta + fallback
QUERY_STRICT = '((ponatinib OR iclusig OR AP24534) OR (asciminib OR scemblix OR ABL001)) AND (("chronic myeloid leukemia") OR CML OR ("chronic myelogenous leukemia"))'
QUERY_FALLBACK = '(ponatinib OR asciminib) AND (CML OR "chronic myeloid leukemia")'

BASE_DIR = Path(__file__).parent
DB_PATH = BASE_DIR / "ponatinib_digest.sqlite"
OUT_DIR = BASE_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)
OUT_MD = OUT_DIR / "weekly_digest.md"
OUT_HTML = OUT_DIR / "weekly_digest.html"

EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPE_PMC_FULLTEXT_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"


def init_db():
    con = sqlite3.connect(DB_PATH)
    con.execute("""
        CREATE TABLE IF NOT EXISTS seen (
            pmid TEXT PRIMARY KEY,
            doi TEXT,
            title TEXT,
            first_seen TEXT
        )
    """)
    con.commit()
    return con


def pubmed_search_last_days(days=15, retmax=300):
    today = datetime.utcnow().date()
    start = today - timedelta(days=days)
    mindate = start.strftime("%Y/%m/%d")
    maxdate = today.strftime("%Y/%m/%d")

    def run_search(term):
        handle = Entrez.esearch(
            db="pubmed",
            term=term,
            mindate=mindate,
            maxdate=maxdate,
            datetype="edat",   # entry date: meglio per rassegna periodica
            retmax=retmax,
            sort="pub+date"
        )
        res = Entrez.read(handle)
        return res["IdList"]

    # Ricerca stretta
    ids_strict = run_search(QUERY_STRICT)
    print(f"Ricerca stretta: trovati {len(ids_strict)} PMID")

    # Fallback se zero
    if len(ids_strict) == 0:
        ids_fallback = run_search(QUERY_FALLBACK)
        print(f"Fallback: trovati {len(ids_fallback)} PMID")
        return ids_fallback

    return ids_strict


def pubmed_fetch_details(pmids):
    if not pmids:
        return []

    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
    records = Entrez.read(handle)

    articles = []
    for art in records["PubmedArticle"]:
        med = art["MedlineCitation"]
        pmid = str(med["PMID"])
        article = med["Article"]

        title = str(article.get("ArticleTitle", "")).strip()

        abstract = ""
        if "Abstract" in article and "AbstractText" in article["Abstract"]:
            abstract = " ".join([str(x) for x in article["Abstract"]["AbstractText"]]).strip()

        doi = ""
        if "ELocationID" in article:
            for eloc in article["ELocationID"]:
                if getattr(eloc, "attributes", {}).get("EIdType") == "doi":
                    doi = str(eloc)

        articles.append({
            "pmid": pmid,
            "doi": doi,
            "title": title,
            "abstract": abstract
        })
    return articles


def filter_new(con, articles):
    new = []
    for a in articles:
        cur = con.execute("SELECT 1 FROM seen WHERE pmid=?", (a["pmid"],))
        if cur.fetchone() is None:
            new.append(a)
    return new


def mark_seen(con, articles):
    now = datetime.utcnow().isoformat()
    for a in articles:
        con.execute(
            "INSERT OR IGNORE INTO seen(pmid, doi, title, first_seen) VALUES (?,?,?,?)",
            (a["pmid"], a["doi"], a["title"], now)
        )
    con.commit()


# -----------------------------
# Utility testo / riassunto gratis
# -----------------------------
def clean_text(text):
    if not text:
        return ""
    text = unescape(text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def split_sentences(text):
    text = clean_text(text)
    parts = re.split(r'(?<=[\.\!\?])\s+', text)
    return [p.strip() for p in parts if p and len(p.strip()) > 20]


def pick_best_sentences(text, max_sentences=3):
    """Riassunto estrattivo gratis: sceglie frasi utili."""
    sents = split_sentences(text)
    if not sents:
        return []

    keywords = [
        "method", "methods", "patients", "population", "cohort", "study", "trial",
        "results", "response", "survival", "safety", "adverse", "conclusion",
        "ponatinib", "asciminib", "cml", "chronic myeloid leukemia", "bcr-abl", "t315i"
    ]

    scored = []
    for i, s in enumerate(sents):
        sl = s.lower()
        score = 0
        score += sum(1 for k in keywords if k in sl)
        if re.search(r"\b\d+(\.\d+)?%?\b", s):
            score += 2
        if 40 <= len(s) <= 260:
            score += 2
        if i < 3:
            score += 1
        scored.append((score, i, s))

    scored.sort(reverse=True)
    top = sorted(scored[:max_sentences], key=lambda x: x[1])
    return [s for _, _, s in top]


def summarize_abstract_free(abstract):
    if not abstract or len(abstract.strip()) < 30:
        return [
            "Tipo di studio/popolazione: non riportato (abstract non disponibile).",
            "Risultato chiave: non riportato.",
            "Rilevanza clinica: non valutabile senza abstract."
        ]

    chosen = pick_best_sentences(abstract, max_sentences=3)

    labels = [
        "Tipo di studio/popolazione",
        "Risultato chiave",
        "Rilevanza clinica"
    ]
    bullets = []
    for i in range(3):
        if i < len(chosen):
            bullets.append(f"{labels[i]}: {chosen[i]}")
        else:
            bullets.append(f"{labels[i]}: non riportato.")
    return bullets


# -----------------------------
# Europe PMC OA full-text (gratis)
# -----------------------------
def europe_pmc_find_pmcid(pmid=None, doi=None):
    query_parts = []
    if pmid:
        query_parts.append(f"EXT_ID:{pmid} AND SRC:MED")
    if doi:
        query_parts.append(f'DOI:"{doi}"')

    if not query_parts:
        return None

    params = {
        "query": " OR ".join(query_parts),
        "format": "json",
        "pageSize": 5
    }

    try:
        r = requests.get(EUROPE_PMC_SEARCH, params=params, timeout=20)
        r.raise_for_status()
        data = r.json()
    except Exception:
        return None

    results = data.get("resultList", {}).get("result", [])

    # Preferisci OA vero
    for rec in results:
        pmcid = rec.get("pmcid")
        is_oa = rec.get("isOpenAccess")
        if pmcid and str(is_oa).upper() == "Y":
            return pmcid

    # fallback (alcuni record hanno pmcid ma flag non chiaro)
    for rec in results:
        pmcid = rec.get("pmcid")
        if pmcid:
            return pmcid

    return None


def extract_section_text_from_nxml(root, section_keywords):
    texts = []
    for sec in root.findall(".//sec"):
        title_el = sec.find("./title")
        sec_title = clean_text("".join(title_el.itertext())) if title_el is not None else ""
        sec_title_l = sec_title.lower()

        if any(k in sec_title_l for k in section_keywords):
            paragraphs = []
            for p in sec.findall(".//p"):
                t = clean_text("".join(p.itertext()))
                if t:
                    paragraphs.append(t)
            if paragraphs:
                texts.append((sec_title, " ".join(paragraphs)))
    return texts


def get_oa_fulltext_summary(pmcid):
    url = EUROPE_PMC_FULLTEXT_XML.format(pmcid=pmcid)

    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        xml_text = r.text
    except Exception:
        return None

    try:
        root = ET.fromstring(xml_text)
    except Exception:
        return None

    methods_sections = extract_section_text_from_nxml(root, ["method", "materials", "patients"])
    results_sections = extract_section_text_from_nxml(root, ["result"])
    discussion_sections = extract_section_text_from_nxml(root, ["discussion", "conclusion"])

    bullets = []

    if methods_sections:
        picks = pick_best_sentences(methods_sections[0][1], max_sentences=1)
        bullets.append(f"Metodi (full text OA): {picks[0] if picks else 'non riportato.'}")
    else:
        bullets.append("Metodi (full text OA): non riportato.")

    if results_sections:
        picks = pick_best_sentences(results_sections[0][1], max_sentences=2)
        bullets.append(f"Risultato principale (full text OA): {picks[0] if len(picks) >= 1 else 'non riportato.'}")
        bullets.append(f"Dettaglio utile (full text OA): {picks[1] if len(picks) >= 2 else 'non riportato.'}")
    else:
        bullets.append("Risultato principale (full text OA): non riportato.")
        bullets.append("Dettaglio utile (full text OA): non riportato.")

    if discussion_sections:
        picks = pick_best_sentences(discussion_sections[0][1], max_sentences=1)
        bullets.append(f"Messaggio clinico (full text OA): {picks[0] if picks else 'non riportato.'}")
    else:
        bullets.append("Messaggio clinico (full text OA): non riportato.")

    return bullets[:4]


# -----------------------------
# Tag semplice: Ponatinib / Asciminib / Entrambi
# -----------------------------
def tag_topic(title, abstract):
    txt = f"{title} {abstract}".lower()

    has_pona = any(k in txt for k in ["ponatinib", "iclusig", "ap24534"])
    has_asci = any(k in txt for k in ["asciminib", "scemblix", "abl001"])

    if has_pona and has_asci:
        return "Entrambi"
    elif has_pona:
        return "Ponatinib"
    elif has_asci:
        return "Asciminib"
    return "Altro"


def build_markdown(articles, days=15):
    today = datetime.utcnow().strftime("%Y-%m-%d")
    lines = []
    lines.append(f"# Ponatinib / Asciminib nella LMC – Rassegna quindicinale (ultimi {days} giorni)")
    lines.append(f"_Generata il: {today}_\n")

    if not articles:
        lines.append("✅ Nessun nuovo articolo nel periodo.\n")
        return "\n".join(lines)

    grouped = {"Ponatinib": [], "Asciminib": [], "Entrambi": [], "Altro": []}
    for a in articles:
        tag = tag_topic(a["title"], a.get("abstract", ""))
        grouped[tag].append(a)

    total = len(articles)
    lines.append(f"## Nuovi articoli: {total}\n")

    for section in ["Ponatinib", "Asciminib", "Entrambi", "Altro"]:
        section_articles = grouped.get(section, [])
        if not section_articles:
            continue

        lines.append(f"## {section} ({len(section_articles)})\n")

        for a in section_articles:
            pmid = a["pmid"]
            doi = a["doi"] if a["doi"] else "DOI non riportato"
            title = a["title"] if a["title"] else "(Titolo non disponibile)"
            abstract = a.get("abstract", "")
            pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            abs_bullets = summarize_abstract_free(abstract)

            pmcid = europe_pmc_find_pmcid(pmid=pmid, doi=a.get("doi"))
            oa_bullets = get_oa_fulltext_summary(pmcid) if pmcid else None

            lines.append(f"### {title}")
            lines.append(f"- PMID: {pmid} | {doi}")
            lines.append(f"- Link PubMed: {pubmed_link}")
            lines.append(f"- Open Access (Europe PMC/PMC): {'sì (' + pmcid + ')' if pmcid else 'non rilevato'}")
            lines.append("")
            lines.append("**Mini-riassunto (da abstract, estrattivo):**")
            for b in abs_bullets:
                lines.append(f"- {b}")

            if oa_bullets:
                lines.append("")
                lines.append("**Riassunto aggiuntivo (full text Open Access, estrattivo):**")
                for b in oa_bullets:
                    lines.append(f"- {b}")

            lines.append("\n---\n")

    return "\n".join(lines)


def markdown_to_simple_html(md_text: str) -> str:
    """
    Conversione semplice Markdown -> HTML (senza librerie extra).
    Non è perfetta, ma è leggibile.
    """
    lines = md_text.splitlines()
    html_lines = []
    in_ul = False

    def close_ul():
        nonlocal in_ul
        if in_ul:
            html_lines.append("</ul>")
            in_ul = False

    for raw in lines:
        line = raw.rstrip()

        if not line.strip():
            close_ul()
            html_lines.append("<p></p>")
            continue

        if line.strip() == "---":
            close_ul()
            html_lines.append("<hr>")
            continue

        if line.startswith("### "):
            close_ul()
            txt = html_lib.escape(line[4:])
            txt = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", txt)
            html_lines.append(f"<h3>{txt}</h3>")
            continue

        if line.startswith("## "):
            close_ul()
            txt = html_lib.escape(line[3:])
            txt = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", txt)
            html_lines.append(f"<h2>{txt}</h2>")
            continue

        if line.startswith("# "):
            close_ul()
            txt = html_lib.escape(line[2:])
            txt = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", txt)
            html_lines.append(f"<h1>{txt}</h1>")
            continue

        if line.startswith("- "):
            if not in_ul:
                html_lines.append("<ul>")
                in_ul = True

            item = html_lib.escape(line[2:])
            item = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", item)
            item = re.sub(
                r'(https://pubmed\.ncbi\.nlm\.nih\.gov/\d+/)',
                r'<a href="\1" target="_blank">\1</a>',
                item
            )
            html_lines.append(f"<li>{item}</li>")
            continue

        close_ul()
        txt = html_lib.escape(line)
        txt = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", txt)
        html_lines.append(f"<p>{txt}</p>")

    close_ul()

    body = "\n".join(html_lines)

    return f"""<!doctype html>
<html lang="it">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Rassegna Ponatinib / Asciminib CML</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      line-height: 1.5;
      margin: 24px auto;
      max-width: 900px;
      padding: 0 16px;
      color: #222;
      background: #fff;
    }}
    h1, h2, h3 {{ line-height: 1.2; }}
    h1 {{ border-bottom: 2px solid #ddd; padding-bottom: 8px; }}
    h2 {{ margin-top: 28px; border-bottom: 1px solid #eee; padding-bottom: 4px; }}
    h3 {{ margin-top: 22px; }}
    ul {{ margin-top: 6px; margin-bottom: 10px; }}
    li {{ margin-bottom: 6px; }}
    hr {{ margin: 22px 0; border: none; border-top: 1px solid #ddd; }}
    a {{ color: #0b57d0; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    p {{ margin: 8px 0; }}
  </style>
</head>
<body>
{body}
</body>
</html>
"""


def main(days=15):
    con = init_db()

    pmids = pubmed_search_last_days(days=days)
    print("Numero PMID trovati:", len(pmids))
    print("Primi PMID:", pmids[:10])

    articles = pubmed_fetch_details(pmids)
    new_articles = filter_new(con, articles)

    print(f"Trovati {len(articles)} articoli, nuovi {len(new_articles)}.")

    md = build_markdown(new_articles, days=days)
    OUT_MD.write_text(md, encoding="utf-8")

    html_text = markdown_to_simple_html(md)
    OUT_HTML.write_text(html_text, encoding="utf-8")

    mark_seen(con, new_articles)

    print(f"Creato file Markdown: {OUT_MD}")
    print(f"Creato file HTML: {OUT_HTML}")


if __name__ == "__main__":
    # TEST: 60 giorni per verificare che trovi articoli.
    # Quando confermi che funziona, rimetti 15.
    main(days=60)
