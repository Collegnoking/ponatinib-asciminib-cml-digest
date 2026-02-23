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

def detect_study_type(title: str, abstract: str) -> str:
    t = f"{title} {abstract}".lower()

    # guideline / consensus
    if "guideline" in t or "consensus" in t or "recommendation" in t:
        return "Linee guida / consenso"

    # systematic review / meta
    if "systematic review" in t or "meta-analysis" in t or "metaanalysis" in t:
        return "Revisione sistematica / meta-analisi"

    # review
    if "review" in t or "narrative review" in t or "special report" in t:
        return "Review"

    # trial
    if "randomized" in t or "phase ii" in t or "phase iii" in t or "clinical trial" in t or "open-label" in t:
        return "Trial clinico"

    # real-world / observational
    if "real-world" in t or "observational" in t or "cohort" in t or "registry" in t or "routine clinical practice" in t:
        return "Real-world / osservazionale"

    # pharmacovigilance
    if "eudravigilance" in t or "disproportionality" in t or "spontaneous reports" in t:
        return "Farmacovigilanza"

    # preclinico / in vitro
    if "cell" in t or "cells" in t or "k562" in t or "mouse" in t or "mice" in t or "in vitro" in t:
        return "Preclinico / laboratorio"

    # case reports
    if "case report" in t or "two cases" in t:
        return "Case report"

    return "Altro"


def impact_score(title: str, abstract: str) -> int:
    """Punteggio semplice e trasparente: serve per ordinare il Top 5."""
    t = f"{title} {abstract}".lower()
    stype = detect_study_type(title, abstract)

    score = 0
    # base per tipo studio
    base = {
        "Linee guida / consenso": 8,
        "Trial clinico": 7,
        "Real-world / osservazionale": 6,
        "Revisione sistematica / meta-analisi": 5,
        "Farmacovigilanza": 5,
        "Review": 4,
        "Case report": 2,
        "Preclinico / laboratorio": 2,
        "Altro": 3
    }
    score += base.get(stype, 3)

    # booster "clinicamente strategici"
    boosters = [
        ("tfr", 2), ("treatment-free", 2), ("treatment free", 2),
        ("mr4", 2), ("mr4.5", 2), ("deep molecular", 2),
        ("t315i", 2),
        ("arterial", 2), ("vascular", 2), ("cardio", 2), ("stroke", 2),
        ("anticoagul", 1),
        ("dose", 1), ("15 mg", 1), ("30 mg", 1), ("reduced dose", 1),
        ("resistance", 1), ("intoleran", 1),
        ("sequenc", 1), ("switch", 1)
    ]
    for key, add in boosters:
        if key in t:
            score += add

    # piccolo malus: troppo preclinico
    if stype == "Preclinico / laboratorio":
        score -= 1

    return max(score, 0)


def clinical_buckets(title: str, abstract: str) -> list[str]:
    """Ritorna 1+ categorie cliniche per raggruppare in modo 'strategico'."""
    t = f"{title} {abstract}".lower()
    buckets = []

    if "t315i" in t or "mutation" in t or "resistance" in t:
        buckets.append("Resistenza / T315I / mutazioni")

    if "tfr" in t or "treatment-free" in t or "treatment free" in t or "discontinu" in t or "stop" in t:
        buckets.append("TFR / stop terapia")

    if "dose" in t or "15 mg" in t or "30 mg" in t or "reduced dose" in t or "consolidation" in t:
        buckets.append("Dose strategy / consolidamento")

    if "vascular" in t or "arterial" in t or "cardio" in t or "stroke" in t or "pulmonary" in t or "thrombo" in t or "anticoagul" in t:
        buckets.append("Safety e comorbidità")

    if "first-line" in t or "second-line" in t or "third-line" in t or "sequenc" in t or "switch" in t:
        buckets.append("Sequencing / linee di terapia")

    if "real-world" in t or "routine clinical practice" in t or "registry" in t:
        buckets.append("Real-world (utilità pratica)")

    if "review" in t or "guideline" in t or "consensus" in t:
        buckets.append("Sintesi (review/linee guida)")

    if not buckets:
        buckets.append("Altro")

    return buckets


def make_executive_takeaways(articles: list[dict]) -> tuple[list[str], list[str], list[str]]:
    """
    Genera 3 takeaway + watchlist + action list, usando semplici conteggi per categoria.
    """
    # conteggi per bucket
    bucket_count = {}
    for a in articles:
        for b in clinical_buckets(a["title"], a.get("abstract", "")):
            bucket_count[b] = bucket_count.get(b, 0) + 1

    # ordina bucket per frequenza
    top_buckets = sorted(bucket_count.items(), key=lambda x: x[1], reverse=True)

    takeaways = []
    if top_buckets:
        # takeaway 1: bucket più frequente
        b1, n1 = top_buckets[0]
        takeaways.append(f"Nel periodo emergono soprattutto temi su **{b1}** (n={n1}).")
    if len(top_buckets) > 1:
        b2, n2 = top_buckets[1]
        takeaways.append(f"Secondo filone: **{b2}** (n={n2}).")
    if len(top_buckets) > 2:
        b3, n3 = top_buckets[2]
        takeaways.append(f"Terzo filone: **{b3}** (n={n3}).")

    while len(takeaways) < 3:
        takeaways.append("Nessun trend dominante: rassegna eterogenea nel periodo.")

    # watchlist (cosa monitorare)
    watchlist = []
    for key in ["TFR / stop terapia", "Safety e comorbidità", "Resistenza / T315I / mutazioni"]:
        if bucket_count.get(key, 0) > 0:
            watchlist.append(f"Monitorare evoluzione su **{key}** (nuovi segnali/dati).")
    if not watchlist:
        watchlist = ["Monitorare nuovi dati su efficacia e safety in CP-CML (linee successive)."]

    # action list (azioni pratiche)
    actions = []
    if bucket_count.get("Safety e comorbidità", 0) > 0:
        actions.append("Rivedere/rafforzare **monitoraggio cardiovascolare e interazioni** (es. anticoagulanti) nei pazienti in TKI.")
    if bucket_count.get("Dose strategy / consolidamento", 0) > 0 or bucket_count.get("TFR / stop terapia", 0) > 0:
        actions.append("Valutare criteri interni per **consolidamento / de-escalation / TFR** nei pazienti in deep response.")
    if bucket_count.get("Resistenza / T315I / mutazioni", 0) > 0:
        actions.append("Verificare percorso di **test mutazionale** e criteri di switch nelle resistenze (incl. T315I).")
    if not actions:
        actions = ["Nessuna azione immediata: mantenere sorveglianza e aggiornare la rassegna al prossimo ciclo."]

    return takeaways[:3], watchlist[:3], actions[:3]
def build_markdown(articles, days=15):
    today = datetime.utcnow().strftime("%Y-%m-%d")
    lines = []
    lines.append(f"# Ponatinib / Asciminib nella LMC – Report strategico quindicinale (ultimi {days} giorni)")
    lines.append(f"_Generata il: {today}_\n")

    if not articles:
        lines.append("✅ Nessun nuovo articolo nel periodo.\n")
        return "\n".join(lines)

    # --- Executive brief ---
    takeaways, watchlist, actions = make_executive_takeaways(articles)

    lines.append("## Executive brief\n")
    lines.append("**3 takeaway (cosa conta davvero):**")
    for t in takeaways:
        lines.append(f"- {t}")
    lines.append("")
    lines.append("**Watchlist (cosa tenere d’occhio):**")
    for w in watchlist:
        lines.append(f"- {w}")
    lines.append("")
    lines.append("**Action list (cosa fare in pratica):**")
    for a in actions:
        lines.append(f"- {a}")
    lines.append("\n---\n")

    # --- Top 5 per impatto ---
    scored = []
    for a in articles:
        s = impact_score(a["title"], a.get("abstract", ""))
        stype = detect_study_type(a["title"], a.get("abstract", ""))
        scored.append((s, stype, a))

    scored.sort(key=lambda x: x[0], reverse=True)
    top5 = scored[:5]

    lines.append("## Top 5 da leggere subito (ordinati per impatto)\n")
    for rank, (s, stype, a) in enumerate(top5, start=1):
        pmid = a["pmid"]
        doi = a["doi"] if a["doi"] else "DOI non riportato"
        title = a["title"] if a["title"] else "(Titolo non disponibile)"
        pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        buckets = ", ".join(clinical_buckets(title, a.get("abstract", "")))

        lines.append(f"### {rank}) {title}")
        lines.append(f"- **Impact score**: {s} | **Tipo**: {stype}")
        lines.append(f"- **Focus clinico**: {buckets}")
        lines.append(f"- PMID: {pmid} | {doi}")
        lines.append(f"- Link PubMed: {pubmed_link}")
        lines.append("")
        # mini-riassunto già esistente
        abs_bullets = summarize_abstract_free(a.get("abstract", ""))
        for b in abs_bullets:
            lines.append(f"- {b}")
        lines.append("\n---\n")

    # --- Raggruppamento per domande cliniche ---
    bucket_map = {}
    for a in articles:
        for b in clinical_buckets(a["title"], a.get("abstract", "")):
            bucket_map.setdefault(b, []).append(a)

    # Ordina bucket per numero articoli
    bucket_order = sorted(bucket_map.items(), key=lambda x: len(x[1]), reverse=True)

    lines.append("## Evidenze per domande cliniche\n")
    for bucket, items in bucket_order:
        lines.append(f"### {bucket} ({len(items)})\n")

        # ordina gli item del bucket per impatto
        items_sorted = sorted(items, key=lambda x: impact_score(x["title"], x.get("abstract", "")), reverse=True)

        for a in items_sorted:
            pmid = a["pmid"]
            doi = a["doi"] if a["doi"] else "DOI non riportato"
            title = a["title"] if a["title"] else "(Titolo non disponibile)"
            abstract = a.get("abstract", "")
            pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            # OA info + riassunti (già nel tuo flusso)
            pmcid = europe_pmc_find_pmcid(pmid=pmid, doi=a.get("doi"))
            oa_bullets = get_oa_fulltext_summary(pmcid) if pmcid else None

            stype = detect_study_type(title, abstract)
            s = impact_score(title, abstract)

            lines.append(f"#### {title}")
            lines.append(f"- Impact score: {s} | Tipo: {stype}")
            lines.append(f"- PMID: {pmid} | {doi}")
            lines.append(f"- Link PubMed: {pubmed_link}")
            lines.append(f"- Open Access (Europe PMC/PMC): {'sì (' + pmcid + ')' if pmcid else 'non rilevato'}")
            lines.append("")
            lines.append("**Mini-riassunto (da abstract, estrattivo):**")
            for b in summarize_abstract_free(abstract):
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
    main(days=30)
