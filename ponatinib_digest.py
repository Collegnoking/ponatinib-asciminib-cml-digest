from datetime import datetime, timedelta
from Bio import Entrez
import sqlite3
from pathlib import Path
import requests
import re
from html import unescape
from xml.etree import ElementTree as ET
import html as html_lib

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, HRFlowable,
    Table, TableStyle, KeepTogether
)
from reportlab.platypus.flowables import HRFlowable


# =========================
# CONFIG
# =========================
Entrez.email = "nicola.battaglia71@gmail.com"  # <-- CAMBIA QUI

QUERY_STRICT = '((ponatinib OR iclusig OR AP24534) OR (asciminib OR scemblix OR ABL001)) AND (("chronic myeloid leukemia") OR CML OR ("chronic myelogenous leukemia"))'
QUERY_FALLBACK = '(ponatinib OR asciminib) AND (CML OR "chronic myeloid leukemia")'

BASE_DIR = Path(__file__).parent
DB_PATH = BASE_DIR / "ponatinib_digest.sqlite"

OUT_DIR = BASE_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)

OUT_MD = OUT_DIR / "weekly_digest.md"
OUT_HTML = OUT_DIR / "weekly_digest.html"
OUT_PDF = OUT_DIR / "weekly_digest.pdf"

EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPE_PMC_FULLTEXT_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"


# =========================
# DB
# =========================
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


# =========================
# PubMed
# =========================
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
            datetype="edat",
            retmax=retmax,
            sort="pub+date"
        )
        res = Entrez.read(handle)
        return res["IdList"]

    ids_strict = run_search(QUERY_STRICT)
    print(f"Ricerca stretta: trovati {len(ids_strict)} PMID")

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

        journal = ""
        try:
            if "Journal" in article and "Title" in article["Journal"]:
                journal = str(article["Journal"]["Title"]).strip()
        except Exception:
            journal = ""

        articles.append({
            "pmid": pmid,
            "doi": doi,
            "title": title,
            "abstract": abstract,
            "journal": journal
        })

    return articles


# =========================
# Text utilities
# =========================
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
    labels = ["Tipo di studio/popolazione", "Risultato chiave", "Rilevanza clinica"]

    bullets = []
    for i in range(3):
        if i < len(chosen):
            bullets.append(f"{labels[i]}: {chosen[i]}")
        else:
            bullets.append(f"{labels[i]}: non riportato.")
    return bullets


# =========================
# Europe PMC OA full-text
# =========================
def europe_pmc_find_pmcid(pmid=None, doi=None):
    query_parts = []
    if pmid:
        query_parts.append(f"EXT_ID:{pmid} AND SRC:MED")
    if doi:
        query_parts.append(f'DOI:"{doi}"')

    if not query_parts:
        return None

    params = {"query": " OR ".join(query_parts), "format": "json", "pageSize": 5}

    try:
        r = requests.get(EUROPE_PMC_SEARCH, params=params, timeout=20)
        r.raise_for_status()
        data = r.json()
    except Exception:
        return None

    results = data.get("resultList", {}).get("result", [])

    for rec in results:
        pmcid = rec.get("pmcid")
        is_oa = rec.get("isOpenAccess")
        if pmcid and str(is_oa).upper() == "Y":
            return pmcid

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


# =========================
# Scoring / buckets
# =========================
def detect_study_type(title: str, abstract: str) -> str:
    t = f"{title} {abstract}".lower()

    if "guideline" in t or "consensus" in t or "recommendation" in t:
        return "Linee guida / consenso"
    if "systematic review" in t or "meta-analysis" in t or "metaanalysis" in t:
        return "Revisione sistematica / meta-analisi"
    if "review" in t or "narrative review" in t or "special report" in t:
        return "Review"
    if "randomized" in t or "phase ii" in t or "phase iii" in t or "clinical trial" in t or "open-label" in t:
        return "Trial clinico"
    if "real-world" in t or "observational" in t or "cohort" in t or "registry" in t or "routine clinical practice" in t:
        return "Real-world / osservazionale"
    if "eudravigilance" in t or "disproportionality" in t or "spontaneous reports" in t:
        return "Farmacovigilanza"
    if "case report" in t or "two cases" in t:
        return "Case report"
    if "cell" in t or "cells" in t or "k562" in t or "mouse" in t or "mice" in t or "in vitro" in t:
        return "Preclinico / laboratorio"

    return "Altro"


def impact_score(title: str, abstract: str) -> int:
    t = f"{title} {abstract}".lower()
    stype = detect_study_type(title, abstract)

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
    score = base.get(stype, 3)

    boosters = [
        ("tfr", 2), ("treatment-free", 2), ("treatment free", 2),
        ("mr4", 2), ("mr4.5", 2), ("deep molecular", 2),
        ("t315i", 2),
        ("arterial", 2), ("vascular", 2), ("cardio", 2), ("stroke", 2),
        ("anticoagul", 1),
        ("dose", 1), ("15 mg", 1), ("30 mg", 1), ("reduced dose", 1), ("consolidation", 1),
        ("resistance", 1), ("intoleran", 1),
        ("sequenc", 1), ("switch", 1)
    ]
    for key, add in boosters:
        if key in t:
            score += add

    if stype == "Preclinico / laboratorio":
        score -= 1

    return max(score, 0)


def clinical_buckets(title: str, abstract: str) -> list[str]:
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
    bucket_count = {}
    for a in articles:
        for b in clinical_buckets(a["title"], a.get("abstract", "")):
            bucket_count[b] = bucket_count.get(b, 0) + 1

    top_buckets = sorted(bucket_count.items(), key=lambda x: x[1], reverse=True)

    takeaways = []
    if top_buckets:
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

    watchlist = []
    for key in ["TFR / stop terapia", "Safety e comorbidità", "Resistenza / T315I / mutazioni"]:
        if bucket_count.get(key, 0) > 0:
            watchlist.append(f"Monitorare evoluzione su **{key}** (nuovi segnali/dati).")
    if not watchlist:
        watchlist = ["Monitorare nuovi dati su efficacia e safety in CP-CML (linee successive)."]

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


# =========================
# Report Markdown
# =========================
def build_markdown(articles, days=15):
    today = datetime.utcnow().strftime("%Y-%m-%d")
    lines = []
    lines.append(f"# Ponatinib / Asciminib nella LMC – Newsletter interna (ultimi {days} giorni)")
    lines.append(f"_Generata il: {today}_\n")

    if not articles:
        lines.append("✅ Nessun nuovo articolo nel periodo.\n")
        return "\n".join(lines)

    takeaways, watchlist, actions = make_executive_takeaways(articles)

    lines.append("## Newsletter interna – sintesi operativa\n")
    lines.append("**Messaggi chiave (da leggere in 60 secondi):**")
    for t in takeaways:
        lines.append(f"- {t}")
    lines.append("")
    lines.append("**Cosa cambia per noi (azioni pratiche):**")
    for a in actions:
        lines.append(f"- {a}")
    lines.append("")
    lines.append("**Safety radar (segnali/interazioni da monitorare):**")
    for w in watchlist:
        lines.append(f"- {w}")
    lines.append("\n---\n")

    scored = []
    for a in articles:
        s = impact_score(a["title"], a.get("abstract", ""))
        stype = detect_study_type(a["title"], a.get("abstract", ""))
        scored.append((s, stype, a))
    scored.sort(key=lambda x: x[0], reverse=True)

    lines.append("## Top 5 da leggere subito (ordinati per impatto)\n")
    for rank, (s, stype, a) in enumerate(scored[:5], start=1):
        pmid = a["pmid"]
        doi = a["doi"] if a["doi"] else "DOI non riportato"
        title = a["title"] if a["title"] else "(Titolo non disponibile)"
        pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        buckets = ", ".join(clinical_buckets(title, a.get("abstract", "")))

        journal = a.get("journal", "").strip() or "Rivista non riportata"

        lines.append(f"### {rank}) {title}")
        lines.append(f"- **Rivista**: {journal} | **Tipo**: {stype}")
        lines.append(f"- **Impact score**: {s}")
        lines.append(f"- **Focus clinico**: {buckets}")
        lines.append(f"- PMID: {pmid} | {doi}")
        lines.append(f"- Link PubMed: {pubmed_link}")
        lines.append("")
        for b in summarize_abstract_free(a.get("abstract", "")):
            lines.append(f"- {b}")
        lines.append("\n---\n")

    bucket_map = {}
    for a in articles:
        for b in clinical_buckets(a["title"], a.get("abstract", "")):
            bucket_map.setdefault(b, []).append(a)

    bucket_order = sorted(bucket_map.items(), key=lambda x: len(x[1]), reverse=True)

    lines.append("## Evidenze per domande cliniche\n")
    for bucket, items in bucket_order:
        lines.append(f"### {bucket} ({len(items)})\n")
        items_sorted = sorted(items, key=lambda x: impact_score(x["title"], x.get("abstract", "")), reverse=True)

        for a in items_sorted:
            pmid = a["pmid"]
            doi = a["doi"] if a["doi"] else "DOI non riportato"
            title = a["title"] if a["title"] else "(Titolo non disponibile)"
            abstract = a.get("abstract", "")
            pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            journal = a.get("journal", "").strip() or "Rivista non riportata"
            stype = detect_study_type(title, abstract)
            s = impact_score(title, abstract)

            pmcid = europe_pmc_find_pmcid(pmid=pmid, doi=a.get("doi"))
            oa_bullets = get_oa_fulltext_summary(pmcid) if pmcid else None

            lines.append(f"### {title}")
            lines.append(f"- **Rivista**: {journal} | **Tipo**: {stype}")
            lines.append(f"- Impact score: {s}")
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


# =========================
# HTML – Design raffinato
# =========================
def markdown_to_simple_html(md_text: str) -> str:
    try:
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
                continue

            if line.strip() == "---":
                close_ul()
                html_lines.append('<div class="divider"></div>')
                continue

            if line.startswith("### "):
                close_ul()
                txt = html_lib.escape(line[4:])
                txt = re.sub(r"\*\*(.*?)\*\*", r"<strong>\1</strong>", txt)
                html_lines.append(f'<h3>{txt}</h3>')
                continue

            if line.startswith("## "):
                close_ul()
                txt = html_lib.escape(line[3:])
                txt = re.sub(r"\*\*(.*?)\*\*", r"<strong>\1</strong>", txt)
                html_lines.append(f'<h2>{txt}</h2>')
                continue

            if line.startswith("# "):
                close_ul()
                txt = html_lib.escape(line[2:])
                html_lines.append(f'<h1>{txt}</h1>')
                continue

            # Badge rivista / tipo
            if line.startswith("- **Rivista**:"):
                close_ul()
                m = re.search(r"\*\*Rivista\*\*:\s*(.*?)\s*\|\s*\*\*Tipo\*\*:\s*(.*)$", line)
                if m:
                    journal = html_lib.escape(m.group(1).strip())
                    stype_raw = m.group(2).strip().lower()

                    if "trial" in stype_raw:
                        cls, label, icon = "badge-trial", "Trial clinico", "🔬"
                    elif "real-world" in stype_raw or "osserv" in stype_raw:
                        cls, label, icon = "badge-rwe", "Real-world", "🏥"
                    elif "farmacovigil" in stype_raw or "safety" in stype_raw:
                        cls, label, icon = "badge-safety", "Farmacovigilanza", "⚠️"
                    elif "meta" in stype_raw:
                        cls, label, icon = "badge-meta", "Meta-analisi", "📊"
                    elif "review" in stype_raw:
                        cls, label, icon = "badge-review", "Review", "📖"
                    elif "linee guida" in stype_raw or "consenso" in stype_raw:
                        cls, label, icon = "badge-guideline", "Linee guida", "📋"
                    elif "case report" in stype_raw:
                        cls, label, icon = "badge-case", "Case Report", "📝"
                    else:
                        cls, label, icon = "badge-other", "Altro", "📄"

                    html_lines.append(
                        f'<div class="meta-row">'
                        f'<span class="badge {cls}">{icon} {label}</span>'
                        f'<span class="journal-name">📰 {journal}</span>'
                        f'</div>'
                    )
                    continue

            if line.startswith("- "):
                if not in_ul:
                    html_lines.append('<ul>')
                    in_ul = True
                item = html_lib.escape(line[2:])
                item = re.sub(r"\*\*(.*?)\*\*", r"<strong>\1</strong>", item)
                item = re.sub(
                    r'(https://pubmed\.ncbi\.nlm\.nih\.gov/\d+/)',
                    r'<a href="\1" target="_blank" rel="noopener">\1</a>',
                    item
                )
                # Evidenzia etichette tipo "Risultato chiave:", "Tipo di studio:"
                item = re.sub(
                    r'^(Tipo di studio/popolazione|Risultato chiave|Rilevanza clinica|'
                    r'Metodi \(full text OA\)|Risultato principale \(full text OA\)|'
                    r'Dettaglio utile \(full text OA\)|Messaggio clinico \(full text OA\))'
                    r'(:)',
                    r'<span class="bullet-label">\1\2</span>',
                    item
                )
                html_lines.append(f'<li>{item}</li>')
                continue

            close_ul()
            txt = html_lib.escape(line)
            txt = re.sub(r"\*\*(.*?)\*\*", r"<strong>\1</strong>", txt)
            # Gestione corsivo "_testo_"
            txt = re.sub(r"_(.*?)_", r"<em>\1</em>", txt)
            html_lines.append(f'<p>{txt}</p>')

        close_ul()
        body = "\n".join(html_lines)

        today_str = datetime.utcnow().strftime("%d %B %Y")

        return f"""<!doctype html>
<html lang="it">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>LMC Newsletter – Ponatinib / Asciminib</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Serif:ital,wght@0,400;0,600;1,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono&display=swap" rel="stylesheet">
  <style>
    /* ── Reset & base ── */
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

    :root {{
      --ink:        #1a1a2e;
      --ink-light:  #4a4a6a;
      --accent:     #1b4f72;
      --accent2:    #0e6655;
      --surface:    #f8f9fb;
      --card:       #ffffff;
      --border:     #dde3ec;
      --red:        #922b21;
      --amber:      #7d6608;
      --green:      #0e6655;
      --blue:       #1a5276;
      --purple:     #6c3483;
      --teal:       #0e7490;
      --radius:     10px;
      --shadow:     0 2px 12px rgba(27,79,114,.10);
    }}

    body {{
      font-family: 'IBM Plex Sans', sans-serif;
      font-size: 15px;
      line-height: 1.7;
      color: var(--ink);
      background: var(--surface);
      padding: 0;
    }}

    /* ── Header banner ── */
    .site-header {{
      background: var(--accent);
      color: #fff;
      padding: 36px 40px 32px;
      position: relative;
      overflow: hidden;
    }}
    .site-header::after {{
      content: '';
      position: absolute;
      right: -60px; top: -60px;
      width: 260px; height: 260px;
      border-radius: 50%;
      background: rgba(255,255,255,.06);
    }}
    .site-header .label {{
      font-size: 11px;
      letter-spacing: .12em;
      text-transform: uppercase;
      opacity: .75;
      margin-bottom: 8px;
    }}
    .site-header h1 {{
      font-family: 'IBM Plex Serif', serif;
      font-size: 26px;
      font-weight: 600;
      line-height: 1.3;
      border: none;
      padding: 0;
      color: #fff;
    }}
    .site-header .date-badge {{
      display: inline-block;
      margin-top: 14px;
      background: rgba(255,255,255,.18);
      border: 1px solid rgba(255,255,255,.3);
      border-radius: 999px;
      padding: 4px 14px;
      font-size: 12px;
      letter-spacing: .04em;
    }}

    /* ── Layout ── */
    .wrapper {{
      max-width: 860px;
      margin: 0 auto;
      padding: 32px 20px 60px;
    }}

    /* ── Section headings ── */
    h2 {{
      font-family: 'IBM Plex Serif', serif;
      font-size: 20px;
      font-weight: 600;
      color: var(--accent);
      margin: 44px 0 16px;
      padding-bottom: 10px;
      border-bottom: 2px solid var(--accent);
      letter-spacing: -.01em;
    }}
    h3 {{
      font-family: 'IBM Plex Sans', sans-serif;
      font-size: 15px;
      font-weight: 600;
      color: var(--ink);
      margin: 28px 0 10px;
      line-height: 1.4;
    }}

    /* ── Callout box (sintesi operativa) ── */
    .wrapper > h2:first-of-type + * {{ }}  /* override se serve */

    /* Rileva il blocco sintesi avvolto nel div "callout" via JS → no JS, usiamo stili diretti */
    /* Usiamo una classe aggiunta nel rendering */
    .callout {{
      background: #eaf3fb;
      border-left: 4px solid var(--accent);
      border-radius: 0 var(--radius) var(--radius) 0;
      padding: 20px 24px;
      margin: 16px 0 28px;
    }}
    .callout-warning {{
      background: #fef9e7;
      border-left-color: #d4ac0d;
    }}
    .callout-danger {{
      background: #fdf2f0;
      border-left-color: var(--red);
    }}
    .callout p {{
      font-weight: 600;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: .06em;
      color: var(--accent);
      margin-bottom: 10px;
    }}
    .callout-warning p {{ color: var(--amber); }}
    .callout-danger p {{ color: var(--red); }}

    /* ── Cards articolo ── */
    .article-card {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: var(--radius);
      padding: 20px 22px 18px;
      margin: 0 0 20px;
      box-shadow: var(--shadow);
      transition: box-shadow .2s;
    }}
    .article-card:hover {{
      box-shadow: 0 4px 20px rgba(27,79,114,.16);
    }}

    /* ── Meta row (badge + journal) ── */
    .meta-row {{
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 8px;
      margin: 8px 0 12px;
    }}
    .badge {{
      display: inline-flex;
      align-items: center;
      gap: 4px;
      padding: 3px 10px;
      border-radius: 999px;
      font-size: 11.5px;
      font-weight: 600;
      letter-spacing: .03em;
      white-space: nowrap;
    }}
    .badge-trial    {{ background: #dbeafe; color: #1e3a8a; }}
    .badge-rwe      {{ background: #d1fae5; color: #064e3b; }}
    .badge-safety   {{ background: #fee2e2; color: #7f1d1d; }}
    .badge-meta     {{ background: #ede9fe; color: #4c1d95; }}
    .badge-review   {{ background: #fef3c7; color: #78350f; }}
    .badge-guideline{{ background: #e0f2fe; color: #0c4a6e; }}
    .badge-case     {{ background: #f3e8ff; color: #581c87; }}
    .badge-other    {{ background: #f1f5f9; color: #475569; }}
    .journal-name {{
      font-size: 12px;
      color: var(--ink-light);
      font-style: italic;
    }}

    /* ── Lists ── */
    ul {{
      padding-left: 20px;
      margin: 8px 0 14px;
    }}
    li {{
      margin-bottom: 7px;
      color: var(--ink);
      font-size: 14px;
    }}
    .bullet-label {{
      font-weight: 600;
      color: var(--accent);
    }}

    /* ── Links ── */
    a {{
      color: var(--accent);
      text-decoration: underline;
      text-underline-offset: 2px;
      word-break: break-all;
      font-family: 'IBM Plex Mono', monospace;
      font-size: 12.5px;
    }}
    a:hover {{ color: var(--accent2); }}

    /* ── Divider ── */
    .divider {{
      height: 1px;
      background: linear-gradient(to right, var(--border), transparent);
      margin: 28px 0;
    }}

    /* ── Paragraphs ── */
    p {{
      margin: 8px 0;
      font-size: 14.5px;
      color: var(--ink);
    }}
    em {{ color: var(--ink-light); font-size: 13px; }}

    /* ── Responsive ── */
    @media (max-width: 600px) {{
      .site-header {{ padding: 24px 18px; }}
      .site-header h1 {{ font-size: 20px; }}
      .wrapper {{ padding: 20px 14px 40px; }}
    }}
  </style>
</head>
<body>

<header class="site-header">
  <div class="label">Newsletter Oncologica Interna</div>
  <h1>Ponatinib / Asciminib nella LMC</h1>
  <span class="date-badge">📅 Aggiornamento del {today_str}</span>
</header>

<div class="wrapper">
{body}
</div>

</body>
</html>
"""
    except Exception as e:
        return "<html><body><pre>" + html_lib.escape(md_text) + "</pre></body></html>"


# =========================
# PDF – ReportLab Platypus
# =========================

# Palette colori
C_ACCENT   = colors.HexColor("#1b4f72")
C_ACCENT2  = colors.HexColor("#0e6655")
C_SURFACE  = colors.HexColor("#f0f4f8")
C_BORDER   = colors.HexColor("#dde3ec")
C_INK      = colors.HexColor("#1a1a2e")
C_INK_LIGHT= colors.HexColor("#4a4a6a")
C_RED      = colors.HexColor("#922b21")
C_AMBER    = colors.HexColor("#7d6608")
C_WHITE    = colors.white

# Badge type → color
BADGE_COLORS = {
    "trial":       colors.HexColor("#dbeafe"),
    "real-world":  colors.HexColor("#d1fae5"),
    "osserv":      colors.HexColor("#d1fae5"),
    "farmacovigil":colors.HexColor("#fee2e2"),
    "safety":      colors.HexColor("#fee2e2"),
    "meta":        colors.HexColor("#ede9fe"),
    "review":      colors.HexColor("#fef3c7"),
    "linee guida": colors.HexColor("#e0f2fe"),
    "consenso":    colors.HexColor("#e0f2fe"),
    "case report": colors.HexColor("#f3e8ff"),
}

def stype_badge_color(stype_raw: str):
    s = stype_raw.lower()
    for k, c in BADGE_COLORS.items():
        if k in s:
            return c
    return colors.HexColor("#f1f5f9")


def build_pdf_styles():
    base = getSampleStyleSheet()

    styles = {}

    styles["Title"] = ParagraphStyle(
        "Title",
        fontName="Helvetica-Bold",
        fontSize=22,
        leading=28,
        textColor=C_WHITE,
        spaceAfter=4,
    )
    styles["Subtitle"] = ParagraphStyle(
        "Subtitle",
        fontName="Helvetica",
        fontSize=11,
        leading=14,
        textColor=colors.HexColor("#cce0f5"),
        spaceAfter=0,
    )
    styles["H2"] = ParagraphStyle(
        "H2",
        fontName="Helvetica-Bold",
        fontSize=14,
        leading=18,
        textColor=C_ACCENT,
        spaceBefore=24,
        spaceAfter=8,
        borderPadding=(0, 0, 4, 0),
    )
    styles["H3"] = ParagraphStyle(
        "H3",
        fontName="Helvetica-Bold",
        fontSize=11,
        leading=15,
        textColor=C_INK,
        spaceBefore=14,
        spaceAfter=4,
    )
    styles["Body"] = ParagraphStyle(
        "Body",
        fontName="Helvetica",
        fontSize=10,
        leading=14,
        textColor=C_INK,
        spaceAfter=4,
    )
    styles["BulletItem"] = ParagraphStyle(
        "BulletItem",
        fontName="Helvetica",
        fontSize=9.5,
        leading=13,
        textColor=C_INK,
        leftIndent=12,
        spaceAfter=3,
        bulletIndent=0,
        bulletFontName="Helvetica",
        bulletFontSize=9,
        bulletText="•",
    )
    styles["BulletLabel"] = ParagraphStyle(
        "BulletLabel",
        fontName="Helvetica-Bold",
        fontSize=9.5,
        leading=13,
        textColor=C_ACCENT,
        leftIndent=12,
        spaceAfter=3,
    )
    styles["Meta"] = ParagraphStyle(
        "Meta",
        fontName="Helvetica-Oblique",
        fontSize=8.5,
        leading=12,
        textColor=C_INK_LIGHT,
        spaceAfter=6,
    )
    styles["CalloutTitle"] = ParagraphStyle(
        "CalloutTitle",
        fontName="Helvetica-Bold",
        fontSize=9,
        leading=12,
        textColor=C_ACCENT,
        spaceAfter=4,
    )
    styles["CalloutBody"] = ParagraphStyle(
        "CalloutBody",
        fontName="Helvetica",
        fontSize=9,
        leading=13,
        textColor=C_INK,
        leftIndent=8,
    )
    styles["Link"] = ParagraphStyle(
        "Link",
        fontName="Courier",
        fontSize=8,
        leading=11,
        textColor=C_ACCENT,
        spaceAfter=2,
    )

    return styles


def make_callout_table(title: str, items: list[str], color=None):
    """Crea un box colorato con titolo e bullet list."""
    bg = color or colors.HexColor("#eaf3fb")
    styles = build_pdf_styles()

    content = [Paragraph(title, styles["CalloutTitle"])]
    for item in items:
        item_clean = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", item)
        content.append(Paragraph(f"• {item_clean}", styles["CalloutBody"]))

    inner = Table([[c] for c in content], colWidths=["100%"])
    inner.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), bg),
        ("LEFTPADDING",  (0, 0), (-1, -1), 14),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
        ("TOPPADDING",   (0, 0), (0, 0), 10),
        ("BOTTOMPADDING",(0, -1),(-1,-1),  10),
        ("TOPPADDING",   (0, 1), (-1, -1), 2),
        ("BOTTOMPADDING",(0, 0), (-1, -2), 2),
        ("LINEAFTER",    (0, 0), (0, -1), 4, C_ACCENT),
        ("BOX",          (0, 0), (-1, -1), 0.5, C_BORDER),
        ("ROWBACKGROUNDS",(0,0),(-1,-1),[bg]),
    ]))
    return inner


def make_article_card(title: str, journal: str, stype: str, score: int,
                      meta_lines: list[str], bullet_lines: list[str],
                      styles: dict) -> Table:
    """Crea una card PDF per un singolo articolo."""
    badge_bg = stype_badge_color(stype)

    rows = []

    # Riga titolo
    rows.append([Paragraph(f"<b>{html_lib.escape(title)}</b>", styles["H3"])])

    # Riga badge
    badge_cell = Table(
        [[
            Paragraph(f"<b>{html_lib.escape(stype)}</b>", ParagraphStyle(
                "Badge",
                fontName="Helvetica-Bold", fontSize=8,
                leading=10, textColor=C_INK,
                backColor=badge_bg,
                borderPadding=(2, 6, 2, 6),
            )),
            Paragraph(f"<i>{html_lib.escape(journal)}</i>", styles["Meta"]),
            Paragraph(f"Score: <b>{score}</b>", styles["Meta"]),
        ]],
        colWidths=[None, None, 60],
    )
    badge_cell.setStyle(TableStyle([
        ("VALIGN",       (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING",  (0,0), (-1,-1), 4),
        ("RIGHTPADDING", (0,0), (-1,-1), 4),
        ("TOPPADDING",   (0,0), (-1,-1), 0),
        ("BOTTOMPADDING",(0,0), (-1,-1), 0),
    ]))
    rows.append([badge_cell])

    # Meta (PMID, DOI, link)
    for ml in meta_lines:
        ml_clean = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", html_lib.escape(ml))
        if "pubmed" in ml.lower():
            rows.append([Paragraph(ml_clean, styles["Link"])])
        else:
            rows.append([Paragraph(ml_clean, styles["Meta"])])

    # Bullets abstract
    for bl in bullet_lines:
        bl_clean = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", html_lib.escape(bl))
        # Evidenzia etichetta
        bl_clean = re.sub(
            r'^(Tipo di studio/popolazione|Risultato chiave|Rilevanza clinica'
            r'|Metodi \(full text OA\)|Risultato principale \(full text OA\)'
            r'|Dettaglio utile \(full text OA\)|Messaggio clinico \(full text OA\))(:)',
            r'<b><font color="#1b4f72">\1\2</font></b>',
            bl_clean
        )
        rows.append([Paragraph(f"• {bl_clean}", styles["BulletItem"])])

    card = Table([[r[0]] for r in rows], colWidths=["100%"])
    card.setStyle(TableStyle([
        ("BACKGROUND",   (0, 0), (-1, -1), colors.white),
        ("BOX",          (0, 0), (-1, -1), 0.8, C_BORDER),
        ("LINEAFTER",    (0, 0), (0, -1),  3,   C_ACCENT),
        ("LEFTPADDING",  (0, 0), (-1, -1), 14),
        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
        ("TOPPADDING",   (0, 0), (0,  0),  10),
        ("BOTTOMPADDING",(0,-1), (-1,-1),  10),
        ("TOPPADDING",   (0, 1), (-1, -1), 3),
        ("BOTTOMPADDING",(0, 0), (-1, -2), 3),
        ("ROUNDEDCORNERS", [6]),
    ]))
    return card


def _header_footer(canvas_obj, doc):
    """Intestazione e piè di pagina su ogni pagina."""
    w, h = A4
    canvas_obj.saveState()

    # Header strip
    canvas_obj.setFillColor(C_ACCENT)
    canvas_obj.rect(0, h - 28, w, 28, fill=1, stroke=0)
    canvas_obj.setFont("Helvetica-Bold", 9)
    canvas_obj.setFillColor(C_WHITE)
    canvas_obj.drawString(1.8*cm, h - 18, "LMC Newsletter  ·  Ponatinib / Asciminib")
    canvas_obj.setFont("Helvetica", 9)
    canvas_obj.drawRightString(w - 1.8*cm, h - 18,
        datetime.utcnow().strftime("%d/%m/%Y"))

    # Footer
    canvas_obj.setFillColor(C_BORDER)
    canvas_obj.rect(0, 0, w, 20, fill=1, stroke=0)
    canvas_obj.setFont("Helvetica", 8)
    canvas_obj.setFillColor(C_INK_LIGHT)
    canvas_obj.drawCentredString(w/2, 6,
        f"Pagina {doc.page}  ·  Uso interno riservato  ·  Generata automaticamente da PubMed/EuropePMC")

    canvas_obj.restoreState()


def markdown_to_simple_pdf(md_text: str, out_path: Path):
    styles = build_pdf_styles()

    doc = SimpleDocTemplate(
        str(out_path),
        pagesize=A4,
        leftMargin=1.8*cm, rightMargin=1.8*cm,
        topMargin=2.2*cm,  bottomMargin=1.8*cm,
        title="LMC Newsletter – Ponatinib / Asciminib",
        author="Sistema automatico",
    )

    story = []
    lines = md_text.splitlines()

    # Stato del parsing
    in_bullet_block = False
    bullet_buffer   = []
    current_section = None   # "messaggi" | "azioni" | "radar" | None
    in_card         = False
    card_title      = ""
    card_journal    = ""
    card_stype      = ""
    card_score      = 0
    card_meta       = []
    card_bullets    = []
    in_top5         = False
    article_rank    = 0

    SECTION_HEADERS = {
        "messaggi chiave": ("messaggi",  "📌 Messaggi chiave (da leggere in 60 secondi)", colors.HexColor("#eaf3fb")),
        "cosa cambia":     ("azioni",    "⚡ Cosa cambia per noi (azioni pratiche)",       colors.HexColor("#fef9e7")),
        "safety radar":    ("radar",     "⚠️ Safety radar (segnali da monitorare)",        colors.HexColor("#fdf2f0")),
    }

    def flush_bullets():
        nonlocal bullet_buffer, current_section, in_bullet_block
        if bullet_buffer and current_section:
            for key, (sid, title, col) in SECTION_HEADERS.items():
                if current_section == sid:
                    story.append(make_callout_table(title, bullet_buffer, col))
                    story.append(Spacer(1, 8))
                    break
        elif bullet_buffer:
            for b in bullet_buffer:
                b_clean = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", html_lib.escape(b))
                story.append(Paragraph(f"• {b_clean}", styles["BulletItem"]))
        bullet_buffer = []
        in_bullet_block = False
        current_section = None

    def flush_card():
        nonlocal in_card, card_title, card_journal, card_stype, card_score, card_meta, card_bullets
        if in_card and card_title:
            card = make_article_card(
                card_title, card_journal, card_stype, card_score,
                card_meta, card_bullets, styles
            )
            story.append(KeepTogether([card, Spacer(1, 10)]))
        in_card = False
        card_title = card_journal = card_stype = ""
        card_score = 0
        card_meta = []
        card_bullets = []

    for raw in lines:
        line = raw.rstrip()

        # Riga vuota
        if not line.strip():
            if in_bullet_block:
                pass  # continua ad accumulare
            continue

        # Separatore
        if line.strip() == "---":
            flush_bullets()
            flush_card()
            story.append(HRFlowable(width="100%", thickness=0.5,
                                    color=C_BORDER, spaceAfter=6, spaceBefore=6))
            continue

        # H1
        if line.startswith("# "):
            flush_bullets()
            flush_card()
            title_text = html_lib.escape(line[2:])
            # Titolo come banner colorato
            title_para = Paragraph(title_text, styles["Title"])
            sub_para   = Paragraph(
                f"Generata il {datetime.utcnow().strftime('%d/%m/%Y')}",
                styles["Subtitle"]
            )
            banner = Table(
                [[title_para], [sub_para]],
                colWidths=["100%"],
            )
            banner.setStyle(TableStyle([
                ("BACKGROUND",   (0,0),(-1,-1), C_ACCENT),
                ("LEFTPADDING",  (0,0),(-1,-1), 20),
                ("RIGHTPADDING", (0,0),(-1,-1), 20),
                ("TOPPADDING",   (0,0),(0,0),   18),
                ("BOTTOMPADDING",(0,-1),(-1,-1),18),
                ("TOPPADDING",   (0,1),(-1,-1),  4),
            ]))
            story.append(banner)
            story.append(Spacer(1, 18))
            continue

        # H2
        if line.startswith("## "):
            flush_bullets()
            flush_card()
            txt = line[3:]
            story.append(Paragraph(html_lib.escape(txt), styles["H2"]))
            story.append(HRFlowable(width="100%", thickness=1.5,
                                    color=C_ACCENT, spaceAfter=6, spaceBefore=0))
            in_top5 = "top 5" in txt.lower()
            continue

        # H3
        if line.startswith("### "):
            flush_bullets()
            flush_card()
            txt = html_lib.escape(line[4:])
            # Se inizia con numero → titolo articolo top5
            if in_top5 and re.match(r'^\d+\)', txt):
                in_card = True
                card_title  = line[4:]
                card_meta   = []
                card_bullets= []
                card_journal = ""
                card_stype  = ""
                card_score  = 0
            elif not in_top5:
                # Bucket header o articolo nelle sezioni per domanda clinica
                in_card = True
                card_title  = line[4:]
                card_meta   = []
                card_bullets= []
                card_journal = ""
                card_stype  = ""
                card_score  = 0
            else:
                story.append(Paragraph(txt, styles["H3"]))
            continue

        # Bullet item
        if line.startswith("- "):
            item = line[2:]

            # Riconosci sezione callout
            item_lower = item.lower()
            for key, (sid, _, _) in SECTION_HEADERS.items():
                if key in item_lower and item.startswith("**"):
                    flush_bullets()
                    current_section = sid
                    in_bullet_block = True
                    break

            # Badge rivista/tipo → meta della card
            if item.startswith("**Rivista**:"):
                m = re.search(r"\*\*Rivista\*\*:\s*(.*?)\s*\|\s*\*\*Tipo\*\*:\s*(.*)$", item)
                if m and in_card:
                    card_journal = m.group(1).strip()
                    card_stype   = m.group(2).strip()
                continue

            # Impact score
            if item.startswith("**Impact score**:") or item.startswith("Impact score:"):
                m = re.search(r"(\d+)", item)
                if m and in_card:
                    card_score = int(m.group(1))
                continue

            # Focus clinico, PMID, link, OA → meta
            if any(item.startswith(p) for p in ["**Focus clinico**:", "PMID:", "Link PubMed:", "Open Access"]):
                if in_card:
                    card_meta.append(item)
                continue

            # Mini-riassunto / full text → bullets card
            if item.startswith("**Mini-riassunto") or item.startswith("**Riassunto aggiuntivo"):
                continue  # etichetta sezione, saltiamo

            # Righe abstract bullet
            bullet_labels = [
                "Tipo di studio", "Risultato chiave", "Rilevanza clinica",
                "Metodi (full text", "Risultato principale (full text",
                "Dettaglio utile (full text", "Messaggio clinico (full text"
            ]
            if in_card and any(item.startswith(l) for l in bullet_labels):
                card_bullets.append(item)
                continue

            # Bullet generico
            if in_bullet_block:
                bullet_buffer.append(item)
            elif in_card:
                card_bullets.append(item)
            else:
                b_clean = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", html_lib.escape(item))
                story.append(Paragraph(f"• {b_clean}", styles["BulletItem"]))
            continue

        # Paragrafo normale
        flush_bullets()
        txt = html_lib.escape(line)
        txt = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", txt)
        txt = re.sub(r"_(.*?)_", r"<i>\1</i>", txt)
        story.append(Paragraph(txt, styles["Body"]))

    # Flush finale
    flush_bullets()
    flush_card()

    doc.build(story, onFirstPage=_header_footer, onLaterPages=_header_footer)


# =========================
# MAIN
# =========================
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
    if not isinstance(html_text, str) or not html_text.strip():
        html_text = "<html><body><pre>" + html_lib.escape(md) + "</pre></body></html>"
    OUT_HTML.write_text(html_text, encoding="utf-8")

    markdown_to_simple_pdf(md, OUT_PDF)

    mark_seen(con, new_articles)

    print(f"Creato file Markdown: {OUT_MD}")
    print(f"Creato file HTML: {OUT_HTML}")
    print(f"Creato file PDF: {OUT_PDF}")


if __name__ == "__main__":
    main(days=15)
