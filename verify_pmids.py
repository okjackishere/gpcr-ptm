#!/usr/bin/env python3
"""
逐条核验 PMID: 抓取 PubMed 标题+摘要, 自动判断该文献是否支持"该残基 + 该PTM类型"。

判定 (verdict):
  直接               摘要/标题点名该残基 且 含该PTM类型关键词
  间接(提及位点)     点名该残基, 但未命中PTM关键词
  间接(未点名位点)   该蛋白该PTM相关(如GRK磷酸化β2AR), 但未点名该残基
  不支持(位点不符)   点名了同位置的其他残基(如 Tyr364 而非 Ser364)
  不支持(非该蛋白)   摘要明确是其他蛋白(如加压素受体)

摘要按 PMID 缓存在 data/pmid_abstract_cache.json, 避免重复请求。
"""
import json
import os
import re
import time

import requests

CACHE_FILE = os.path.join(os.path.dirname(__file__), "data", "pmid_abstract_cache.json")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RATE_LIMIT = 0.35

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}

FULL = {
    "A": "alanine", "R": "arginine", "N": "asparagine", "D": "aspartic acid",
    "C": "cysteine", "Q": "glutamine", "E": "glutamic acid", "G": "glycine",
    "H": "histidine", "I": "isoleucine", "L": "leucine", "K": "lysine",
    "M": "methionine", "F": "phenylalanine", "P": "proline", "S": "serine",
    "T": "threonine", "W": "tryptophan", "Y": "tyrosine", "V": "valine",
}

PTM_KEYWORDS = {
    "phosphorylation": ["phosphorylat", "phosphoserine", "phosphothreonine",
                        "phosphotyrosine", "phospho", "pS(", "pT(", "pY(",
                        "kinase", "grk", "pka", "capk", "desensitiz"],
    "glycosylation": ["glycosyl", "glycan", "n-glyc", "glycosylat"],
    "palmitoylation": ["palmitoyl", "palmitate", "acyl"],
    "ubiquitination": ["ubiquitin", "ubiquitylat"],
}

OTHER_PROTEINS = re.compile(
    r"\b(avp|vasopressin|muscarinic|rhodopsin|dopamine|serotonin|histamine|"
    r"opioid|glucagon|mGlu|GABA|metabotropic|olfactory)\b", re.I)


def _pubmed_url(pmid):
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


STOPWORDS = {"receptor", "protein", "factor", "subunit", "type", "family",
             "member", "alpha", "beta", "gamma", "delta", "human", "mouse",
             "receptors", "protein-coupled"}


def protein_keywords_from_entry(entry):
    """从 UniProt entry 提取该蛋白的检索关键词 (基因名 + 蛋白名特征词)。

    刻意过滤 receptor/protein 等通用词, 避免把任何 GPCR 论文都判定为"该蛋白",
    否则会掩盖"非该蛋白"的误归属 (如加压素受体论文挂在 β2-AR 位点上)。
    """
    kws = []
    for g in entry.get("genes", []):
        gn = (g.get("geneName") or {}).get("value") or ""
        if gn:
            kws.append(gn)
    pd = entry.get("proteinDescription", {})
    name = (pd.get("recommendedName") or {}).get("fullName") or {}
    val = name.get("value", "") if isinstance(name, dict) else str(name)
    if val:
        for t in re.split(r"[\s\-/]+", val.lower()):
            if len(t) > 3 and t not in STOPWORDS:
                kws.append(t)
    kws = sorted({k.lower() for k in kws if k})
    return kws or ["receptor"]


# ---------------------------------------------------------------------------
# 摘要抓取与缓存
# ---------------------------------------------------------------------------
def load_cache():
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE) as f:
                return json.load(f)
        except Exception:
            return {}
    return {}


def save_cache(cache):
    os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
    with open(CACHE_FILE, "w") as f:
        json.dump(cache, f, ensure_ascii=False)


def _parse_efetch_xml(xml_text):
    """efetch (retmode=xml) -> {title, abstract}. 若解析失败返回 None。"""
    try:
        from xml.etree import ElementTree as ET
        root = ET.fromstring(xml_text)
        art = root.find(".//PubmedArticle/MedlineCitation/Article")
        if art is None:
            return None
        title = art.findtext("ArticleTitle") or ""
        abstract = " ".join(t.text or "" for t in art.findall(".//Abstract/AbstractText"))
        return {"title": title, "abstract": abstract}
    except Exception:
        return None


def _parse_efetch_text(txt):
    """NCBI efetch (rettype=abstract, retmode=text) -> {title, abstract}."""
    lines = [ln.strip() for ln in txt.splitlines()]
    title, abstract = "", ""
    i = 0
    # skip citation header lines (e.g. "1. J Biol Chem. ..." or empty)
    while i < len(lines):
        if not lines[i]:
            i += 1
            continue
        if re.match(r"^\d+\.\s", lines[i]):
            i += 1
            continue
        if lines[i].startswith("Comment in"):
            i += 1
            continue
        break
    body = lines[i:]
    # title = first non-empty, non-affiliation line
    for ln in body:
        if not ln:
            continue
        if ln.startswith("Author information:") or ln.startswith("(") and "Department" in ln:
            continue
        title = ln
        break
    # abstract = between title and "Author information:"
    j = body.index(title) + 1 if title in body else 0
    parts = []
    for ln in body[j:]:
        if ln.startswith("Author information:"):
            break
        if ln.startswith("DOI:") or ln.startswith("PMCID:") or ln.startswith("PMID:"):
            break
        parts.append(ln)
    abstract = " ".join(parts).strip()
    return {"title": title, "abstract": abstract}


def fetch_abstracts(pmids, progress=None):
    cache = load_cache()
    missing = [p for p in pmids if p not in cache]
    total = len(missing)
    for i, pid in enumerate(missing, 1):
        if progress:
            progress(f"核验文献 {i}/{total}: PMID {pid} (首次抓取摘要, 约15-30秒)")
        info = {"title": "", "abstract": ""}
        # 1) XML (reliable title+abstract parsing)
        try:
            r = requests.get(EFETCH, params={"db": "pubmed", "id": pid,
                                             "retmode": "xml"}, timeout=40)
            if r.status_code == 200:
                info = _parse_efetch_xml(r.text) or info
        except Exception:
            pass
        # 2) fallback to plain-text abstract
        if not info["abstract"] and not info["title"]:
            try:
                r = requests.get(EFETCH, params={"db": "pubmed", "id": pid,
                                                 "rettype": "abstract",
                                                 "retmode": "text"}, timeout=40)
                if r.status_code == 200:
                    info = _parse_efetch_text(r.text) or info
            except Exception:
                pass
        cache[pid] = info
        time.sleep(RATE_LIMIT)
    if missing:
        save_cache(cache)
    return cache


# ---------------------------------------------------------------------------
# 判定逻辑
# ---------------------------------------------------------------------------
def _residue_patterns(res, pos):
    three = AA3.get(res)
    full = FULL.get(res)
    pats = [rf"\b{res}\s*{pos}\b", rf"\b{res}{pos}\b"]
    if three:
        pats += [
            rf"\b{three}\s*[-–—]?\s*{pos}\b",          # Ser-355 / Ser 355
            rf"\b{three}\s*\(\s*{pos}\s*\)",            # Ser(355)
            rf"\b{three}\s*\(\s*[^)]*?\b{pos}\b",       # Ser(345,346)
            rf"\b{three}\d+\s*,\s*{pos}\b",             # Ser355,356
            rf"\b{three}\d+\s+and\s+{pos}\b",           # Ser355 and 356
            rf"\b{three.lower()}\s*[-–—]?\s*{pos}\b",   # ser-355
        ]
    if full:
        pats += [
            rf"\b{full}\s*[-–—]?\s*{pos}\b",            # tyrosine-141
            rf"\b{full}s?\s*[-–—]?\s*{pos}\b",          # serines 396
            rf"\b{full}\s+(?:residues?|residue)\s*[- ]?\s*{pos}\b",  # serine residues 355
        ]
    pats += [
        rf"\bp{res}\s*\(\s*{pos}\s*\)",
        rf"\bp{res}\s*\(\s*[^)]*?\b{pos}\b",            # pS(355,356)
        rf"\bp{res}{pos}\b",
        rf"[A-Z]\d+\s*/\s*{pos}\b",                     # Y132/141
    ]
    return pats


def _excerpt(text, span, width=90, cap=160):
    """围绕匹配位点截取一段可读证据: 从句子边界开始, 总长不超过 cap。"""
    s, e = span
    lo = max(0, s - width)
    if lo > 0:
        cut = max(text.rfind(". ", lo, s), text.rfind("; ", lo, s),
                  text.rfind("\n", lo, s))
        if cut != -1:
            lo = cut + 2
    hi = min(len(text), e + width)
    out = re.sub(r"\s+", " ", text[lo:hi]).strip()
    if len(out) > cap:
        out = out[:cap].rsplit(" ", 1)[0] + "…"
    return out


def verdict_for(text, ptm_type, residue, pos, protein_keywords=None):
    """Return (verdict, evidence_excerpt)."""
    if protein_keywords is None:
        protein_keywords = ["g-protein", "gpcr", "receptor"]
    low = text.lower()
    ptm_kw = PTM_KEYWORDS.get(ptm_type, [])
    ptm_ok = any(k in low for k in ptm_kw)
    protein_ok = any(k in low for k in protein_keywords)

    # 1) target residue mentioned?
    for pat in _residue_patterns(residue, pos):
        m = re.search(pat, text)
        if m:
            ev = _excerpt(text, m.span())
            return ("直接" if ptm_ok else "间接(提及位点)"), ev

    # 2) a different residue at the SAME position? (e.g. Tyr364 for Ser364)
    for other, three in AA3.items():
        if other == residue:
            continue
        full = FULL.get(other)
        pats = [rf"\b{other}\s*{pos}\b", rf"\b{three}\s*[-–—]?\s*{pos}\b"]
        if full:
            pats.append(rf"\b{full}\s*[-–—]?\s*{pos}\b")
        for pat in pats:
            m = re.search(pat, text)
            if m:
                return ("不支持(位点不符)", f"摘要为 {three}{pos}, 而非 {AA3[residue]}{pos}")

    # 3) protein context
    if not protein_ok:
        if OTHER_PROTEINS.search(low):
            return "不支持(非该蛋白)", _excerpt(text, (0, 0), 140)
        return "间接(未点名位点)", _excerpt(text, (0, 0), 140)
    return "间接(未点名位点)", _excerpt(text, (0, 0), 140)


# ---------------------------------------------------------------------------
# 集成入口
# ---------------------------------------------------------------------------
def attach_verification(records, protein_keywords=None, progress=None):
    """为每条含 PMID 的记录附加 pmid_verification 字段。"""
    pmids = sorted({p for r in records for p in r.get("pmids", []) if p})
    if not pmids:
        return records
    abstracts = fetch_abstracts(pmids, progress=progress)
    for r in records:
        plist = r.get("pmids", [])
        if not plist:
            continue
        ver = []
        for p in plist:
            info = abstracts.get(p) or {"title": "", "abstract": ""}
            text = f"{info.get('title', '')}\n{info.get('abstract', '')}".strip()
            if not text:
                verdict, evidence = "无摘要", info.get("title", "")
            else:
                verdict, evidence = verdict_for(text, r["ptm_type"],
                                                r["residue"], r["position"],
                                                protein_keywords)
            ver.append({"pmid": p, "url": _pubmed_url(p),
                        "verdict": verdict, "evidence": evidence})
        r["pmid_verification"] = ver
    return records


def summarize(record):
    """紧凑的核验汇总, 用于控制台: '直接x2 间接x3 不支持x1:8557631'"""
    ver = record.get("pmid_verification", [])
    if not ver:
        return ""
    d = {"直接": 0, "间接": 0, "不支持": 0, "其他": 0}
    bad = []
    for v in ver:
        vd = v["verdict"]
        if vd == "直接":
            d["直接"] += 1
        elif vd.startswith("间接"):
            d["间接"] += 1
        elif vd.startswith("不支持"):
            d["不支持"] += 1
            bad.append(v["pmid"])
        else:
            d["其他"] += 1
    s = f"✓x{d['直接']} △x{d['间接']} ✗x{d['不支持']}"
    if bad:
        s += f" (✗:{','.join(bad)})"
    return s
