#!/usr/bin/env python3
"""
PTM data sources: iPTMnet live API + dbPTM bulk flat file.
"""
import os
import re
import requests

IPTMNET_API = "https://research.bioinformatics.udel.edu/iptmnet/api"

# dbPTM current host (user's network may access it; this sandbox is IP-blocked).
DBPTM_DEFAULT_URL = ("https://awi.cuhk.edu.cn/dbPTM/download/dbPTM_all.txt")

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
DBPTM_FILE = os.path.join(DATA_DIR, "dbptm.txt")

UNIPROT_AC_RE = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$')
POSITION_RE = re.compile(r'^([A-Z])(\d+)$')  # e.g. S345, N6


# ---------------------------------------------------------------------------
# iPTMnet
# ---------------------------------------------------------------------------
def map_iptmnet_type(ptm_type):
    """Map iPTMnet nomenclature onto our PTM types. None if out of scope.

    C-glycosylation (C-mannosylation on Trp) is a distinct PTM from
    N-glycosylation and is deliberately NOT folded into the glycosylation slot.
    """
    if not ptm_type:
        return None
    t = ptm_type.lower()
    if "phospho" in t:
        return "phosphorylation"
    if "ubiquitin" in t:
        return "ubiquitination"
    if t == "n-glycosylation":
        return "glycosylation"
    if "palm" in t:
        return "palmitoylation"
    # O/S-glycosylation, sumoylation, nitrosylation, acetylation,
    # methylation, myristoylation are out of scope for this tool.
    return None


def fetch_iptmnet_sites(accession):
    """Live query iPTMnet for PTM sites of a protein. Returns list of dicts."""
    url = f"{IPTMNET_API}/v1/{accession}/substrate"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return []
    payload = r.json()
    raw = payload.get(accession, []) or payload.get("substrates", [])
    sites = []
    for item in raw:
        ptm_type = map_iptmnet_type(item.get("ptm_type"))
        if not ptm_type:
            continue
        site = item.get("site")
        m = POSITION_RE.match(site) if site else None
        pos = int(m.group(2)) if m else item.get("position")
        residue = (m.group(1) if m else item.get("residue")) or ""
        if not pos:
            continue
        pmids = item.get("pmids", []) or []
        sources = [s.get("name", "") for s in item.get("sources", []) if s.get("name")]
        sites.append({
            "ptm_type": ptm_type,
            "position": pos,
            "residue": residue,
            "score": float(item.get("score", 1)),
            "source": "iPTMnet",
            "evidence_sources": sources,
            "pmids": pmids,
        })
    return sites


# ---------------------------------------------------------------------------
# dbPTM bulk flat file
# ---------------------------------------------------------------------------
def _detect_columns(header):
    """Auto-detect column indices from a dbPTM header row."""
    col = {}
    for i, h in enumerate(header):
        hl = h.strip().lower()
        if "uniprot" in hl or hl in ("ac", "accession", "acc"):
            col["accession"] = i
        elif hl in ("position", "pos", "site"):
            col["position"] = i
        elif hl in ("residue", "aa", "res"):
            col["residue"] = i
        elif "ptm" in hl or "modification" in hl or "type" in hl:
            col["ptm_type"] = i
        elif "pmid" in hl or "pubmed" in hl:
            col["pmid"] = i
    return col


def parse_dbptm_tsv(path=None):
    """Parse dbPTM flat file -> {accession: {pos: {ptm_type, source}}}.

    Format is flexible: header row with recognizable column names, otherwise
    heuristic detection (UniProt AC pattern, integer/S123 position, PTM keyword).
    """
    if path is None:
        path = DBPTM_FILE
    if not os.path.exists(path):
        return {}

    result = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip() and not ln.startswith("#")]

    header = None
    for ln in lines[:8]:
        parts = ln.split("\t")
        if len(parts) >= 3 and any(k in ln.lower() for k in ("uniprot", "position", "ptm", "modification")):
            header = parts
            break

    header_idx = 0
    if header is not None:
        # find index of header line in lines
        for i, ln in enumerate(lines):
            if ln.split("\t") == header:
                header_idx = i
                break

    for ln in lines[header_idx:]:
        parts = ln.split("\t")
        if len(parts) < 2:
            continue

        acc = pos = residue = ptm = pmid = None

        if header is not None:
            ci = _detect_columns(header)
            if "accession" in ci and ci["accession"] < len(parts):
                acc = parts[ci["accession"]].strip()
            if "position" in ci and ci["position"] < len(parts):
                ps = parts[ci["position"]].strip()
                pos = int(ps) if ps.isdigit() else None
            if "residue" in ci and ci["residue"] < len(parts):
                residue = parts[ci["residue"]].strip()
            if "ptm_type" in ci and ci["ptm_type"] < len(parts):
                ptm = parts[ci["ptm_type"]].strip()
            if "pmid" in ci and ci["pmid"] < len(parts):
                pmid = parts[ci["pmid"]].strip()
        else:
            for p in parts:
                p = p.strip()
                if acc is None and UNIPROT_AC_RE.match(p):
                    acc = p
                elif ptm is None and any(k in p.lower() for k in
                                         ("phospho", "glyco", "ubiqui", "palm", "myrist", "sumo")):
                    ptm = p
            # find residue+position token like S355 or separate tokens
            for p in parts:
                p = p.strip()
                if p == acc:
                    continue  # accession token (e.g. P07550) looks like P+7550
                m = POSITION_RE.match(p)
                if m and pos is None:
                    residue, pos = m.group(1), int(m.group(2))
                elif p.isdigit() and pos is None and (residue is None or len(p) <= 5):
                    pos = int(p)

        if not acc or pos is None:
            continue

        ptm_type = None
        if ptm:
            pl = ptm.lower()
            if "phospho" in pl:
                ptm_type = "phosphorylation"
            elif "glyco" in pl:
                ptm_type = "glycosylation"
            elif "ubiqui" in pl:
                ptm_type = "ubiquitination"
            elif "palm" in pl or "myrist" in pl:
                ptm_type = "palmitoylation"
        if not ptm_type:
            continue

        result.setdefault(acc, {})[pos] = {
            "ptm_type": ptm_type,
            "residue": residue or "",
            "source": "dbPTM",
            "pmid": pmid or "",
        }

    return result


def download_dbptm(url=None, path=None):
    if url is None:
        url = os.environ.get("DBPTM_URL", DBPTM_DEFAULT_URL)
    if path is None:
        path = DBPTM_FILE
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"[*] Downloading dbPTM flat file from {url} ...")
    r = requests.get(url, timeout=120, stream=True)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
    print(f"[*] Saved to {path}")


def parse_iptmnet_bulk(path=None):
    """Parse iPTMnet bulk ptm.txt -> {accession: {pos: {ptm_type, source}}}.

    Columns (tab): PTM_TYPE, source, UniProt_AC, protein_name, organism,
    position (e.g. K652), ..., ..., source_entry, PMID.
    """
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "data", "iptmnet_ptm.txt")
    if not os.path.exists(path):
        return {}

    result = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            ptm_raw = parts[0].strip().lower()
            acc = parts[2].strip()
            site = parts[5].strip()
            m = POSITION_RE.match(site)
            if not m or not acc:
                continue
            pos = int(m.group(2))
            residue = m.group(1)

            ptm_type = None
            if "phospho" in ptm_raw:
                ptm_type = "phosphorylation"
            elif "glyco" in ptm_raw:
                ptm_type = "glycosylation"
            elif "ubiqui" in ptm_raw:
                ptm_type = "ubiquitination"
            elif "palm" in ptm_raw or "myrist" in ptm_raw:
                ptm_type = "palmitoylation"
            if not ptm_type:
                continue

            pmid = parts[9].strip() if len(parts) > 9 else ""
            result.setdefault(acc, {})[pos] = {
                "ptm_type": ptm_type,
                "residue": residue,
                "source": "iPTMnet-bulk",
                "pmid": pmid,
            }

    return result
