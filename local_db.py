#!/usr/bin/env python3
"""
Local flat-file database parser for PTM data.
Supports:
  - UniProt Swiss-Prot flat file (uniprot_sprot_human.dat)
    * now also parses /evidence=ECO:... so curated experimental sites are not
      silently downgraded
  - User-supplied extensions TSV (data/ptm_extensions.tsv)
"""
import os
import re

PTM_TYPES = {"MOD_RES", "CARBOHYD", "LIPID"}


def classify_ft_desc(ftype, desc):
    desc = desc.lower()
    if ftype == "MOD_RES":
        if "phosphoserine" in desc or "phosphothreonine" in desc or "phosphotyrosine" in desc:
            return "phosphorylation"
        if "ubiquitylated" in desc or "ubiquitin" in desc:
            return "ubiquitination"
        # acetyl-lysine is NOT ubiquitination; out of scope here
        if "acetyl" in desc:
            return None
        if "glycosyl" in desc or "carbohydrate" in desc:
            return "glycosylation"
    elif ftype == "CARBOHYD":
        return "glycosylation"
    elif ftype == "LIPID":
        if "palm" in desc:
            return "palmitoylation"
    return None


def parse_uniprot_flat_file(path=None):
    """Parse UniProt flat file.

    Returns {accession: {position: {"ptm_type": ..., "evidence": set_of_ECO}}}
    """
    if path is None:
        from download_data import FLAT_FILE
        path = FLAT_FILE

    if not os.path.exists(path):
        return {}

    ptm_map = {}
    current_accs = []
    feature = None  # (ftype, start, end)
    desc_parts = []
    ev_codes = set()

    def flush_feature():
        nonlocal feature, desc_parts, ev_codes
        if feature is None:
            return
        ftype, start, end = feature
        desc = " ".join(desc_parts)
        ptm = classify_ft_desc(ftype, desc)
        if ptm and start:
            for acc in current_accs:
                bucket = ptm_map.setdefault(acc, {})
                for pos in range(start, end + 1):
                    bucket[pos] = {"ptm_type": ptm, "evidence": set(ev_codes)}
        feature = None
        desc_parts = []
        ev_codes = set()

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("AC"):
                for part in line.split()[1:]:
                    part = part.rstrip(";")
                    if part:
                        current_accs.append(part)
            elif line.startswith("//"):
                flush_feature()
                current_accs = []
            elif line.startswith("FT"):
                parts = line[2:].split(None, 3)
                if parts and parts[0] in PTM_TYPES:
                    flush_feature()
                    try:
                        start = int(parts[1])
                    except (ValueError, IndexError):
                        start = None
                    try:
                        end = int(parts[2])
                    except (ValueError, IndexError):
                        end = start
                    feature = (parts[0], start, end)
                    desc_parts = []
                    if len(parts) > 3:
                        desc_parts.append(parts[3])
                elif feature is not None:
                    cont = line[2:].strip()
                    if cont.startswith("/note="):
                        cont = cont[len("/note="):].strip().strip('"')
                    elif cont.startswith("/evidence="):
                        m = re.search(r"ECO:\d{7}", cont)
                        if m:
                            ev_codes.add(m.group(0))
                        continue
                    desc_parts.append(cont)

    flush_feature()
    return ptm_map


def parse_extensions_tsv(path=None):
    """Parse user-supplied TSV: accession\\tposition\\tptm_type\\tsource"""
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "data", "ptm_extensions.tsv")
    if not os.path.exists(path):
        return {}
    result = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            acc, pos_str, ptm = parts[0], parts[1], parts[2]
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            result.setdefault(acc, {})[pos] = ptm
    return result
