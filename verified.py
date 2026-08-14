"""
Extract PTM sites with three-tier evidence grading.

  L1 "Verified"  已查证·实验证据  : UniProt 人工实验证据 (ECO:0000269, 带 PMID)
  L2 "Supported" 有支持·未定论    : UniProt 相似性/推断 (ECO:0000255/0305/7744),
                                     iPTMnet / dbPTM 文献或高通量单一检出
  L3 "Predicted" 预测(未被查证)   : 规则预测 (predict.py), 不含数据库记录

合并规则: 每个 (ptm_type, position) 只有一条; 证据等级取最高层, 来源/PMID 合并。
"""
from collections import defaultdict

ECO_EXPERIMENTAL = "ECO:0000269"
ECO_INFERRED = {"ECO:0000255", "ECO:0000250", "ECO:0000305", "ECO:0007744"}

LAYER_SCORE = {"Verified": 1.0, "Supported": 0.75, "Predicted": 0.5}


def classify_feature(feature):
    """Map a UniProt feature onto our PTM types. Out-of-scope -> None."""
    ftype = feature.get("type", "")
    desc = feature.get("description", "") or ""
    desc = desc.lower()

    if ftype == "Modified residue":
        if "phosphoserine" in desc or "phosphothreonine" in desc or "phosphotyrosine" in desc:
            return "phosphorylation"
        if "ubiquitylated" in desc or "ubiquitin" in desc:
            return "ubiquitination"
        # 乙酰化 != 泛素化; 不在本工具范围, 跳过而非错分
        if "acetyl" in desc:
            return None
        if "glycosyl" in desc:
            return "glycosylation"
    elif ftype == "Glycosylation":
        return "glycosylation" if "n-linked" in desc or "glycosylation" in desc else None
    elif ftype == "Lipidation":
        if "palm" in desc:
            return "palmitoylation"
    return None


def get_context(seq, pos, window=5):
    start = max(0, pos - 1 - window)
    end = min(len(seq), pos + window)
    return seq[start:end]


def _eco_codes(evidences):
    codes = set()
    for ev in evidences:
        code = ev.get("evidenceCode")
        if code:
            codes.add(code)
    return codes


def _ev_pmids(evidences):
    pmids = []
    for ev in evidences:
        if ev.get("source") == "PubMed" and ev.get("id"):
            pmids.append(ev["id"])
    return pmids


def _ev_sources(evidences):
    srcs = []
    for ev in evidences:
        s = ev.get("source")
        if s and s not in srcs:
            srcs.append(s)
    return srcs


def _record(ptm_type, pos, residue, region, layer, layer_detail, score, source, pmids, context):
    return {
        "ptm_type": ptm_type,
        "position": pos,
        "residue": residue,
        "region": region,
        "layer": layer,
        "layer_detail": layer_detail,
        "score": score,
        "source": source,
        "pmids": pmids,
        "context": context,
    }


def merge_records(records):
    """Deduplicate by (ptm_type, position), keeping the strongest layer."""
    order = {"Verified": 2, "Supported": 1, "Predicted": 0}
    merged = {}
    for r in records:
        key = (r["ptm_type"], r["position"])
        cur = merged.get(key)
        if cur is None or order[r["layer"]] > order[cur["layer"]]:
            merged[key] = dict(r)
            if cur is not None:
                merged[key]["source"] = _join(cur["source"], r["source"])
                merged[key]["pmids"] = _join(cur["pmids"], r["pmids"])
        else:
            cur["source"] = _join(cur["source"], r["source"])
            cur["pmids"] = _join(cur["pmids"], r["pmids"])
    return list(merged.values())


def _join(a, b):
    if not a:
        return b
    if not b:
        return a
    if isinstance(a, list) and isinstance(b, list):
        seen = []
        for x in list(a) + list(b):
            if x not in seen:
                seen.append(x)
        return seen
    seen = []
    for x in str(a).split(",") + str(b).split(","):
        x = x.strip()
        if x and x not in seen:
            seen.append(x)
    return ", ".join(seen)


def extract_verified(entry, iptmnet_data, topology, progress=None):
    seq = entry["sequence"]["value"]
    acc = entry.get("primaryAccession")
    records = []

    # ------------------------------------------------------------------
    # 1. UniProt live features (authoritative; evidence codes available)
    # ------------------------------------------------------------------
    experimental = set()  # (type, pos) known experimentally -> L1 corroboration
    for f in entry.get("features", []):
        ptm_type = classify_feature(f)
        if not ptm_type:
            continue
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value") or start
        if not start:
            continue
        evs = f.get("evidences", [])
        codes = _eco_codes(evs)
        pmids = _ev_pmids(evs)
        is_exp = ECO_EXPERIMENTAL in codes
        for pos in range(start, end + 1):
            if pos < 1 or pos > len(seq):
                continue
            if is_exp:
                experimental.add((ptm_type, pos))
                records.append(_record(
                    ptm_type, pos, seq[pos - 1], topology.get(pos, "unknown"),
                    "Verified",
                    f"UniProt实验证据; PMID:{','.join(pmids) or 'NA'}",
                    1.0, "UniProt", pmids, get_context(seq, pos)))
            else:
                records.append(_record(
                    ptm_type, pos, seq[pos - 1], topology.get(pos, "unknown"),
                    "Supported",
                    f"UniProt推断/相似性({','.join(sorted(codes)) or 'NA'})",
                    0.75, "UniProt", pmids, get_context(seq, pos)))

    # ------------------------------------------------------------------
    # 2. Offline UniProt Swiss-Prot flat file (with /evidence=ECO parsing)
    # ------------------------------------------------------------------
    try:
        from local_db import parse_uniprot_flat_file
        flat = parse_uniprot_flat_file()
        if acc and acc in flat:
            for pos, info in flat[acc].items():
                if not (1 <= pos <= len(seq)):
                    continue
                ptm_type = info.get("ptm_type")
                if not ptm_type:
                    continue
                is_exp = ECO_EXPERIMENTAL in (info.get("evidence") or set())
                layer = "Verified" if is_exp else "Supported"
                records.append(_record(
                    ptm_type, pos, seq[pos - 1], topology.get(pos, "unknown"),
                    layer,
                    "Swiss-Prot实验证据" if is_exp else "Swiss-Prot推断",
                    1.0 if is_exp else 0.7, "UniProt-flat", [], get_context(seq, pos)))
    except Exception:
        pass

    # ------------------------------------------------------------------
    # 3. iPTMnet live API (score 1-4 = 其自身置信度; PMID = 文献支持)
    # ------------------------------------------------------------------
    for item in iptmnet_data:
        pos = item.get("position")
        ptm_type = item.get("ptm_type")
        if not ptm_type or not pos or not (1 <= pos <= len(seq)):
            continue
        pmids = item.get("pmids", []) or []
        n_src = len(item.get("evidence_sources", []) or [])
        iptm_score = item.get("score", 1)
        if pmids:
            detail = f"iPTMnet: {len(pmids)}篇文献/{(n_src)}库"
            score = 0.8 if iptm_score >= 3 else 0.7
        else:
            detail = "iPTMnet: 高通量单一检出(无PMID), 需实验确认"
            score = 0.6
        records.append(_record(
            ptm_type, pos, seq[pos - 1], topology.get(pos, "unknown"),
            "Supported", detail, score, "iPTMnet", pmids, get_context(seq, pos)))

    # ------------------------------------------------------------------
    # 4. dbPTM (SQLite 按蛋白查询) / iPTMnet-bulk (离线 bulk 文件)
    # ------------------------------------------------------------------
    def _append_bulk(sites, label):
        for pos, info in sites.items():
            if not (1 <= pos <= len(seq)):
                continue
            # dbPTM 2025 的 PMID 列是分号分隔的多个值, 拆开逐条核验
            pmids = [x for x in (info.get("pmid") or "").split(";") if x]
            records.append(_record(
                info["ptm_type"], pos, seq[pos - 1], topology.get(pos, "unknown"),
                "Supported",
                f"{label}数据库记录" + (f"; PMID:{','.join(pmids)}" if pmids else ""),
                0.7 if pmids else 0.6,
                label, pmids, get_context(seq, pos)))

    try:
        from ptm_sources import query_dbptm
        _append_bulk(query_dbptm(acc), "dbPTM")
    except Exception:
        pass

    try:
        from ptm_sources import parse_iptmnet_bulk
        try:
            data = parse_iptmnet_bulk()
        except Exception:
            data = {}
        if acc in data:
            _append_bulk(data[acc], "iPTMnet-bulk")
    except Exception:
        pass

    # iPTMnet 位点若与 UniProt 实验证据吻合, 合并后自动升级为 Verified
    records = merge_records(records)
    records.sort(key=lambda r: (-LAYER_SCORE[r["layer"]], -r["score"]))

    # 文献逐条核验: 为每条带 PMID 的记录附加 pmid_verification
    try:
        from verify_pmids import attach_verification, protein_keywords_from_entry
        records = attach_verification(records, protein_keywords=protein_keywords_from_entry(entry),
                                      progress=progress)
    except Exception:
        pass

    return records
