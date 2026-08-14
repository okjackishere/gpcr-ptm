#!/usr/bin/env python3
"""
PTM data sources: iPTMnet live API + dbPTM bulk flat file.
"""
import gzip
import io
import os
import re
import sqlite3
import tarfile

import requests

IPTMNET_API = "https://research.bioinformatics.udel.edu/iptmnet/api"

# dbPTM 2025 (Chung et al., NAR 2025) 已从 CUHK 迁至 NYCU Biomics Lab:
# 旧地址 awi.cuhk.edu.cn 整站 403 不可用; 新站按 PTM 类型分文件,
# 每个文件是 tar.gz 归档, 内含一个无表头的 5 列 TSV
# (蛋白名 / UniProt AC / 位置 / PTM类型 / PMID / 序列片段)。
# 本工具只取范围内四类。
DBPTM_BASE = "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment"
DBPTM_FILES = (
    "Phosphorylation.gz",
    "N-linked Glycosylation.gz",
    "S-palmitoylation.gz",
    "Ubiquitination.gz",
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
DBPTM_FILE = os.path.join(DATA_DIR, "dbptm.txt")  # 旧版全量文本, 仅用于一次性迁移
DBPTM_DB = os.path.join(DATA_DIR, "dbptm.db")     # SQLite 库, 查询用的单一数据源

UNIPROT_AC_RE = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$')
POSITION_RE = re.compile(r'^([A-Z])(\d+)$')  # e.g. S345, N6
PMID_LIST_RE = re.compile(r'\d{4,9}(;\d{4,9})*')  # 单个或分号分隔的多个 PMID (4位以上含早期ID)


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


def _iter_dbptm_rows(lines):
    """行解析核心: 逐行产出 (accession, pos, ptm_type, residue, pmid)。

    Format is flexible: header row with recognizable column names, otherwise
    heuristic detection (UniProt AC pattern, integer/S123 position, PTM keyword).
    """
    lines = list(lines)
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
                elif PMID_LIST_RE.fullmatch(p):
                    # 分号分隔的多 PMID (dbPTM 2025 磷酸化格式, 如 17925438;8094205)
                    if pmid is None:
                        pmid = p
                elif p.isdigit():
                    # 纯数字: 位置在前(残基数<=5位), 位置已定后的数字视为单 PMID
                    if pos is None and (residue is None or len(p) <= 5):
                        pos = int(p)
                    elif pmid is None and 4 <= len(p) <= 9:
                        pmid = p

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

        yield acc, pos, ptm_type, residue or "", pmid or ""


def parse_dbptm_tsv(path=None):
    """Parse dbPTM 文本 -> {accession: {pos: info}}。仅用于旧 txt 的一次性迁移。"""
    if path is None:
        path = DBPTM_FILE
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip() and not ln.startswith("#")]
    result = {}
    for acc, pos, ptm_type, residue, pmid in _iter_dbptm_rows(lines):
        result.setdefault(acc, {})[pos] = {
            "ptm_type": ptm_type,
            "residue": residue,
            "source": "dbPTM",
            "pmid": pmid,
        }
    return result


# ---------------------------------------------------------------------------
# dbPTM SQLite 库 (标准库 sqlite3, 按 accession 索引, 查询毫秒级)
# ---------------------------------------------------------------------------
def _build_dbptm_db(rows, db_path=None):
    """把 (acc, pos, ptm_type, residue, pmid) 行流写入 SQLite, 返回位点总数。"""
    if db_path is None:
        db_path = DBPTM_DB
    if os.path.exists(db_path):
        os.remove(db_path)
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("""CREATE TABLE sites (
            accession TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ptm_type TEXT NOT NULL,
            residue TEXT DEFAULT '',
            pmid TEXT DEFAULT '',
            PRIMARY KEY (accession, pos)
        )""")
        conn.executemany("INSERT OR REPLACE INTO sites VALUES (?,?,?,?,?)", rows)
        conn.commit()
        return conn.execute("SELECT COUNT(*) FROM sites").fetchone()[0]
    finally:
        conn.close()


def _migrate_txt_to_db():
    """旧版全量 dbptm.txt 一次性转为 SQLite; 成功后删除 txt (保持单一数据源)。"""
    print("[*] 检测到旧版 data/dbptm.txt, 正在一次性迁移为 SQLite (约 30-60 秒)...")
    with open(DBPTM_FILE, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip() and not ln.startswith("#")]
        n = _build_dbptm_db(_iter_dbptm_rows(lines))
    os.remove(DBPTM_FILE)
    print(f"[*] 迁移完成: {DBPTM_DB} ({n} 位点); 已删除旧文本 {DBPTM_FILE}")


def query_dbptm(accession, db_path=None):
    """按 accession 查询 dbPTM 位点, 返回 {pos: {ptm_type, residue, source, pmid}}。

    只有旧 txt 没有 db 时自动迁移一次; 两者皆缺返回空 (正常降级, 不报错)。
    """
    if db_path is None:
        db_path = DBPTM_DB
    if not os.path.exists(db_path):
        if not os.path.exists(DBPTM_FILE):
            return {}
        try:
            _migrate_txt_to_db()
        except Exception as e:
            print(f"[!] dbPTM txt->db 迁移失败, 本轮跳过 dbPTM 数据: {e}")
            return {}
    if not os.path.exists(db_path):
        return {}
    conn = sqlite3.connect(db_path)
    try:
        rows = conn.execute(
            "SELECT pos, ptm_type, residue, pmid FROM sites WHERE accession = ?",
            (accession,)).fetchall()
    finally:
        conn.close()
    return {pos: {"ptm_type": t, "residue": r, "source": "dbPTM", "pmid": p}
            for pos, t, r, p in rows}


def _extract_text(raw):
    """dbPTM 数据是 tar.gz 归档(内含无表头 TSV); 兼容纯 gzip 文本与纯文本。"""
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    # tar 魔数 "ustar" 固定在 257 偏移处; 非归档则按纯文本处理
    if len(raw) > 262 and raw[257:262] == b"ustar":
        with tarfile.open(fileobj=io.BytesIO(raw)) as tf:
            for m in tf.getmembers():
                if m.isfile():
                    return tf.extractfile(m).read().decode("utf-8", errors="replace")
    return raw.decode("utf-8", errors="replace")


def download_dbptm(url=None, db_path=None):
    """下载 dbPTM 实验数据并直接解析入 SQLite 库 (data/dbptm.db)。

    默认从 dbPTM 2025 新站 (NYCU) 下载范围内四类文件;
    显式给 url (参数或环境变量 DBPTM_URL) 时按单文件下载, 兼容自定义镜像。
    """
    if url is None:
        url = os.environ.get("DBPTM_URL")
    if db_path is None:
        db_path = DBPTM_DB
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    urls = [url] if url else [f"{DBPTM_BASE}/{name}" for name in DBPTM_FILES]

    def rows():
        for u in urls:
            label = u.rsplit("/", 1)[-1]
            print(f"[*] Downloading dbPTM data: {label} ...")
            r = requests.get(u, timeout=300, stream=True)
            r.raise_for_status()
            raw = b"".join(r.iter_content(chunk_size=1024 * 1024))
            text = _extract_text(raw)
            if not text.strip():
                raise ValueError(f"dbPTM 下载内容为空: {u}")
            yield from _iter_dbptm_rows(text.splitlines())

    n = _build_dbptm_db(rows(), db_path)
    print(f"[*] dbPTM ready: {db_path} ({n} 位点)")


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
