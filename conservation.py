"""
Cross-species conservation via global pairwise alignment (Needleman-Wunsch).

Rationale: functional PTM sites (e.g. GRK/PKA phospho sites, the C-tail
palmitoylation Cys, N-glycosylation sequons) are usually conserved across
orthologs. We align each ortholog to the human reference and read off the
residue at the aligned column, which is far more reliable than sliding-window
matches (which spuriously align fast-evolving tails).

Limitation: orthologs are looked up by exact gene name in UniProt, so only
taxa sharing the gene name are found (here mostly mammals). This gives high
resolution for divergent sites but saturates at ~1.0 for conserved ones.
"""
import requests

MATCH = 1
MISMATCH = -1
GAP = -2

TAXA = {
    "mouse": "10090",
    "rat": "10116",
    "dog": "9615",
    "cow": "9913",
    "macaque": "9544",
    "zebrafish": "7955",
    "fruit_fly": "7227",
    "worm": "6239",
}


def fetch_orthologs(gene):
    sequences = {}
    for org, taxon in TAXA.items():
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": f"gene:{gene} AND organism_id:{taxon} AND reviewed:true",
            "fields": "accession,sequence",
            "format": "json",
            "size": 1,
        }
        try:
            r = requests.get(url, params=params, timeout=20)
            if r.status_code != 200:
                continue
            data = r.json()
            results = data.get("results", [])
            if results:
                seq = results[0].get("sequence", {}).get("value", "")
                if seq and len(seq) > 50:
                    sequences[org] = seq
        except Exception:
            continue
    return sequences


def global_align(a, b):
    """Needleman-Wunsch global alignment.

    Returns a dict mapping 0-based index in `a` (human) to the aligned 0-based
    index in `b` (ortholog), or None when the human residue aligns to a gap.
    """
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        score[i][0] = GAP * i
    for j in range(m + 1):
        score[0][j] = GAP * j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (MATCH if a[i - 1] == b[j - 1] else MISMATCH)
            score[i][j] = max(diag, score[i - 1][j] + GAP, score[i][j - 1] + GAP)

    # traceback
    mapping = {}
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i - 1][j - 1] + (MATCH if a[i - 1] == b[j - 1] else MISMATCH):
            mapping[i - 1] = j - 1
            i, j = i - 1, j - 1
        elif i > 0 and score[i][j] == score[i - 1][j] + GAP:
            mapping[i - 1] = None
            i -= 1
        else:
            j -= 1
    return mapping


def compute_conservation(predictions, sequence, record):
    gene = record.get("gene", "")
    if not gene:
        return predictions

    orthologs = fetch_orthologs(gene)
    if not orthologs:
        return predictions

    maps = []
    for _org, seq in orthologs.items():
        if len(seq) >= 0.5 * len(sequence):
            maps.append((seq, global_align(sequence, seq)))
    if not maps:
        return predictions

    for p in predictions:
        pos = p["position"]  # 1-based
        idx = pos - 1
        total, conserved = 0, 0
        for _seq, mapping in maps:
            j = mapping.get(idx)
            if j is None:
                continue  # human residue falls in a gap -> skip
            total += 1
            if _seq[j] == p["residue"]:
                conserved += 1
        p["conservation"] = round(conserved / total, 3) if total else 0.5

    return predictions
