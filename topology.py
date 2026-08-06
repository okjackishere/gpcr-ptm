"""
Assigns each residue (1-based) a region label:
  N-term, TM1..TM7, ECL1..3, ICL1..3, C-tail, unknown.

Priority order:
  1. UniProt "Topological domain" descriptions (extracellular/cytoplasmic) --
     the most reliable, used verbatim (first = N-term, last = C-tail).
  2. UniProt "Transmembrane" features.
  3. Gap inference between consecutive TMs (ICL after odd TM, ECL after even TM).
  4. If the entry has NO transmembrane annotation at all, a Kyte-Doolittle
     hydrophobicity-window fallback assigns N-term/C-tail only (loops are left
     unknown, so PTM prediction stays conservative rather than guessing).
"""

KD_WINDOW = 19
KD_THRESHOLD = 1.2
MIN_TM_RUN = 15

# Kyte-Doolittle hydrophobicity (fractional scale)
HYDRO = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2,
}


def _desc_of(feature):
    desc = feature.get("description", "")
    if isinstance(desc, dict):
        desc = desc.get("value", "")
    return str(desc).lower()


def _kd_window_topology(seq):
    """Kyte-Doolittle fallback: label N-term / C-tail only, rest unknown."""
    n = len(seq)
    kd = [HYDRO.get(c, 0.0) for c in seq]
    win = [None] * n
    for i in range(n):
        lo = max(0, i - KD_WINDOW // 2)
        hi = min(n, i + KD_WINDOW // 2 + 1)
        win[i] = sum(kd[lo:hi]) / (hi - lo)

    # hydrophobic runs -> candidate TM segments
    runs, cur = [], None
    for i, v in enumerate(win):
        if v > KD_THRESHOLD:
            cur = [i, i] if cur is None else [cur[0], i]
        else:
            if cur is not None:
                runs.append(cur)
                cur = None
    if cur:
        runs.append(cur)
    runs = [r for r in runs if r[1] - r[0] + 1 >= MIN_TM_RUN]
    if not runs:
        return None

    topo = {i: "unknown" for i in range(1, n + 1)}
    tm_marks = set()
    for s, e in runs:
        for i in range(s, e + 1):
            topo[i + 1] = "TM"
            tm_marks.add(i)
    first = min(s for s, _ in runs)
    for i in range(0, first):
        topo[i + 1] = "N-term"
    # Label a C-tail only if enough hydrophobic runs were found to be
    # plausibly a 7TM GPCR; otherwise keep the rest unknown (conservative:
    # prediction then stays silent instead of flooding a fake tail).
    if len(runs) >= 4:
        last = max(e for _, e in runs)
        for i in range(last + 1, n):
            topo[i + 1] = "C-tail"
    return topo


def build_topology_map(entry):
    features = entry.get("features", [])
    seq = entry["sequence"]["value"]
    seq_len = entry["sequence"]["length"]
    topology = {i: "unknown" for i in range(1, seq_len + 1)}

    tms = []      # (start, end)
    topo_domains = []  # (start, end, label)

    for f in features:
        ftype = f.get("type", "")
        loc = f.get("location", {})
        s = loc.get("start", {}).get("value")
        e = loc.get("end", {}).get("value")
        if not s or not e:
            continue
        if ftype.lower() == "transmembrane":
            tms.append((s, e))
        elif ftype.lower() == "topological domain":
            desc = _desc_of(f)
            if any(k in desc for k in ("cytoplasmic", "cytosol", "inside")):
                topo_domains.append((s, e, "ICL"))
            elif any(k in desc for k in ("extracellular", "luminal", "periplasmic", "outside")):
                topo_domains.append((s, e, "ECL"))

    # No annotation at all -> hydrophobicity fallback
    if not tms and not topo_domains:
        fb = _kd_window_topology(seq)
        if fb:
            return fb
        return topology

    # Transmembrane regions first
    for s, e in tms:
        for i in range(s, e + 1):
            topology[i] = "TM"

    # Topological domains: first -> N-term, last -> C-tail, middle verbatim
    if topo_domains:
        topo_domains.sort()
        ordered = []
        for idx, (s, e, lab) in enumerate(topo_domains):
            label = "N-term" if idx == 0 else ("C-tail" if idx == len(topo_domains) - 1 else lab)
            ordered.append((s, e, label))
        for s, e, lab in ordered:
            for i in range(s, e + 1):
                if topology[i] == "unknown":
                    topology[i] = lab

    # Gap inference between known TMs
    tms.sort()
    if tms:
        first_tm, last_tm = tms[0][0], tms[-1][1]
        for i in range(1, first_tm):
            if topology[i] == "unknown":
                topology[i] = "N-term"
        for i in range(last_tm + 1, seq_len + 1):
            if topology[i] == "unknown":
                topology[i] = "C-tail"
        for idx, (s, e) in enumerate(tms):
            if idx + 1 >= len(tms):
                break
            gap_start, gap_end = e + 1, tms[idx + 1][0] - 1
            label = "ICL" if (idx + 1) % 2 == 1 else "ECL"  # after TM1->ICL1, TM2->ECL1
            for i in range(gap_start, gap_end + 1):
                if topology[i] == "unknown":
                    topology[i] = label

    return topology
