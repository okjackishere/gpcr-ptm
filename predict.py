"""
Rule-based PTM prediction, restricted to sites with a real consensus motif AND
the correct membrane topology. No bare-residue ("all S/T in the tail") outputs.

Caveats (scientific):
  - A consensus motif is necessary but NOT sufficient: e.g. only ~30-40% of
    N-glycosylation sequons are actually occupied; phospho motifs are weak.
    Scores are for ranking only, never probabilities.
  - GRK has no reliable linear motif; the "GRK-candidate" rule below is only a
    weak positional hint (juxtamembrane intracellular Ser/Thr with nearby
    acidic residues), flagged low confidence.
  - Ubiquitination has no simple linear consensus -> intentionally NOT predicted.
"""

ALLOWED_PHOSPHO = ("ICL", "C-tail")
ALLOWED_GLYCO = ("N-term", "ECL")
ALLOWED_PALM = ("ICL", "C-tail")

# position i is a 1-based S/T. Indices below are 0-based into `sequence`:
#   residue at 1-based p  ==  sequence[p-1]
MOTIF_BASE = {
    "PKA": 0.75,
    "CDK/MAPK": 0.65,
    "proline-directed": 0.60,
    "PKC": 0.55,
    "CK2": 0.50,
    "GRK-candidate": 0.50,
    "CK1": 0.40,
}

MOTIF_DESC = {
    "PKA": "PKA: [RK][RK]x[S/T]",
    "CDK/MAPK": "CDK/MAPK: [S/T]Px[RK]",
    "proline-directed": "proline-directed: [S/T]P",
    "PKC": "PKC: [S/T]x[R/K]",
    "CK2": "CK2: [S/T]xx[D/E]",
    "CK1": "CK1: [S/T]xx[pS/pT] (需引物磷酸化)",
    "GRK-candidate": "GRK无可靠线性motif, 仅位置提示(低置信)",
}


def _count(seq, chars):
    return sum(1 for c in seq if c in chars)


def predict_phospho(sequence, topology):
    preds = []
    seq_len = len(sequence)

    for i in range(1, seq_len + 1):  # i = 1-based S/T position
        if sequence[i - 1] not in "ST":
            continue
        region = topology.get(i, "unknown")
        if region not in ALLOWED_PHOSPHO:
            continue

        motifs, reasons = [], []
        pkc_hydro = False
        idx = i - 1  # 0-based index of the S/T

        # PKA: [RK][RK]x[S/T]  -> basic at 1-based i-3,i-2
        if i >= 4 and sequence[idx - 3] in "RK" and sequence[idx - 2] in "RK":
            motifs.append("PKA")
        # PKC: [S/T]x[R/K]  -> basic at i+1; hydrophobic at i-1 boosts confidence
        if i < seq_len and sequence[idx + 1] in "RK":
            motifs.append("PKC")
            if i >= 2 and sequence[idx - 1] in "IVLMFWY":
                pkc_hydro = True
                reasons.append("i-1为疏水残基(PKC偏好)")
        # CK2: [S/T]xx[D/E] -> acidic at i+2 or i+3
        if any(c in "DE" for c in sequence[idx + 2:idx + 4]):
            motifs.append("CK2")
        # CK1 priming: [S/T]xx[pS/pT] -> S/T at i+3
        if idx + 3 < seq_len and sequence[idx + 3] in "ST":
            motifs.append("CK1")
        # proline-directed: [S/T]P;  [S/T]Px[RK] is a better CDK/MAPK hit
        if i < seq_len and sequence[idx + 1] == "P":
            motifs.append("CDK/MAPK" if (i + 1 < seq_len and sequence[idx + 2] in "RK")
                          else "proline-directed")

        # GRK: no linear motif; weak positional hint = intracellular S/T with
        # >=2 acidic residues within +/-7 and >=2 S/T within +/-10
        up = sequence[max(0, idx - 7):idx]
        down = sequence[idx + 1:min(seq_len, idx + 8)]
        acidic = _count(up + down, "DE")
        st_cluster = _count(sequence[max(0, idx - 10):min(seq_len, idx + 11)], "ST")
        grk_hint = acidic >= 2 and st_cluster >= 2

        for m in motifs:
            base = MOTIF_BASE[m]
            if m == "PKC" and pkc_hydro:
                base += 0.05
            r = [f"{MOTIF_DESC[m]} (区域:{region})"] + list(reasons)
            if grk_hint:
                base = max(base, 0.5)
                r.append("邻近酸性残基>=2 且 S/T簇集 (GRK候选提示)")
                m = f"GRK+{m}" if m != "GRK-candidate" else m
            preds.append(_mk(sequence, topology, "phosphorylation", i, m, base, r))

        if grk_hint and not motifs:
            preds.append(_mk(
                sequence, topology, "phosphorylation", i, "GRK-candidate", 0.5,
                [f"{MOTIF_DESC['GRK-candidate']} (区域:{region})",
                 "邻近酸性残基>=2 且 S/T簇集"]))

    return preds


def predict_glyco(sequence, topology):
    preds = []
    seq_len = len(sequence)
    for i in range(seq_len - 2):
        if sequence[i] != "N":
            continue
        if sequence[i + 1] == "P":
            continue
        if sequence[i + 2] not in "ST":
            continue
        pos = i + 1
        region = topology.get(pos, "unknown")
        if region not in ALLOWED_GLYCO:
            continue
        preds.append(_mk(sequence, topology, "glycosylation", pos,
                         f"N-{sequence[i+1]}-{sequence[i+2]}", 0.6,
                         ["N-X-S/T sequon (X!=P)", f"区域:{region}",
                          "注意: 约30-40% sequon实际不占用, 需实验确认"]))
    return preds


def predict_palm(sequence, topology, tm7_end=None):
    preds = []
    seq_len = len(sequence)
    if tm7_end is None:
        for pos in range(len(topology), 0, -1):
            if topology.get(pos) == "TM":
                tm7_end = pos
                break
    for i in range(1, seq_len + 1):
        if sequence[i - 1] != "C":
            continue
        region = topology.get(i, "unknown")
        if region not in ALLOWED_PALM:
            continue
        win = sequence[max(0, i - 11):min(seq_len, i + 10)]
        kr = _count(win, "KR")
        juxta = tm7_end is not None and 0 < (i - tm7_end) <= 25
        if kr < 2 and not juxta:
            continue
        base, reasons = 0.5, [f"区域:{region}"]
        if juxta:
            base += 0.05
            reasons.append(f"TM7后{(i-tm7_end)}aa的膜旁区(经典棕榈酰化位置)")
        if kr >= 2:
            reasons.append(f"+-10aa内{kr}个碱性K/R残基(利于膜表面酰化)")
        preds.append(_mk(sequence, topology, "palmitoylation", i,
                         "juxtamembrane+Cys" if juxta else "Cys+K/R-rich",
                         base, reasons))
    return preds


def _ctx(seq, pos, w=5):
    start = max(0, pos - 1 - w)
    end = min(len(seq), pos + w)
    return seq[start:end]


def _mk(sequence, topology, ptm_type, pos, motif, base, reasons):
    return {
        "ptm_type": ptm_type,
        "position": pos,
        "residue": sequence[pos - 1],
        "region": topology.get(pos, "unknown"),
        "layer": "Predicted",
        "motif": motif,
        "base": round(base, 3),
        "reasons": reasons,
        "context": _ctx(sequence, pos),
        "conservation": None,
        "score": 0.0,
    }


def predict_ptms(ptm_type, sequence, topology, entry):
    if ptm_type == "phosphorylation":
        return predict_phospho(sequence, topology)
    if ptm_type == "glycosylation":
        return predict_glyco(sequence, topology)
    if ptm_type == "palmitoylation":
        return predict_palm(sequence, topology)
    return []
