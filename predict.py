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
  - Ubiquitination: no simple linear consensus, but a conservative predictor
    flags intracellular Lys near a PPxY motif (Nedd4 WW-domain docking) or a
    Pro-rich stretch; always Low confidence, a weak hint only. PPxY is often
    on adaptors (e.g. ARRDC3) rather than the GPCR itself, so adaptor-mediated
    ubiquitination is missed -> flagged "needs experimental confirmation".
"""

ALLOWED_PHOSPHO = ("ICL", "C-tail")
ALLOWED_GLYCO = ("N-term", "ECL")
ALLOWED_PALM = ("ICL", "C-tail")
ALLOWED_UBIQ = ("ICL", "C-tail")

# position i is a 1-based S/T. Indices below are 0-based into `sequence`:
#   residue at 1-based p  ==  sequence[p-1]
MOTIF_BASE = {
    "PKA": 0.75,
    "CDK/MAPK": 0.65,
    "proline-directed": 0.60,
    "pXpp": 0.60,
    "PKC": 0.55,
    "CK2": 0.50,
    "GRK-candidate": 0.50,
    "CK1": 0.40,
}

MOTIF_DESC = {
    "PKA": "PKA: [RK][RK]x[S/T]",
    "CDK/MAPK": "CDK/MAPK: [S/T]Px[RK]",
    "proline-directed": "proline-directed: [S/T]P",
    "pXpp": "pXpp: [S/T]x[S/T][S/T] (Isaikina 2023, arrestin2募集保守motif)",
    "PKC": "PKC: [S/T]x[R/K]",
    "CK2": "CK2: [S/T]xx[D/E]",
    "CK1": "CK1: [S/T]xx[pS/pT] (需引物磷酸化)",
    "GRK-candidate": "GRK无可靠线性motif, 仅位置提示(低置信)",
}

# 每条规则的支撑文献(可点击核验)。来源: AMiner 一手文献检索 + searchPro 交叉验证。
# 每条 = {short, title, url}；url 优先 DOI，无 DOI 用 AMiner 链接。
REFS = {
    "pXpp": [{"short": "Isaikina 2023·Mol Cell",
              "title": "A Key GPCR Phosphorylation Motif Discovered in Arrestin2·CCR5 Phosphopeptide Complexes",
              "url": "https://doi.org/10.1016/j.molcel.2023.05.002"}],
    "barcode": [{"short": "Latorraca 2020·Cell",
                 "title": "How GPCR Phosphorylation Patterns Orchestrate Arrestin-Mediated Signaling",
                 "url": "https://doi.org/10.1016/j.cell.2020.11.014"}],
    "distal_proximal": [{"short": "Sente 2018·Nat Struct Mol Biol",
                         "title": "Molecular Mechanism of Modulating Arrestin Conformation by GPCR Phosphorylation",
                         "url": "https://doi.org/10.1038/s41594-018-0071-3"}],
    "GRK": [{"short": "Ribas 2007·BBA",
             "title": "The G Protein-Coupled Receptor Kinase (GRK) Interactome",
             "url": "https://doi.org/10.1016/j.bbamem.2006.09.019"}],
    "kinase": [{"short": "Xue 2005·NAR (GPS)",
                "title": "GPS: a Comprehensive WWW Server for Phosphorylation Sites Prediction",
                "url": "https://doi.org/10.1093/nar/gki393"}],
    "glycosylation": [{"short": "Rodriguez 1995·JBC",
                       "title": "Role of N-Glycosylation for Functional Expression of the Human PAF Receptor",
                       "url": "https://doi.org/10.1074/jbc.270.42.25178"}],
    "palmitoylation": [{"short": "Hussain 2018·Anal Biochem",
                        "title": "SPalmitoylC-PseAAC: identifying S-palmitoylation sites in proteins",
                        "url": "https://doi.org/10.1016/j.ab.2018.12.019"}],
    "ubiquitination": [
        {"short": "Kennedy & Marchese 2015",
         "title": "Regulation of GPCR Trafficking by Ubiquitin",
         "url": "https://doi.org/10.1016/bs.pmbts.2015.02.005"},
        {"short": "Min 2012",
         "title": "In silico identification and in vitro validation of Nedd4-mediated GPCR ubiquitination",
         "url": "https://www.aminer.cn/pub/56d87a47dabfae2eee362f64"}],
}


def _count(seq, chars):
    return sum(1 for c in seq if c in chars)


def _tm7_end(topology):
    """最后一个 TM 残基的 1-based 位置；无 TM 注释返回 None。"""
    for pos in range(len(topology), 0, -1):
        if topology.get(pos) == "TM":
            return pos
    return None


def _is_pxpp(seq, idx, n):
    """候选位点(idx, 0-based S/T) 是否属于 [S/T]-X-[S/T]-[S/T] 4-mer 中的某个磷酸位。
    pXpp motif: Isaikina 2023 (Mol Cell), GPCR arrestin2 募集保守 motif;
    p = 磷酸化 S/T (非脯氨酸), X = 任意残基。"""
    s = "ST"
    # 候选=pos0: [ST]-X-[ST]-[ST]
    if idx + 3 < n and seq[idx + 2] in s and seq[idx + 3] in s:
        return True
    # 候选=pos2: [ST]-X-[ST]-[ST]
    if 2 <= idx < n - 1 and seq[idx - 2] in s and seq[idx + 1] in s:
        return True
    # 候选=pos3: [ST]-X-[ST]-[ST]
    if idx >= 3 and seq[idx - 3] in s and seq[idx - 1] in s:
        return True
    return False


def _has_ppxy(win):
    """窗口内是否含 PPxY (P-P-x-Y) 4-mer，Nedd4 WW 结构域识别 motif。"""
    for j in range(len(win) - 3):
        if win[j] == "P" and win[j + 1] == "P" and win[j + 3] == "Y":
            return True
    return False


def _dedup_refs(refs):
    """按 url 去重文献列表，保持顺序。"""
    seen, out = set(), []
    for x in refs:
        if x["url"] not in seen:
            seen.add(x["url"])
            out.append(x)
    return out


def _motif_refs(motif_name):
    """激酶 motif 的支撑文献: pXpp 用 Isaikina 2023, 其余通用激酶共识用 GPS 2005。"""
    if "pXpp" in motif_name:
        return REFS["pXpp"]
    return REFS["kinase"]


def predict_phospho(sequence, topology, group=None):
    preds = []
    seq_len = len(sequence)
    tm7 = _tm7_end(topology)
    from gpcr_classify import pxpp_filter

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
        # pXpp: [S/T]-X-[S/T]-[S/T], Isaikina 2023 (Mol Cell) arrestin2 募集保守 motif
        # 按判组过滤(而非调权): 原文表明某些组的特定区域 pXpp 缺失, 故直接不报——
        # 调权只改排序不改命中集合, 过滤才能实质移除假阳性
        if _is_pxpp(sequence, idx, seq_len):
            grp = group.get("group") if group else None
            if pxpp_filter(grp, region):
                motifs.append("pXpp")
                # 在 reasons 里标注判组依据（通俗易懂）
                grp_caveat = group.get("caveat", "") if group else ""
                if grp == "B":
                    reasons.append(f"该pXpp命中保留：受体在Table S2中被归类为class B（C-tail富集组），C-tail正是该组pXpp的典型分布区（{grp_caveat}）")
                elif grp == "A":
                    reasons.append(f"该pXpp命中保留：受体在Table S2中被归类为class A（ICL3富集组），ICL3正是该组pXpp的典型分布区（{grp_caveat}）")
                elif grp == "magenta":
                    reasons.append(f"该受体pXpp分布因亚型而异（magenta类），两区pXpp均暂保留，需结合具体亚型判断（{grp_caveat}）")

        # GRK: no linear motif; weak positional hint = intracellular S/T with
        # >=2 acidic residues within +/-7 and >=2 S/T within +/-10
        up = sequence[max(0, idx - 7):idx]
        down = sequence[idx + 1:min(seq_len, idx + 8)]
        acidic = _count(up + down, "DE")
        st_cluster = _count(sequence[max(0, idx - 10):min(seq_len, idx + 11)], "ST")
        grk_hint = acidic >= 2 and st_cluster >= 2

        # C-tail / ICL positional context: distal vs proximal (Sente 2018,
        # Nat Struct Mol Biol) and multi-S/T barcode cluster (Latorraca 2020, Cell)
        ctx_reasons, ctx_refs = [], []
        if region == "C-tail" and tm7 is not None:
            d = i - tm7
            ctx_reasons.append(
                f"距TM7 {d}aa（{'近端proximal' if d <= 15 else '远端distal'}区）")
            ctx_refs += REFS["distal_proximal"]
        if st_cluster >= 3:
            ctx_reasons.append(f"多磷酸化簇(±10内{st_cluster}个S/T)")
            ctx_refs += REFS["barcode"]

        for m in motifs:
            base = MOTIF_BASE[m]
            if m == "PKC" and pkc_hydro:
                base += 0.05
            r = [f"{MOTIF_DESC[m]} (区域:{region})"] + list(reasons) + ctx_reasons
            refs = _motif_refs(m) + ctx_refs
            if grk_hint:
                base = max(base, 0.5)
                r.append("邻近酸性残基>=2 且 S/T簇集 (GRK候选提示)")
                refs = refs + REFS["GRK"]
                m = f"GRK+{m}" if m != "GRK-candidate" else m
            preds.append(_mk(sequence, topology, "phosphorylation", i, m, base, r,
                              _dedup_refs(refs)))

        if grk_hint and not motifs:
            preds.append(_mk(
                sequence, topology, "phosphorylation", i, "GRK-candidate", 0.5,
                [f"{MOTIF_DESC['GRK-candidate']} (区域:{region})",
                 "邻近酸性残基>=2 且 S/T簇集"] + ctx_reasons,
                _dedup_refs(REFS["GRK"] + ctx_refs)))

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
                          "注意: 约30-40% sequon实际不占用, 需实验确认"],
                         REFS["glycosylation"]))
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
                         base, reasons, REFS["palmitoylation"]))
    return preds


def predict_ubiquitination(sequence, topology):
    """保守版泛素化预测: 胞内 Lys + 邻近 PPxY(Nedd4 WW域直接docking强信号) 或 Pro-rich(弱信号)。
    无线性 consensus, 恒 Low 置信, 仅弱提示。PPxY 常在 adaptor(ARRDC3)上而非 GPCR 自身,
    adaptor 介导的泛素化会漏报, 需实验确认。"""
    preds = []
    seq_len = len(sequence)
    for i in range(1, seq_len + 1):
        if sequence[i - 1] != "K":
            continue
        region = topology.get(i, "unknown")
        if region not in ALLOWED_UBIQ:
            continue
        idx = i - 1
        win15 = sequence[max(0, idx - 15):min(seq_len, idx + 16)]
        win10 = sequence[max(0, idx - 10):min(seq_len, idx + 11)]
        if _has_ppxy(win15):
            preds.append(_mk(sequence, topology, "ubiquitination", i,
                             "K+PPxY(Nedd4候选)", 0.50,
                             ["K + 胞内区",
                              "邻近PPxY motif(±15内, Nedd4 WW域直接docking信号)",
                              "无可靠线性motif, 弱位置提示, 需实验确认"],
                             REFS["ubiquitination"]))
        else:
            pro = _count(win10, "P")
            if pro < 2:
                continue
            preds.append(_mk(sequence, topology, "ubiquitination", i,
                             "K+Pro-rich(弱)", 0.40,
                             ["K + 胞内区",
                              f"邻近Pro-rich(±10内{pro}个P, Nedd4 WW结构域识别提示)",
                              "无可靠线性motif, 弱位置提示, 需实验确认"],
                             REFS["ubiquitination"]))
    return preds


def _ctx(seq, pos, w=5):
    start = max(0, pos - 1 - w)
    end = min(len(seq), pos + w)
    return seq[start:end]


def _mk(sequence, topology, ptm_type, pos, motif, base, reasons, refs=None):
    return {
        "ptm_type": ptm_type,
        "position": pos,
        "residue": sequence[pos - 1],
        "region": topology.get(pos, "unknown"),
        "layer": "Predicted",
        "motif": motif,
        "base": round(base, 3),
        "reasons": reasons,
        "refs": refs or [],
        "context": _ctx(sequence, pos),
        "conservation": None,
        "score": 0.0,
    }


def predict_ptms(ptm_type, sequence, topology, entry, group=None):
    if ptm_type == "phosphorylation":
        return predict_phospho(sequence, topology, group=group)
    if ptm_type == "glycosylation":
        return predict_glyco(sequence, topology)
    if ptm_type == "palmitoylation":
        return predict_palm(sequence, topology)
    if ptm_type == "ubiquitination":
        return predict_ubiquitination(sequence, topology)
    return []
