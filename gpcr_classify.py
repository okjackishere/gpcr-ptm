"""
GPCR arrestin-interaction class classification (Class A / B / magenta).

依据: Isaikina et al. 2023, Molecular Cell 的 Table S2 + Figure 7B。
这里的 class A/B 指 arrestin 相互作用行为(A=transient/偏好arrestin3,
B=stable/偏好arrestin2), 不是 GPCRdb 的 GRAFS 分子分类。

判组两层次:
  B(首选): Table S2 内置映射表 (accession → class), 一手实验数据
  A(回退): 最长 ICL3 残基数 (原文阈值 <~5aa 或 >100aa)

用于 predict_phospho 的 pXpp 区域过滤。
"""

# Table S2: accession → class (Isaikina 2023 Mol Cell)
# class A = red (ICL3富集, C-tail pXpp 几乎缺失)
# class B = blue (C-tail富集, ICL3 短)
# CXCR3 = magenta (CXCR3A=B, CXCR3B=A, 亚型依赖)
ARRESTIN_CLASS = {
    # --- Class A (ICL3 富集组) ---
    "P35368": "A",  # α1B-ADR
    "P07550": "A",  # β2-ADR
    "P35372": "A",  # μOR
    "P25101": "A",  # ETAR
    "P21728": "A",  # D1R
    "P14416": "A",  # D2R
    "P28223": "A",  # 5HT2aR
    "P41595": "A",  # 5HT2bR
    "P11229": "A",  # M1R
    "P08172": "A",  # M2R
    "P32745": "A",  # SST3R
    "P35346": "A",  # SST5R
    # CXCR3B 也是 A, 但 CXCR3 整体归 magenta (见下)
    # --- Class B (C-tail 富集组) ---
    "P49682": "magenta",  # CXCR3 (CXCR3A=B, CXCR3B=A, 亚型依赖)
    "P61073": "B",  # CXCR4
    "P51681": "B",  # CCR5
    "P32248": "B",  # CCR7
    "O00590": "B",  # ACKR2
    "P25106": "B",  # ACKR3
    "P30518": "B",  # V2R
    "P30989": "B",  # NTR1
    "Q9P296": "B",  # C5aR2
    "P08100": "B",  # Rho
    "P30556": "B",  # AT1R
    "P30559": "B",  # OTR
    "P30874": "B",  # SST2R
    "P34981": "B",  # TRH1R
}

# 原文长度阈值 (逐字引自 Isaikina 2023 正文)
ICL3_SHORT_MAX = 5    # "< ~5aa" → B 倾向 (ICL3 极短)
ICL3_LONG_MIN = 100   # ">100 aa" → A 倾向 (ICL3 极长)


def _max_icl3_length(topology):
    """返回 topology 中最长 ICL3 段的残基数; 无 ICL3 返回 0。"""
    # topology: {pos: region_label}; 找连续的 "ICL" 段最长者
    max_len = cur = 0
    for pos in range(1, len(topology) + 2):
        if topology.get(pos) == "ICL":
            cur += 1
            max_len = max(max_len, cur)
        else:
            cur = 0
    return max_len


def classify_gpcr_group(entry, topology):
    """判组: 返回 dict {group, confidence, caveat}。

    group: "A" / "B" / "magenta" / None(不定)
    confidence: "table_s2" / "icl3_length" / "unknown"
    caveat: 人类可读的判组依据说明
    """
    accession = entry.get("primaryAccession") or entry.get("accession") or ""

    # 层次 B: Table S2 映射表 (首选, 一手实验数据)
    if accession in ARRESTIN_CLASS:
        cls = ARRESTIN_CLASS[accession]
        labels = {"A": "class A (ICL3富集, arrestin3偏好/瞬时)",
                  "B": "class B (C-tail富集, arrestin2偏好/稳定)",
                  "magenta": "亚型依赖 (CXCR3A=B / CXCR3B=A)"}
        return {"group": cls, "confidence": "table_s2",
                "caveat": f"Table S2 实验归类: {labels[cls]}"}

    # 层次 A: ICL3 长度回退 (仅当未在 Table S2)
    icl3_len = _max_icl3_length(topology)
    if icl3_len > 0 and icl3_len <= ICL3_SHORT_MAX:
        return {"group": "B", "confidence": "icl3_length",
                "caveat": f"ICL3极短({icl3_len}aa), 按原文第一组倾向(B组)"}
    if icl3_len >= ICL3_LONG_MIN:
        return {"group": "A", "confidence": "icl3_length",
                "caveat": f"ICL3极长({icl3_len}aa), 按原文第二组倾向(A组)"}

    # 不定: 原文未给 5-100aa 区间的长度判据
    return {"group": None, "confidence": "unknown",
            "caveat": f"未在Table S2, ICL3={icl3_len}aa(中间地带), 无法判组"}


def pxpp_filter(group, region):
    """根据判组结果决定某区域的 pXpp 是否应保留。

    返回 True=保留(正常报), False=过滤(不报)。
    原则: 原文说"该区域 pXpp 缺失"的, 直接过滤而非调权。
    """
    if group == "B":
        # C-tail 富集组: ICL3 的 pXpp 非典型 → 过滤
        return region == "C-tail"
    if group == "A":
        # ICL3 富集组: C-tail 的 pXpp 几乎缺失 → 过滤
        return region == "ICL"
    # magenta / None: 两区都保留
    return True
