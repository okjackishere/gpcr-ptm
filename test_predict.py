"""
Minimal tests for predict.py rule enhancements (pure assert, no pytest).
Run:  python test_predict.py   (prints "all tests passed" on success)

Covers: pXpp motif scan, C-tail distal/proximal annotation, barcode cluster
reason, conservative ubiquitination (PPxY strong / Pro-rich weak / topology
gate), and a regression check that new motifs are registered in MOTIF_BASE
(no KeyError when predict_ptms runs all four PTM types).
"""
from predict import (
    predict_phospho, predict_ubiquitination, predict_ptms,
    MOTIF_BASE,
)
from gpcr_classify import classify_gpcr_group, pxpp_filter, ARRESTIN_CLASS


def _topo(seq_len, tm_end, tail_start):
    """Build {pos: region}: 1..tm_end=TM, tail_start..seq_len=C-tail."""
    topo = {}
    for i in range(1, tm_end + 1):
        topo[i] = "TM"
    for i in range(tail_start, seq_len + 1):
        topo[i] = "C-tail"
    return topo


def test_pxpp_detection_and_proximal_and_cluster():
    # 20 TM (M) + "AAASASSAAA": S at pos 24, 26, 27 -> S-A-S-S pXpp 4-mer
    seq = "M" * 20 + "AAASASSAAA"           # length 30; S@24(idx23),26,27
    topo = _topo(30, 20, 21)
    preds = predict_phospho(seq, topo)
    pxpp = [p for p in preds if "pXpp" in p["motif"]]
    assert pxpp, "pXpp motif should fire on the S-A-S-S cluster"
    pos = {p["position"] for p in pxpp}
    assert {24, 26, 27} <= pos, f"pXpp positions wrong: {pos}"
    # C-tail + tm7=20 -> dist 4/6/7 -> proximal; cluster of 3 S/T -> barcode reason
    for p in pxpp:
        joined = " ".join(p["reasons"])
        assert "近端proximal" in joined, f"missing proximal: {p['reasons']}"
        assert "多磷酸化簇" in joined, f"missing cluster reason: {p['reasons']}"
    print("  ok: pXpp fires at 24/26/27; proximal + cluster reasons present")


def test_pxpp_negative_isolated():
    # isolated single S, no cluster -> no pXpp
    seq = "M" * 20 + "AAASAAAAAAAA"        # single S @ pos 24
    topo = _topo(31, 20, 21)
    preds = predict_phospho(seq, topo)
    assert not any("pXpp" in p["motif"] for p in preds), "isolated S must not hit pXpp"
    print("  ok: isolated S does not hit pXpp")


def test_distal_annotation():
    # TM 1-20, then 20 A, then S-A-S-S at pos 41/43/44 -> dist 21 -> distal
    seq = "M" * 20 + "A" * 20 + "SASS"     # length 44; S@41,43,44
    topo = _topo(44, 20, 21)
    preds = predict_phospho(seq, topo)
    pxpp = [p for p in preds if "pXpp" in p["motif"]]
    assert pxpp and 41 in {p["position"] for p in pxpp}
    joined = " ".join(pxpp[0]["reasons"])
    assert "远端distal" in joined, f"missing distal: {pxpp[0]['reasons']}"
    print("  ok: far C-tail site annotated 远端distal")


def test_ubi_ppxy_strong():
    # C-tail "PPAYK": K near PPxY -> strong Nedd4-docking signal, base 0.50
    seq = "M" * 20 + "PPAYKAAA"            # K@25 (idx24); PPxY @ 21-24
    topo = _topo(28, 20, 21)
    preds = predict_ubiquitination(seq, topo)
    assert len(preds) == 1 and preds[0]["position"] == 25, preds
    p = preds[0]
    assert p["base"] == 0.50, p
    assert "PPxY" in p["motif"] and "PPxY" in " ".join(p["reasons"]), p
    print("  ok: K near PPxY flagged base=0.50 (Nedd4 docking)")


def test_ubi_prorich_weak():
    # "PPPPK": Pro-rich but no PPxY (no Y) -> weak signal, base 0.40
    seq = "M" * 20 + "PPPPKAAA"            # K@25; 4 Pro, no Tyr -> no PPxY
    topo = _topo(28, 20, 21)
    preds = predict_ubiquitination(seq, topo)
    assert len(preds) == 1 and preds[0]["position"] == 25, preds
    p = preds[0]
    assert p["base"] == 0.40, p
    assert "Pro-rich" in p["motif"], p
    print("  ok: K in Pro-rich (no PPxY) flagged base=0.40")


def test_ubi_negative_no_signal():
    # "AAAAK": no Pro, no PPxY -> not flagged
    seq = "M" * 20 + "AAAAKAAA"
    topo = _topo(28, 20, 21)
    preds = predict_ubiquitination(seq, topo)
    assert preds == [], f"should be empty: {preds}"
    print("  ok: K without Pro/PPxY not flagged")


def test_ubi_topology_gate_extracellular():
    # K near PPxY but in an extracellular region -> not flagged
    seq = "M" * 20 + "PPAYKAAA"
    topo = {i: "TM" for i in range(1, 21)}
    for i in range(21, 29):
        topo[i] = "ECL"                    # extracellular -> out of scope
    preds = predict_ubiquitination(seq, topo)
    assert preds == [], f"extracellular K must not be flagged: {preds}"
    print("  ok: extracellular K not flagged (topology gate)")


def test_regression_no_keyerror_all_types():
    # new motifs (pXpp) must be registered in MOTIF_BASE, else KeyError here
    assert "pXpp" in MOTIF_BASE, "pXpp missing from MOTIF_BASE"
    seq = "M" * 20 + "AAASASSAAAPPAYKPPP"  # mixed phospho + ubi signals
    topo = _topo(37, 20, 21)
    for ptm in ("phosphorylation", "glycosylation", "palmitoylation", "ubiquitination"):
        out = predict_ptms(ptm, seq, topo, {})   # must not raise
        assert isinstance(out, list)
    print("  ok: predict_ptms runs all 4 types without KeyError")


# ---- pXpp 判组过滤测试 (Isaikina 2023 Table S2 + ICL3 长度回退) ----

def _topo_with_icl3(seq_len, tm_end, icl3_start, icl3_end, tail_start):
    """建拓扑: 1..tm_end=TM, icl3 区间=ICL, tail_start..end=C-tail。"""
    topo = {}
    for i in range(1, tm_end + 1):
        topo[i] = "TM"
    for i in range(icl3_start, icl3_end + 1):
        topo[i] = "ICL"
    for i in range(tail_start, seq_len + 1):
        topo[i] = "C-tail"
    # 其余填 unknown
    for i in range(1, seq_len + 1):
        topo.setdefault(i, "unknown")
    return topo


def test_table_s2_mapping():
    """Table S2 映射表命中: ADRB2=A, CCR5=B, CXCR3=magenta。"""
    assert ARRESTIN_CLASS["P07550"] == "A", "ADRB2 应为 class A"
    assert ARRESTIN_CLASS["P51681"] == "B", "CCR5 应为 class B"
    assert ARRESTIN_CLASS["P49682"] == "magenta", "CXCR3 应为 magenta"
    print("  ok: Table S2 映射正确 (ADRB2=A, CCR5=B, CXCR3=magenta)")


def test_classify_table_s2_hit():
    """Table S2 命中时用一手实验数据, confidence=table_s2。"""
    g = classify_gpcr_group({"primaryAccession": "P07550"}, {})
    assert g["group"] == "A" and g["confidence"] == "table_s2", g
    g = classify_gpcr_group({"primaryAccession": "P51681"}, {})
    assert g["group"] == "B" and g["confidence"] == "table_s2", g
    print("  ok: Table S2 命中返回正确 class + confidence")


def test_classify_icl3_length_fallback():
    """未在 Table S2 时用 ICL3 长度回退: 极短→B, 极长→A, 中间→None。"""
    # ICL3=3aa (极短 → B 倾向)
    topo_short = _topo_with_icl3(40, 20, 25, 27, 30)
    g = classify_gpcr_group({"primaryAccession": "XXXXXX"}, topo_short)
    assert g["group"] == "B" and g["confidence"] == "icl3_length", g
    # ICL3=120aa (极长 → A 倾向)
    topo_long = _topo_with_icl3(200, 20, 25, 144, 150)
    g = classify_gpcr_group({"primaryAccession": "XXXXXX"}, topo_long)
    assert g["group"] == "A" and g["confidence"] == "icl3_length", g
    # ICL3=30aa (中间 → None)
    topo_mid = _topo_with_icl3(60, 20, 25, 54, 55)
    g = classify_gpcr_group({"primaryAccession": "XXXXXX"}, topo_mid)
    assert g["group"] is None and g["confidence"] == "unknown", g
    print("  ok: ICL3 长度回退正确 (<5→B, >100→A, 中间→None)")


def test_pxpp_filter_class_a_blocks_ctail():
    """Class A 组: C-tail 的 pXpp 应被过滤, ICL3 的保留。"""
    # seq: TM(1-20) + ICL3(21-32, 含SASS@25-28) + C-tail(33+, 含TSTT@37-40)
    seq = "M" * 20 + "AAAASASSAAAA" + "AAAATSTTAAAA"
    topo = _topo_with_icl3(len(seq), 20, 21, 32, 33)
    group = {"group": "A", "confidence": "table_s2", "caveat": "test"}
    preds = predict_phospho(seq, topo, group=group)
    pxpp = [p for p in preds if "pXpp" in p["motif"]]
    ctail_pxpp = [p for p in pxpp if p["region"] == "C-tail"]
    icl_pxpp = [p for p in pxpp if p["region"] == "ICL"]
    assert len(ctail_pxpp) == 0, f"class A C-tail pXpp 应被过滤: {ctail_pxpp}"
    assert len(icl_pxpp) > 0, f"class A ICL3 pXpp 应保留: {pxpp}"
    print(f"  ok: class A 过滤 C-tail pXpp ({len(ctail_pxpp)} C-tail, {len(icl_pxpp)} ICL3)")


def test_pxpp_filter_class_b_blocks_icl():
    """Class B 组: ICL3 的 pXpp 应被过滤, C-tail 的保留。"""
    seq = "M" * 20 + "AAAASASSAAAA" + "AAAATSTTAAAA"
    topo = _topo_with_icl3(len(seq), 20, 21, 32, 33)
    group = {"group": "B", "confidence": "table_s2", "caveat": "test"}
    preds = predict_phospho(seq, topo, group=group)
    pxpp = [p for p in preds if "pXpp" in p["motif"]]
    ctail_pxpp = [p for p in pxpp if p["region"] == "C-tail"]
    icl_pxpp = [p for p in pxpp if p["region"] == "ICL"]
    assert len(icl_pxpp) == 0, f"class B ICL3 pXpp 应被过滤: {icl_pxpp}"
    assert len(ctail_pxpp) > 0, f"class B C-tail pXpp 应保留: {pxpp}"
    print(f"  ok: class B 过滤 ICL3 pXpp ({len(ctail_pxpp)} C-tail, {len(icl_pxpp)} ICL3)")


def test_pxpp_filter_none_keeps_both():
    """无判组(None/magenta): 两区 pXpp 都保留 (现状不变)。"""
    seq = "M" * 20 + "AAAASASSAAAA" + "AAAATSTTAAAA"
    topo = _topo_with_icl3(len(seq), 20, 21, 32, 33)
    # group=None
    preds = predict_phospho(seq, topo, group={"group": None, "confidence": "unknown", "caveat": ""})
    pxpp = [p for p in preds if "pXpp" in p["motif"]]
    ctail_pxpp = [p for p in pxpp if p["region"] == "C-tail"]
    icl_pxpp = [p for p in pxpp if p["region"] == "ICL"]
    assert len(ctail_pxpp) > 0 and len(icl_pxpp) > 0, f"None 组两区都应保留: {pxpp}"
    # magenta
    preds2 = predict_phospho(seq, topo, group={"group": "magenta", "confidence": "table_s2", "caveat": ""})
    pxpp2 = [p for p in preds2 if "pXpp" in p["motif"]]
    ctail2 = [p for p in pxpp2 if p["region"] == "C-tail"]
    icl2 = [p for p in pxpp2 if p["region"] == "ICL"]
    assert len(ctail2) > 0 and len(icl2) > 0, f"magenta 两区都应保留: {pxpp2}"
    print("  ok: None/magenta 两区 pXpp 都保留")


def test_ck1_priming_direction():
    """CK1 共识 pS/pT-x-x-S/T: 潜在引物须在候选位点上游 i-3 (Venerando 2014)。
    修正前代码检查 i+3, 方向相反 (把引物位点错当靶位点)。"""
    # S@21 为潜在引物, S@24 在其 +3 下游 -> CK1 应命中 24 (而非 21)
    seq = "M" * 20 + "SAASTAAAA"
    preds = predict_phospho(seq, _topo(30, 20, 21))
    ck1 = {p["position"] for p in preds if "CK1" in p["motif"]}
    assert 24 in ck1, f"CK1 应命中引物下游的 24: {ck1}"
    assert 21 not in ck1, f"引物位点自身(21)不应因下游 S/T 命中 CK1: {ck1}"
    # 反向验证: T@24 上游(21)无 S/T -> 24 不命中; S@27 上游 24 有 T -> 27 命中
    seq2 = "M" * 20 + "AAATAASAAA"
    preds2 = predict_phospho(seq2, _topo(30, 20, 21))
    ck1_2 = {p["position"] for p in preds2 if "CK1" in p["motif"]}
    assert 24 not in ck1_2, f"上游无引物的 24 不应命中 CK1 (旧 i+3 检查的误报): {ck1_2}"
    assert 27 in ck1_2, f"CK1 应命中有上游引物(T@24)的 27: {ck1_2}"
    print("  ok: CK1 引物方向为 i-3 (24 引物→27 命中; 上游无引物不误报)")


def main():
    tests = [
        test_pxpp_detection_and_proximal_and_cluster,
        test_pxpp_negative_isolated,
        test_distal_annotation,
        test_ubi_ppxy_strong,
        test_ubi_prorich_weak,
        test_ubi_negative_no_signal,
        test_ubi_topology_gate_extracellular,
        test_regression_no_keyerror_all_types,
        test_ck1_priming_direction,
        test_table_s2_mapping,
        test_classify_table_s2_hit,
        test_classify_icl3_length_fallback,
        test_pxpp_filter_class_a_blocks_ctail,
        test_pxpp_filter_class_b_blocks_icl,
        test_pxpp_filter_none_keeps_both,
    ]
    for t in tests:
        print(f"[{t.__name__}]")
        t()
    print("\nall tests passed")


if __name__ == "__main__":
    main()
