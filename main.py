#!/usr/bin/env python3
"""
GPCR-PTM: Predict and verify PTM sites on GPCRs.
Usage: python main.py ADRB2
       python main.py P07550
"""
import argparse
import sys

from resolve import resolve_input
from fetch import fetch_uniprot, fetch_iptmnet
from topology import build_topology_map
from verified import extract_verified
from predict import predict_ptms
from conservation import compute_conservation
from scoring import score_predictions
from output import print_report, save_json, write_verification_md, render_html

MODELS = {
    "phosphorylation": "Phosphorylation",
    "glycosylation": "N-Glycosylation",
    "palmitoylation": "Palmitoylation",
    # 泛素化无可靠线性consensus, 不做规则预测; 仅由数据库给出 L1/L2
}


class AnalysisError(Exception):
    """分析过程的可预期失败(如无法解析输入、UniProt获取失败)。"""


def analyze(gpcr_input, refresh=False, progress=None, record=None):
    """核心分析流水线。返回 (all_results, record, entry, steps)。

    progress: 可选回调 progress(消息字符串), 用于网页端展示进度。
    record:   可选, 若已解析过输入(如网页端做了缓存检查)则直接传入, 避免二次解析。
    """
    steps = []

    def emit(msg):
        steps.append(msg)
        print(f"[*] {msg}")
        if progress:
            progress(msg)

    # Offline data: default uses existing files; only --refresh triggers download
    from download_data import report_status
    if refresh:
        from download_data import download
        emit("刷新离线数据库文件(约120MB, 可能需要数分钟)...")
        try:
            download()
        except Exception as e:
            raise AnalysisError(f"离线数据下载失败: {e}")
    else:
        report_status()

    emit(f"正在解析输入: {gpcr_input}")
    if record is None:
        record = resolve_input(gpcr_input)
    if not record:
        raise AnalysisError(f"未找到 '{gpcr_input}'，请检查基因名 / UniProt ID / 蛋白名称")

    acc = record["accession"]
    name = record["name"]
    gene = record["gene"]
    emit(f"目标蛋白: {gene} ({name}) | {acc}")

    emit("正在获取 UniProt / iPTMnet 数据...")
    try:
        entry = fetch_uniprot(acc)
    except Exception as e:
        raise AnalysisError(f"UniProt 获取失败: {e}")

    iptmnet_data = fetch_iptmnet(acc)
    emit(f"iPTMnet 位点: {len(iptmnet_data)} 个")
    sequence = entry["sequence"]["value"]
    topology = build_topology_map(entry)
    emit(f"序列长度: {len(sequence)} | 拓扑注释区域: "
         f"{sum(1 for v in topology.values() if v != 'unknown')}")

    emit("正在提取已查证位点并核验文献...")
    verified = extract_verified(entry, iptmnet_data, topology, progress=emit)

    predictions = []
    for key, label in MODELS.items():
        emit(f"正在预测 {label}...")
        preds = predict_ptms(key, sequence, topology, entry)
        predictions.extend(preds)

    emit("正在计算跨物种保守性...")
    try:
        gene = entry["genes"][0]["geneName"]["value"]
    except Exception:
        gene = record.get("gene", "")
    predictions = compute_conservation(predictions, sequence, {"gene": gene})

    emit("正在评分排序...")
    predictions = score_predictions(predictions, verified)

    # Remove predictions that overlap with verified/homology
    known = {(v["ptm_type"], v["position"]) for v in verified}
    predictions = [p for p in predictions if (p["ptm_type"], p["position"]) not in known]

    for p in predictions:
        p["layer"] = "Predicted"

    all_results = verified + predictions
    all_results.sort(key=lambda x: ({"Verified": 0, "Supported": 1, "Predicted": 2}.get(x.get("layer", "Predicted"), 2), -x.get("score", 0)))

    n1 = sum(1 for r in all_results if r["layer"] == "Verified")
    n2 = sum(1 for r in all_results if r["layer"] == "Supported")
    n3 = sum(1 for r in all_results if r["layer"] == "Predicted")
    emit(f"完成: L1已查证 {n1} / L2有支持 {n2} / L3预测 {n3}")
    return all_results, record, entry, steps


def run(gpcr_input, refresh=False, output=None, verbose=False):
    try:
        all_results, record, entry, _steps = analyze(gpcr_input, refresh=refresh)
    except AnalysisError as e:
        print(f"[!] {e}")
        sys.exit(1)

    print_report(all_results, record, entry, verbose=verbose)
    save_json(all_results, record, output)
    write_verification_md(all_results, record)
    render_html(all_results, record, entry)


def main():
    parser = argparse.ArgumentParser(description="GPCR PTM predictor and verifier")
    parser.add_argument("gpcr", help="Gene name, UniProt ID, or common name (e.g., ADRB2, P07550, beta-2 adrenergic receptor)")
    parser.add_argument("--output", "-o", help="Output JSON file path")
    parser.add_argument("--refresh", action="store_true", help="Force re-download of offline database files")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="在终端同时打印每条 PMID 的链接与摘要证据(默认仅打印紧凑摘要)")
    args = parser.parse_args()
    run(args.gpcr, refresh=args.refresh, output=args.output, verbose=args.verbose)


if __name__ == "__main__":
    main()
