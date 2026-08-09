import json

LAYER_LABELS = {
    "Verified": "L1 已查证·实验证据",
    "Supported": "L2 有支持·未定论",
    "Predicted": "L3 预测·未被查证",
}


def _ver_summary(r):
    """紧凑核验汇总: '✓2 △3 ✗0' ; 有✗时附PMID"""
    ver = r.get("pmid_verification", [])
    if not ver:
        return ""
    d = {"直接": 0, "间接": 0, "不支持": 0, "其他": 0}
    bad = []
    for v in ver:
        vd = v["verdict"]
        if vd == "直接":
            d["直接"] += 1
        elif vd.startswith("间接"):
            d["间接"] += 1
        elif vd.startswith("不支持"):
            d["不支持"] += 1
            bad.append(v["pmid"])
        else:
            d["其他"] += 1
    s = f"✓{d['直接']} △{d['间接']} ✗{d['不支持']}"
    if bad:
        s += f"[✗:{','.join(bad)}]"
    return s


def print_report(results, record, entry, verbose=False):
    print()
    print("=" * 100)
    print(f"GPCR PTM Report: {record['gene']} ({record['name']}) | {record['accession']}")
    print("=" * 100)

    from collections import defaultdict
    layers = defaultdict(list)
    for r in results:
        layers[r["layer"]].append(r)

    for layer in ["Verified", "Supported", "Predicted"]:
        items = layers.get(layer, [])
        if not items:
            continue
        print(f"\n--- {LAYER_LABELS.get(layer, layer)} ({len(items)}) ---")
        if layer == "Predicted":
            print(f"{'PTM':<18} {'Pos':>5} {'Res':<5} {'Region':<12} {'Score':>6} {'Conf':<8} {'Motif':<20} Context")
            print("-" * 100)
            for r in items:
                print(f"{r['ptm_type']:<18} {r['position']:>5} {r['residue']:<5} {r['region']:<12} "
                      f"{r.get('score', '-'):>6} {r.get('confidence', '-'):<8} {r.get('motif', '-'):<20} {r.get('context', '')}")
        else:
            print(f"{'PTM':<18} {'Pos':>5} {'Res':<5} {'Region':<12} {'Score':>6} {'核验':<14} Evidence")
            print("-" * 100)
            for r in items:
                detail = r.get("layer_detail", "") or ""
                print(f"{r['ptm_type']:<18} {r['position']:>5} {r['residue']:<5} {r['region']:<12} "
                      f"{r.get('score', '-'):>6} {_ver_summary(r):<14} {detail}")
                if verbose:
                    ver = r.get("pmid_verification", [])
                    if ver:
                        for v in ver:
                            ev = (v.get("evidence") or "").replace("\n", " ")
                            print(f"    · PMID {v['pmid']}  [{v['verdict']}]  {v['url']}")
                            print(f"      {ev}")
                    else:
                        print("    · (无 PMID, 高通量单一检出, 需实验确认)")

    # Sequence overview
    sequence = entry.get("sequence", {}).get("value", "")
    if sequence:
        print("\n--- Sequence Overview (score>=0.5 marked) ---")
        marks = {}
        for r in results:
            if r.get("score", 0) >= 0.5:
                marks[r["position"]] = r["ptm_type"][0].upper()
        line = ""
        for i, aa in enumerate(sequence, 1):
            m = marks.get(i, "")
            line += f"{aa}{m} "
            if i % 20 == 0:
                print(f"{i-20:>4}: {line}")
                line = ""
        if line:
            print(f"{len(sequence)-len(line.split())+1:>4}: {line}")
        print()


def save_json(results, record, path=None):
    if path is None:
        path = f"{record['accession']}_ptm.json"
    output = {
        "accession": record["accession"],
        "gene": record["gene"],
        "name": record["name"],
        "layers": LAYER_LABELS,
        "results": results,
    }
    with open(path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"[*] JSON saved to {path}")


def write_verification_md(results, record, path=None):
    """把带链接/摘要证据/判定的核验表写入 Markdown 文件。"""
    if path is None:
        path = f"{record['accession']}_verification.md"

    lines = []
    lines.append(f"# {record['accession']} ({record['gene']}) PTM 位点文献核验表")
    lines.append("")
    lines.append("判定: **直接**=摘要点名该残基且含该PTM关键词; **间接**=该蛋白该PTM相关但未点名残基; "
                 "**不支持**=摘要点名同位置其他残基或其他蛋白。")
    lines.append("")
    lines.append("| 位点 | PTM | 层 | PMID | 链接 | 判定 | 摘要证据 |")
    lines.append("|---|---|---|---|---|---|---|")

    def esc(s):
        return s.replace("|", "/").replace("\n", " ")

    for r in results:
        ver = r.get("pmid_verification", [])
        pos = r["position"]
        if ver:
            for v in ver:
                lines.append(
                    f"| {pos} | {esc(r['ptm_type'])} | {r['layer']} | {v['pmid']} "
                    f"| [pubmed/{v['pmid']}]({v['url']}) | {v['verdict']} | {esc(v['evidence'])} |")
        elif r.get("layer") == "Predicted":
            refs_short = "; ".join(x["short"] for x in r.get("refs", []))
            lines.append(
                f"| {pos} | {esc(r['ptm_type'])} | Predicted | (无PMID) | — | 规则预测 "
                f"| motif={r.get('motif', '')}, 非数据库记录, 规则文献: {refs_short} |")
        else:
            lines.append(
                f"| {pos} | {esc(r['ptm_type'])} | {r['layer']} | (无PMID) | — | 仅高通量 "
                f"| iPTMnet/PSP 单一检出, 无文献可核验 |")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[*] Verification table saved to {path}")


# ---------------------------------------------------------------------------
# HTML 报告 (自包含, 浏览器双击打开即可阅读)
# ---------------------------------------------------------------------------
_CSS = """
:root{
  --bg:#e9edf5; --card:#ffffff; --ink:#1c2230; --muted:#5b6577; --line:#e2e7f0;
  --blue:#2b5cf0; --violet:#7c3aed; --green:#16a34a; --amber:#d97706; --red:#dc2626;
  --shadow:0 1px 2px rgba(16,24,40,.06),0 4px 12px rgba(16,24,40,.06);
  --shadow-lg:0 8px 30px rgba(16,24,40,.14);
}
*{box-sizing:border-box}
body{margin:0;color:var(--ink);
     background:
       radial-gradient(1100px 540px at 50% -180px,#dfe6f5 0%,rgba(223,230,245,0) 60%),
       linear-gradient(180deg,#e9edf5 0%,#e4e9f2 100%);
     background-attachment:fixed;
     font:15px/1.6 "Segoe UI","PingFang SC","Microsoft YaHei",system-ui,sans-serif;}
.wrap{max-width:1080px;margin:0 auto;padding:28px 20px 60px;}
/* ---------- top search bar (web) ---------- */
.topbar{display:flex;align-items:center;gap:16px;background:var(--card);
        border:1px solid var(--line);border-radius:14px;padding:12px 14px;
        box-shadow:var(--shadow);margin-bottom:18px;}
.topbar .qform{flex:1;display:flex;gap:10px;}
.topbar input{flex:1;border:1px solid var(--line);border-radius:10px;padding:10px 14px;
        font-size:14px;outline:none;}
.topbar input:focus{border-color:var(--blue);box-shadow:0 0 0 3px rgba(43,92,240,.12);}
.topbar button{background:linear-gradient(135deg,#2563eb,#4f46e5);color:#fff;border:none;
        border-radius:10px;padding:10px 22px;font-size:14px;font-weight:600;cursor:pointer;}
.topbar button:hover{opacity:.92;}
.topbar .home{color:var(--blue);font-weight:700;text-decoration:none;white-space:nowrap;}
/* ---------- header ---------- */
header{background:linear-gradient(135deg,#1d4ed8 0%,#4f46e5 55%,#7c3aed 100%);
       color:#fff;border-radius:18px;padding:30px 34px 24px;
       box-shadow:var(--shadow-lg);position:relative;overflow:hidden;}
header::after{content:"";position:absolute;right:-80px;top:-80px;width:260px;height:260px;
       border-radius:50%;background:rgba(255,255,255,.08);}
.kicker{font-size:12px;letter-spacing:3px;text-transform:uppercase;opacity:.85;}
h1{margin:6px 0 2px;font-size:26px;font-weight:700;}
h1 .acc{background:rgba(255,255,255,.18);padding:2px 12px;border-radius:999px;
        font-size:15px;font-weight:600;vertical-align:3px;margin-left:10px;}
.meta{font-size:13px;opacity:.9;}
.stats{display:flex;gap:12px;flex-wrap:wrap;margin-top:20px;}
.stat{background:rgba(255,255,255,.12);border:1px solid rgba(255,255,255,.18);
      border-radius:12px;padding:10px 16px;min-width:110px;backdrop-filter:blur(2px);}
.stat b{display:block;font-size:22px;line-height:1.2;}
.stat small{font-size:12px;opacity:.9;}
.stat.bad{border-color:#fecaca;background:rgba(239,68,68,.25);}
/* ---------- sections ---------- */
section{margin-top:34px;}
h2{font-size:18px;font-weight:700;display:flex;align-items:center;gap:10px;
   margin:0 0 14px;padding-bottom:8px;border-bottom:2px solid var(--line);}
h2 .bar{width:5px;height:18px;border-radius:3px;display:inline-block;}
h2 .n{color:var(--muted);font-weight:600;font-size:13px;}
h2.l1 .bar{background:var(--green);} h2.l2 .bar{background:var(--amber);} h2.l3 .bar{background:var(--violet);}
/* ---------- site card ---------- */
.site{background:var(--card);border:1px solid var(--line);border-radius:14px;
      box-shadow:var(--shadow);margin:14px 0;overflow:hidden;}
.site-row{display:flex;align-items:center;flex-wrap:wrap;gap:10px;padding:14px 18px;}
.site details{display:block;}
.site summary.site-row{cursor:pointer;list-style:none;}
.site summary.site-row::-webkit-details-marker{display:none;}
.site summary.site-row::after{content:"▾";margin-left:8px;font-size:13px;color:var(--muted);
        transition:transform .18s ease;}
.site details[open] summary.site-row::after{transform:rotate(180deg);}
.site details[open] summary.site-row{border-bottom:1px solid var(--line);}
.site summary.site-row:hover{background:#f8faff;}
.site .hint{margin-left:auto;font-size:11.5px;color:#9aa3b2;font-weight:600;letter-spacing:.03em;}
.site .detail{padding:12px 18px 4px;}
.site .pred{margin:0;padding:12px 18px;}
.res-badge{width:38px;height:38px;border-radius:10px;display:flex;align-items:center;
           justify-content:center;font-weight:800;font-size:18px;color:#fff;
           background:var(--blue);box-shadow:var(--shadow);}
.res-badge.phos{background:var(--red);} .res-badge.glyco{background:var(--blue);}
.res-badge.palm{background:var(--amber);} .res-badge.ubi{background:var(--violet);}
.pos{font-size:17px;font-weight:800;}
.ptm-tag{color:#334155;font-weight:600;font-size:13px;background:#eef2f7;
         padding:3px 10px;border-radius:999px;}
.region-pill{font-size:12px;color:#475467;background:#f4f6fb;border:1px solid var(--line);
             padding:2px 9px;border-radius:999px;}
.score-pill{font-size:12px;font-weight:700;color:#1d4ed8;background:#e8efff;
            padding:2px 9px;border-radius:999px;}
.src-pill{font-size:12px;color:#6b7280;background:#f9fafb;border:1px solid var(--line);
          padding:2px 9px;border-radius:999px;}
.detail{font-size:12.5px;color:var(--muted);padding:6px 18px 0;}
.motif{color:var(--violet);}
.ver{width:100%;border-collapse:collapse;font-size:13.5px;}
.ver th{text-align:left;color:var(--muted);font-size:11.5px;text-transform:uppercase;
        letter-spacing:.06em;background:#fafbfe;padding:8px 18px;border-bottom:1px solid var(--line);}
.ver td{padding:9px 18px;border-bottom:1px solid #f0f2f7;vertical-align:top;}
.ver tr:last-child td{border-bottom:none;}
.ver tr:hover td{background:#f8faff;}
.ver .ev{color:#4b5563;font-size:13px;}
.ver a{color:var(--blue);text-decoration:none;font-weight:600;}
.ver a:hover{text-decoration:underline;}
.pill{display:inline-block;border-radius:999px;padding:2px 11px;font-size:12px;font-weight:700;white-space:nowrap;}
.pill-ok{background:#e8f7ee;color:#15803d;}
.pill-mid{background:#fdf0d9;color:#b45309;}
.pill-bad{background:#fdeaea;color:#b91c1c;}
.pred{padding:10px 18px;font-size:14px;color:#374151;line-height:1.7;}
.pred b{color:var(--ink);}
.pred a{color:var(--blue);text-decoration:none;font-weight:600;margin-right:6px;}
.pred a:hover{text-decoration:underline;}
ul.reasons{margin:8px 0 8px;padding-left:20px;}
ul.reasons li{margin:5px 0;color:#3f4654;line-height:1.65;}
ul.reasons li::marker{color:var(--blue);}
/* ---------- sequence ---------- */
pre.seq{background:#0f172a;color:#cbd5e1;border-radius:14px;padding:16px 18px;
        font-size:13px;line-height:1.75;overflow-x:auto;font-family:ui-monospace,Consolas,monospace;}
pre.seq .ln{color:#64748b;user-select:none;margin-right:12px;}
.mk{color:#f8fafc;font-weight:800;padding:0 1px;border-radius:3px;}
.mk.phos{background:#dc2626;} .mk.glyco{background:#2563eb;} .mk.palm{background:#d97706;}
.mk.ubi{background:#7c3aed;} .mk.other{background:#64748b;}
/* ---------- legend / footer ---------- */
.legend{display:flex;gap:8px;flex-wrap:wrap;margin:16px 0 4px;}
.legend .pill{margin-right:4px;}
footer{margin-top:40px;padding-top:16px;border-top:1px solid var(--line);
       color:var(--muted);font-size:12.5px;line-height:1.9;}
@media(max-width:640px){ .wrap{padding:16px 12px;} h1{font-size:21px;} }
"""


def build_html(results, record, entry, with_search=False):
    """构建自包含 HTML 报告字符串。with_search=True 时顶部带继续查询框(网页用)。"""

    def h(s):
        return (str(s).replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;").replace('"', "&quot;"))

    from collections import Counter, defaultdict
    cnt = Counter(r["layer"] for r in results)
    n_pm = sum(len(r.get("pmid_verification", [])) for r in results)
    n_bad = sum(1 for r in results for v in r.get("pmid_verification", [])
                if v["verdict"].startswith("不支持"))

    body = []

    # ---------- top search bar (web) ----------
    if with_search:
        body.append("""
    <div class="topbar">
      <form action="/" method="get" class="qform">
        <input type="text" name="q" autofocus
               placeholder="输入基因名 / UniProt ID / 蛋白名称，如 ADRB2、P07550、beta-2 adrenergic receptor">
        <button type="submit">查询</button>
      </form>
      <a class="home" href="/">GPCR-PTM</a>
    </div>""")

    # ---------- header ----------
    body.append(f"""
    <header>
      <div class="kicker">GPCR Post-Translational Modification 报告</div>
      <h1>{h(record['gene'])} <span class="acc">{h(record['accession'])}</span></h1>
      <div class="meta">{h(record['name'])} · 序列 {entry.get('sequence', {}).get('length', 0)} aa</div>
      <div class="stats">
        <div class="stat"><small>L1 已查证</small><b>{cnt.get('Verified', 0)}</b></div>
        <div class="stat"><small>L2 有支持</small><b>{cnt.get('Supported', 0)}</b></div>
        <div class="stat"><small>L3 预测</small><b>{cnt.get('Predicted', 0)}</b></div>
        <div class="stat"><small>核验文献</small><b>{n_pm}</b></div>
        <div class="stat bad"><small>异常归属</small><b>{n_bad}</b></div>
      </div>
    </header>""")

    # ---------- legend ----------
    body.append("""
    <div class="legend">
      <span class="pill pill-ok">直接</span>摘要点名该残基且含该PTM关键词
      <span class="pill pill-mid">间接</span>该蛋白该PTM相关但未点名残基
      <span class="pill pill-bad">不支持</span>摘要点名同位置其他残基或其他蛋白
    </div>""")

    # ---------- layers ----------
    layer_meta = {
        "Verified": ("l1", "已查证 · 实验证据", "experimental"),
        "Supported": ("l2", "有支持 · 未定论", "supported"),
        "Predicted": ("l3", "预测 · 未被查证", "predicted"),
    }
    res2tag = {"phosphorylation": "phos", "glycosylation": "glyco",
               "palmitoylation": "palm", "ubiquitination": "ubi"}

    for layer in ["Verified", "Supported", "Predicted"]:
        items = [r for r in results if r["layer"] == layer]
        if not items:
            continue
        cls, title, _ = layer_meta[layer]
        body.append(f'<section><h2 class="{cls}"><span class="bar"></span>{title}'
                    f' <span class="n">({len(items)} 个位点)</span></h2>')

        for r in items:
            ver = r.get("pmid_verification", [])
            tag = res2tag.get(r["ptm_type"], "other")
            conf = r.get("confidence") or r.get("source") or ""
            extra = f'<span class="motif">motif: {h(r["motif"])}</span>' if r.get("motif") else ""

            hint = (f"文献证据 {len(ver)} 条" if ver
                    else ("预测依据" if r.get("layer") == "Predicted" else "说明"))
            body.append(f"""
            <div class="site">
              <details>
                <summary class="site-row">
                  <span class="res-badge {tag}">{h(r['residue'])}</span>
                  <span class="pos">{r['position']}</span>
                  <span class="ptm-tag">{h(r['ptm_type'])}</span>
                  <span class="region-pill">{h(r['region'])}</span>
                  <span class="score-pill">得分 {r.get('score', '-')}</span>
                  <span class="src-pill">{h(conf)}</span>
                  <span class="hint">{hint}</span>
                </summary>
                <div class="detail">{h(r.get('layer_detail', ''))} {extra}</div>""")

            if ver:
                body.append("""
                <table class="ver">
                  <tr><th>PMID</th><th>链接</th><th>判定</th><th>摘要证据</th></tr>""")
                for v in ver:
                    vd = v["verdict"]
                    pill = "pill-ok" if vd == "直接" else ("pill-bad" if vd.startswith("不支持") else "pill-mid")
                    body.append(
                        f'<tr><td style="font-weight:700">{h(v["pmid"])}</td>'
                        f'<td><a href="{h(v["url"])}" target="_blank">pubmed/{h(v["pmid"])}</a></td>'
                        f'<td><span class="pill {pill}">{h(vd)}</span></td>'
                        f'<td class="ev">{h(v["evidence"])}</td></tr>')
                body.append("</table>")
            elif r.get("layer") == "Predicted":
                rlist = "".join(f"<li>{h(x)}</li>" for x in r.get("reasons", []))
                body.append(f'<div class="pred"><b>规则预测</b>'
                            f'<ul class="reasons">{rlist}</ul></div>')
                refs = r.get("refs", [])
                if refs:
                    rlinks = " ".join(
                        f'<a href="{h(x["url"])}" target="_blank" title="{h(x["title"])}">{h(x["short"])}</a>'
                        for x in refs)
                    body.append(f'<div class="pred"><b>文献依据</b>（可点击核验）: {rlinks}</div>')
            else:
                body.append('<div class="pred">无 PMID，高通量单一检出，需实验确认。</div>')
            body.append("</details></div>")
        body.append("</section>")

    # ---------- sequence ----------
    sequence = entry.get("sequence", {}).get("value", "")
    if sequence:
        marks = {}
        for r in results:
            if r.get("score", 0) >= 0.5:
                marks[r["position"]] = res2tag.get(r["ptm_type"], "other")
        rows = []
        line = ""
        start = 1
        for i, aa in enumerate(sequence, 1):
            mk = marks.get(i)
            if mk:
                line += f'<span class="mk {mk}">{aa}</span> '
            else:
                line += f"{aa} "
            if i % 30 == 0:
                rows.append(f'<span class="ln">{start:>4}</span> {line}')
                line = ""
                start = i + 1
        if line:
            rows.append(f'<span class="ln">{start:>4}</span> {line}')
        body.append(f"""<section><h2 class="l3"><span class="bar"></span>序列标注
                    <span class="n">(score ≥ 0.5 的位点高亮)</span></h2>
                    <pre class="seq">""" + "\n".join(rows) + "</pre></section>")

    # ---------- footer ----------
    body.append(f"""
    <footer>
      说明：L1/L2 位点来自 UniProt、iPTMnet、dbPTM 等数据库；每条 PMID 均抓取 PubMed 摘要自动核验。
      预测位点(L3)基于共识 motif + 膜拓扑 + 跨物种保守性，未经实验验证，仅用于排序参考。
      数据抓取时间见运行日志。本报告为离线生成。
    </footer>""")

    return f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{h(record['accession'])} ({h(record['gene'])}) GPCR PTM 报告</title>
<style>{_CSS}</style>
</head>
<body>
<div class="wrap">
{''.join(body)}
</div>
</body>
</html>
"""


def render_html(results, record, entry, path=None):
    """生成 HTML 报告并写入文件 (CLI 用)。"""
    if path is None:
        path = f"{record['accession']}_report.html"
    html = build_html(results, record, entry, with_search=False)
    with open(path, "w") as f:
        f.write(html)
    print(f"[*] HTML report saved to {path}  (浏览器双击打开)")
