#!/usr/bin/env python3
"""
GPCR-PTM 网页服务 (Flask)

用法:
  ./venv/bin/python webapp.py                 # 默认 http://127.0.0.1:8000
  ./venv/bin/python webapp.py --port 9000 --host 0.0.0.0   # 局域网可访问
"""
import argparse
import threading
import uuid

from flask import Flask, request, jsonify, Response, abort

from main import analyze, AnalysisError
from output import build_html
from resolve import resolve_input

app = Flask(__name__)

JOBS = {}   # job_id -> {state, steps, html, error}
CACHE = {}  # accession -> html (已生成报告缓存)


def _run_job(job_id, q, record):
    j = JOBS[job_id]
    try:
        def progress(msg):
            j["steps"] = j["steps"] + [msg]

        results, rec, entry, _steps = analyze(q, progress=progress, record=record)
        progress("正在生成报告…")
        html = build_html(results, rec, entry, with_search=True)
        CACHE[rec["accession"]] = html
        j["state"] = "done"
        j["html"] = html
        j["steps"].append("完成")
    except AnalysisError as e:
        j["state"] = "error"
        j["error"] = str(e)
    except Exception as e:
        j["state"] = "error"
        j["error"] = f"分析失败: {e}"


@app.route("/")
def index():
    return _LANDING


@app.route("/api/submit", methods=["POST"])
def submit():
    data = request.get_json(silent=True) or {}
    q = (data.get("q") or "").strip()
    if not q:
        return jsonify({"error": "请输入基因名 / UniProt ID / 蛋白名称"}), 400
    if len(q) > 200:
        return jsonify({"error": "输入过长(≤200字符)"}), 400

    try:
        record = resolve_input(q)
    except Exception:
        record = None
    if record is None:
        return jsonify({"error": f"未找到 '{q}'，请检查基因名 / UniProt ID / 蛋白名称"}), 404

    job_id = uuid.uuid4().hex
    if record["accession"] in CACHE:
        JOBS[job_id] = {"state": "done", "steps": ["命中缓存，直接返回"],
                        "html": CACHE[record["accession"]], "error": None}
    else:
        JOBS[job_id] = {"state": "running", "steps": [f"已提交查询: {q}"],
                        "html": None, "error": None}
        threading.Thread(target=_run_job, args=(job_id, q, record), daemon=True).start()
    return jsonify({"job": job_id})


@app.route("/api/status")
def status():
    job = JOBS.get(request.args.get("job", ""))
    if not job:
        return jsonify({"error": "任务不存在或已过期"}), 404
    return jsonify({"state": job["state"], "steps": job["steps"],
                    "error": job.get("error")})


@app.route("/api/result")
def result():
    job = JOBS.get(request.args.get("job", ""))
    if not job or job["state"] != "done" or not job.get("html"):
        return abort(404)
    return Response(job["html"], mimetype="text/html")


_LANDING = """<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GPCR-PTM 查询</title>
<style>
  :root{--blue:#2563eb;--violet:#7c3aed;--line:#e6eaf2;--ink:#1c2230;--muted:#667085;}
  *{box-sizing:border-box}
  body{margin:0;background:#f2f4f8;color:var(--ink);
       font:15px/1.6 "Segoe UI","PingFang SC","Microsoft YaHei",system-ui,sans-serif;}
  .wrap{max-width:760px;margin:0 auto;padding:48px 20px;}
  header{background:linear-gradient(135deg,#1d4ed8,#4f46e5 55%,#7c3aed);
         color:#fff;border-radius:18px;padding:34px 36px;box-shadow:0 8px 30px rgba(16,24,40,.16);}
  .kicker{font-size:12px;letter-spacing:3px;text-transform:uppercase;opacity:.85;}
  h1{margin:6px 0 4px;font-size:26px;}
  .sub{font-size:13px;opacity:.9;}
  .card{background:#fff;border:1px solid var(--line);border-radius:16px;
        padding:22px;margin-top:20px;box-shadow:0 1px 3px rgba(16,24,40,.06);}
  form{display:flex;gap:10px;}
  input{flex:1;border:1px solid var(--line);border-radius:10px;padding:13px 16px;
        font-size:15px;outline:none;}
  input:focus{border-color:var(--blue);box-shadow:0 0 0 3px rgba(37,99,235,.12);}
  button{background:linear-gradient(135deg,#2563eb,#4f46e5);color:#fff;border:none;
         border-radius:10px;padding:13px 28px;font-size:15px;font-weight:600;cursor:pointer;}
  button:hover{opacity:.92;}
  .examples{margin-top:14px;font-size:12.5px;color:var(--muted);}
  .examples b{color:var(--blue);cursor:pointer;margin-right:12px;}
  .examples b:hover{text-decoration:underline;}
  #panel{background:#fff;border:1px solid var(--line);border-radius:16px;
         padding:22px;margin-top:20px;box-shadow:0 1px 3px rgba(16,24,40,.06);}
  .spinner{width:22px;height:22px;border:3px solid #dbe3f0;border-top-color:var(--blue);
           border-radius:50%;animation:spin .8s linear infinite;display:inline-block;vertical-align:-4px;}
  @keyframes spin{to{transform:rotate(360deg)}}
  #panel h3{margin:0 0 12px;font-size:15px;display:flex;align-items:center;gap:10px;}
  ol{margin:0;padding-left:22px;color:#374151;font-size:14px;}
  li{margin:3px 0;}
  li.done{color:#16a34a;} li.cur{color:var(--ink);font-weight:600;}
  .err{background:#fdeaea;border:1px solid #f5c2c2;color:#b91c1c;border-radius:10px;
       padding:12px 16px;margin-top:14px;font-size:14px;}
</style>
</head>
<body>
<div class="wrap">
  <header>
    <div class="kicker">GPCR Post-Translational Modification</div>
    <h1>GPCR-PTM 查询</h1>
    <div class="sub">输入基因名 / UniProt ID / 蛋白名称，输出已查证与预测的翻译后修饰位点</div>
  </header>
  <div class="card">
    <form id="qf">
      <input id="q" name="q" autofocus
             placeholder="例如 ADRB2、P07550、beta-2 adrenergic receptor">
      <button type="submit">查询</button>
    </form>
    <div class="examples">试试:
      <b onclick="setQ('ADRB2')">ADRB2</b>
      <b onclick="setQ('OPRM1')">OPRM1</b>
      <b onclick="setQ('DRD2')">DRD2</b>
    </div>
  </div>
  <div id="panel" hidden>
    <h3><span class="spinner"></span> 分析中…</h3>
    <ol id="steps"></ol>
    <div id="err" class="err" hidden></div>
  </div>
</div>
<script>
(function(){
  var params = new URLSearchParams(location.search);
  var q0 = params.get('q');
  document.getElementById('qf').addEventListener('submit', function(e){
    e.preventDefault();
    start(document.getElementById('q').value.trim());
  });
  window.setQ = function(v){ document.getElementById('q').value = v; start(v); };
  if (q0){ document.getElementById('q').value = q0; start(q0); }

  function start(q){
    if(!q) return;
    var panel = document.getElementById('panel'), err = document.getElementById('err');
    panel.hidden = false; err.hidden = true;
    var list = document.getElementById('steps');
    list.innerHTML = '<li class="cur">提交查询…</li>';
    fetch('/api/submit',{method:'POST',headers:{'Content-Type':'application/json'},
                         body:JSON.stringify({q:q})})
      .then(function(r){ return r.json().then(function(d){ return {ok:r.ok, d:d}; }); })
      .then(function(res){
        if(!res.ok || res.d.error){ showErr(res.d.error||'提交失败'); return; }
        poll(res.d.job);
      })
      .catch(function(){ showErr('网络请求失败，请重试'); });
  }

  function poll(id){
    var list = document.getElementById('steps');
    var t = setInterval(function(){
      fetch('/api/status?job='+id).then(function(r){ return r.json(); }).then(function(d){
        if(d.error){ clearInterval(t); showErr(d.error); return; }
        list.innerHTML = d.steps.map(function(s,i){
          var cls = (i === d.steps.length-1) ? 'cur' : 'done';
          return '<li class="'+cls+'">'+esc(s)+'</li>';
        }).join('');
        if(d.state==='done'){ clearInterval(t); location.href='/api/result?job='+id; }
        else if(d.state==='error'){ clearInterval(t); showErr(d.error||'分析失败'); }
      }).catch(function(){});
    }, 1200);
  }

  function showErr(m){
    document.getElementById('panel').hidden = true;
    var e = document.getElementById('err');
    e.textContent = m; e.hidden = false;
  }
  function esc(s){ var d=document.createElement('div'); d.textContent=s; return d.innerHTML; }
})();
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="GPCR-PTM 网页服务")
    parser.add_argument("--host", default="127.0.0.1", help="监听地址(默认 127.0.0.1；0.0.0.0 允许局域网访问)")
    parser.add_argument("--port", type=int, default=8000, help="监听端口(默认 8000)")
    args = parser.parse_args()
    print(f"[*] GPCR-PTM 网页服务已启动: http://{args.host}:{args.port}")
    print("[*] 浏览器打开上述地址开始查询; Ctrl+C 停止")
    app.run(host=args.host, port=args.port, threaded=True)


if __name__ == "__main__":
    main()
