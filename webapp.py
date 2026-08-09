#!/usr/bin/env python3
"""
GPCR-PTM 网页服务 (零依赖, 仅用 Python 标准库)

无需 pip install 任何 Web 框架。只要本机有 Python 3.8+ 即可一键启动。

用法:
  python webapp.py                              # 默认 http://127.0.0.1:8000
  python webapp.py --port 9000 --host 0.0.0.0   # 局域网可访问
  python webapp.py --no-browser                 # 不自动打开浏览器
"""
import argparse
import json
import threading
import time
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

from main import analyze, AnalysisError
from output import build_html
from resolve import resolve_input

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


# ---------- API 业务逻辑 (与 HTTP 层解耦, 便于测试) ----------
def api_submit(payload):
    """处理 /api/submit。返回 (status_code, dict)。"""
    q = ((payload or {}).get("q") or "").strip()
    if not q:
        return 400, {"error": "请输入基因名 / UniProt ID / 蛋白名称"}
    if len(q) > 200:
        return 400, {"error": "输入过长(≤200字符)"}

    try:
        record = resolve_input(q)
    except Exception:
        record = None
    if record is None:
        return 404, {"error": f"未找到 '{q}'，请检查基因名 / UniProt ID / 蛋白名称"}

    job_id = uuid.uuid4().hex
    if record["accession"] in CACHE:
        JOBS[job_id] = {"state": "done", "steps": ["命中缓存，直接返回"],
                        "html": CACHE[record["accession"]], "error": None}
    else:
        JOBS[job_id] = {"state": "running", "steps": [f"已提交查询: {q}"],
                        "html": None, "error": None}
        threading.Thread(target=_run_job, args=(job_id, q, record), daemon=True).start()
    return 200, {"job": job_id}


def api_status(job_id):
    """处理 /api/status。返回 (status_code, dict)。"""
    job = JOBS.get(job_id or "")
    if not job:
        return 404, {"error": "任务不存在或已过期"}
    return 200, {"state": job["state"], "steps": job["steps"],
                 "error": job.get("error")}


def api_result(job_id):
    """处理 /api/result。返回 (status_code, body_or_None, content_type)。"""
    job = JOBS.get(job_id or "")
    if not job or job["state"] != "done" or not job.get("html"):
        return 404, None, None
    return 200, job["html"], "text/html; charset=utf-8"


# ---------- HTTP 层: 标准 http.server ----------
class Handler(BaseHTTPRequestHandler):
    # 静默默认日志刷屏, 仅保留错误
    def log_message(self, fmt, *args):
        pass

    def _send_json(self, code, obj):
        body = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_html(self, code, html, content_type="text/html; charset=utf-8"):
        body = html.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urlparse(self.path)
        path, qs = parsed.path, parse_qs(parsed.query)
        if path == "/":
            self._send_html(200, _LANDING)
        elif path == "/api/status":
            code, obj = api_status(qs.get("job", [""])[0])
            self._send_json(code, obj)
        elif path == "/api/result":
            code, body, ctype = api_result(qs.get("job", [""])[0])
            if code == 200:
                self._send_html(code, body, ctype)
            else:
                self.send_error(code)
        else:
            self.send_error(404)

    def do_POST(self):
        parsed = urlparse(self.path)
        if parsed.path != "/api/submit":
            self.send_error(404)
            return
        try:
            length = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(length) if length else b"{}"
            payload = json.loads(raw.decode("utf-8") or "{}")
        except Exception:
            payload = {}
        code, obj = api_submit(payload)
        self._send_json(code, obj)


_LANDING = """<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GPCR-PTM 查询</title>
<style>
  :root{--blue:#2563eb;--violet:#7c3aed;--line:#e2e7f0;--ink:#1c2230;--muted:#5b6577;}
  *{box-sizing:border-box}
  body{margin:0;color:var(--ink);
       background:
         radial-gradient(1100px 540px at 50% -180px,#dfe6f5 0%,rgba(223,230,245,0) 60%),
         linear-gradient(180deg,#e9edf5 0%,#e4e9f2 100%);
       background-attachment:fixed;
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
  li{margin:3px 0;display:flex;align-items:center;gap:8px;}
  li.done{color:#16a34a;}
  li.cur{color:var(--ink);font-weight:600;}
  li .tick{width:16px;text-align:center;flex:0 0 auto;}
  li.done .tick{color:#16a34a;}
  li.cur .tick .mini{width:13px;height:13px;border:2px solid #dbe3f0;border-top-color:var(--blue);
        border-radius:50%;animation:spin .8s linear infinite;display:inline-block;vertical-align:-2px;}
  li.wait{color:var(--muted);}
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

  // ---- 逐行匀速揭示播放器 ----
  // 后端 steps 产出快慢不均(拉序列/调API/文献核验耗时差异大), 直接整批渲染会
  // 出现"先憋几秒再一次蹦多行". 这里把"后端到达"与"前端揭示"解耦:
  //   - 后端轮询只负责把 steps 攒进 buffer (生产者)
  //   - 节拍器每 STEP_MS 揭示一行 (消费者), 揭示到与后端齐平时就停在该行等待
  //   - done 后也要等剩余行按节拍播完再跳转, 保证节奏恒定
  var STEP_MS = 700;
  var buffer = [];        // 后端已到达的 steps 文本
  var state = 'running';  // running | done | error
  var errorMsg = null;
  var revealed = 0;       // 已揭示行数 (0..buffer.length)
  var pollTimer = null, tickTimer = null;
  var curJobId = null;    // 当前任务 id, 供 tick 跳转用

  function poll(id){
    curJobId = id;
    pollTimer = setInterval(function(){
      fetch('/api/status?job='+id).then(function(r){ return r.json(); }).then(function(d){
        if(d.error){ state='error'; errorMsg=d.error; stop(); showErr(d.error); return; }
        buffer = d.steps || buffer;
        if(d.state==='done'){ state='done'; stopPoll(); }
        else if(d.state==='error'){ state='error'; errorMsg=d.error||'分析失败'; stop(); showErr(errorMsg); }
      }).catch(function(){});
    }, 1000);

    tickTimer = setInterval(tick, STEP_MS);
    tick();   // 立即播首行
  }

  function stopPoll(){ if(pollTimer){ clearInterval(pollTimer); pollTimer=null; } }
  function stop(){ stopPoll(); /* tickTimer 由 tick 自己在播完后停 */ }

  function tick(){
    var list = document.getElementById('steps');
    // 还有未揭示的后端行 -> 揭示一行 (匀速)
    if(revealed < buffer.length){ revealed++; }
    render();

    if(state==='error'){ if(tickTimer){ clearInterval(tickTimer); tickTimer=null; } return; }
    // 已全部播完且后端也结束 -> 跳转结果页
    if(state==='done' && revealed >= buffer.length){
      if(tickTimer){ clearInterval(tickTimer); tickTimer=null; }
      // 让"完成"这行在屏幕上停一拍再跳, 避免最后一闪而过
      setTimeout(function(){ location.href='/api/result?job='+curJobId; }, STEP_MS);
    }
  }

  function render(){
    var list = document.getElementById('steps');
    var html = '';
    for(var i=0;i<revealed;i++){
      html += '<li class="done"><span class="tick">✓</span><span>'+esc(buffer[i])+'</span></li>';
    }
    // 等待中的"下一行": 后端还没给或正在给 -> 显示 spinner + 占位.
    // buffer 还空时显示"提交查询…", 之后统一"正在处理…" (动态行如 PMID 核验无法预知文案)
    var waiting = null;
    if(state!=='done'){ waiting = (revealed===0) ? '提交查询…' : '正在处理…'; }
    if(waiting){
      html += '<li class="cur"><span class="tick"><span class="mini"></span></span><span>'+esc(waiting)+'</span></li>';
    }
    list.innerHTML = html;
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
    parser = argparse.ArgumentParser(description="GPCR-PTM 网页服务 (零依赖)")
    parser.add_argument("--host", default="127.0.0.1", help="监听地址(默认 127.0.0.1；0.0.0.0 允许局域网访问)")
    parser.add_argument("--port", type=int, default=8000, help="监听端口(默认 8000)")
    parser.add_argument("--no-browser", action="store_true",
                        help="启动后不自动打开浏览器(默认自动打开)")
    args = parser.parse_args()
    print(f"[*] GPCR-PTM 网页服务已启动: http://{args.host}:{args.port}")
    print("[*] 零依赖模式 (Python 标准库 http.server)")
    print("[*] Ctrl+C 停止")
    if not args.no_browser:
        url = (f"http://127.0.0.1:{args.port}" if args.host in ("0.0.0.0", "::")
               else f"http://{args.host}:{args.port}")
        threading.Thread(target=lambda: (time.sleep(1.5), webbrowser.open(url)),
                         daemon=True).start()
    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n[*] 已停止")
        httpd.shutdown()


if __name__ == "__main__":
    main()
