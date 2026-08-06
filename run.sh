#!/usr/bin/env bash
# GPCR-PTM 一键运行 (Linux / macOS / Git-Bash / WSL)
# 用法:
#   bash run.sh                # 首次自动建虚拟环境+装依赖, 然后启动
#   bash run.sh --port 9000    # 透传 webapp 参数
set -e
cd "$(dirname "$0")"

# 选一个可用的 python
if command -v python3 >/dev/null 2>&1; then PY=python3; else PY=python; fi

# 1) 创建虚拟环境 (若已有则复用)
if [ ! -x venv/bin/python ]; then
  echo "[*] 首次运行: 创建虚拟环境 venv ..."
  # --without-pip: 某些系统(如 Debian/Ubuntu)缺 ensurepip, 会导致 -m venv 失败;
  # 用该选项确保创建成功, pip 随后由 get-pip.py 引导安装。
  if ! "$PY" -m venv --without-pip venv 2>/dev/null && [ ! -x venv/bin/python ]; then
    echo "[!] 虚拟环境创建失败。请先安装系统包: sudo apt install python3-venv"
    exit 1
  fi
fi
PY=venv/bin/python

# 2) 确保 venv 里有 pip
if ! "$PY" -m pip --version >/dev/null 2>&1; then
  echo "[*] venv 缺少 pip, 正在引导安装 ..."
  TMPDIRX="$(mktemp -d)"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL https://bootstrap.pypa.io/get-pip.py -o "$TMPDIRX/get-pip.py"
  elif command -v wget >/dev/null 2>&1; then
    wget -q https://bootstrap.pypa.io/get-pip.py -O "$TMPDIRX/get-pip.py"
  else
    "$PY" -c "import urllib.request; urllib.request.urlretrieve('https://bootstrap.pypa.io/get-pip.py', '$TMPDIRX/get-pip.py')"
  fi
  "$PY" "$TMPDIRX/get-pip.py"
  rm -rf "$TMPDIRX"
fi

# 3) 安装依赖 (缺才装)
if ! "$PY" -c "import flask, requests" >/dev/null 2>&1; then
  echo "[*] 安装依赖 (flask, requests) ..."
  "$PY" -m pip install -r requirements.txt
fi

# 4) 启动网页服务
echo "[*] 启动网页服务 ..."
exec "$PY" webapp.py "$@"
