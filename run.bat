@echo off
REM GPCR-PTM 一键运行 (Windows 原生, 无需 Git-Bash / WSL)
REM 用法: 双击 run.bat 或命令行执行 run.bat
setlocal
cd /d "%~dp0"

REM 1) 选一个可用的 python
set "PY=python"
where py >nul 2>nul
if not errorlevel 1 set "PY=py -3"

REM 2) 创建虚拟环境 (若已有则复用)
if not exist "venv\Scripts\python.exe" (
  echo [*] First run: creating virtualenv venv ...
  %PY% -m venv venv
)
set "PY=venv\Scripts\python.exe"

REM 3) 确保 venv 里有 pip
"%PY%" -m pip --version >nul 2>nul
if errorlevel 1 (
  echo [*] pip not found, trying ensurepip ...
  "%PY%" -m ensurepip --upgrade
  if errorlevel 1 (
    echo [!] ensurepip failed. Install Python with "Add to PATH" and pip checked, then retry.
    pause
    exit /b 1
  )
)

REM 4) 安装依赖 (缺才装)
"%PY%" -c "import flask, requests" >nul 2>nul
if errorlevel 1 (
  echo [*] Installing dependencies (flask, requests) ...
  "%PY%" -m pip install -r requirements.txt
)

REM 5) 启动网页服务
echo [*] Starting web server ...
"%PY%" webapp.py %*
endlocal
