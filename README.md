# GPCR-PTM

Predict and verify **post-translational modification (PTM) sites** on G-protein-coupled receptors (GPCRs).

Given a GPCR (gene name, UniProt accession, or protein name), GPCR-PTM outputs:

- **Verified PTM sites** — curated from UniProt / iPTMnet / dbPTM with a three-tier evidence grading
- **Predicted PTM sites** — rule-based candidates from consensus motifs + membrane topology + cross-species conservation
- **Per-PMID literature verification** — for every cited paper it fetches the PubMed abstract and auto-judges whether it actually supports *that residue + that PTM type* (verdict: **direct / indirect / unsupported**), with clickable links and abstract excerpts

An interactive **Flask web interface** (with async progress bar) and a **CLI** are both provided.

---

## Features

| Layer | Meaning |
|------|---------|
| **L1 Verified** | UniProt experimental evidence (`ECO:0000269`, with PMID) |
| **L2 Supported** | Inferred / by-similarity (UniProt `ECO:0000255/0305/7744`), iPTMnet / dbPTM literature or single high-throughput hits |
| **L3 Predicted** | Rule-based candidates (consensus motif + topology + conservation), **not** in any database |

PTM types covered:

- Phosphorylation (S/T — PKA, PKC, CK1, CK2, CDK/MAPK, proline-directed, GRK positional hint)
- N-glycosylation (N-X-S/T sequon, X≠P, extracellular only)
- Palmitoylation (juxtamembrane / K/R-rich Cys in the C-tail)
- Ubiquitination is reported **from databases only** (no reliable linear consensus for prediction)

Output formats:

- **Interactive HTML report** — collapsible site cards, color-coded verdict badges, clickable PubMed links, sequence overview with highlighted sites
- **JSON** — machine-readable, includes per-PMID verification data
- **Markdown** — the verification table with links + abstract evidence

---

## Quick Start

### Linux / macOS

```bash
bash run.sh            # first run installs everything automatically, then starts the web server
# open http://127.0.0.1:8000
```

`run.sh` automatically: creates a `venv` (using `--without-pip` to survive Debian/Ubuntu systems without `ensurepip`), bootstraps `pip` via `get-pip.py`, installs `requirements.txt`, and starts the server.

### Windows

Double-click **`run.bat`** (native, no Git-Bash / WSL required), then open http://127.0.0.1:8000.

### Manual install

```bash
python3 -m venv --without-pip venv
curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
venv/bin/python get-pip.py
venv/bin/pip install -r requirements.txt
venv/bin/python webapp.py            # web interface
# or
venv/bin/python main.py ADRB2        # CLI
```

Requirements: **Python 3.8+**, `requests`, `flask` (see `requirements.txt`).

---

## Usage

### Web interface

1. Start the server: `bash run.sh` (or `./venv/bin/python webapp.py --host 127.0.0.1 --port 8000`)
2. Open http://127.0.0.1:8000
3. Type a GPCR, e.g. `ADRB2`, `P07550`, `OPRM1`, or `beta-2 adrenergic receptor`
4. Watch the progress steps (first run verifies literature per PMID, ~15–30 s; results are then cached)
5. The collapsible report appears; you can search again from the top bar

### CLI

```bash
python main.py ADRB2                          # gene name
python main.py P07550                         # UniProt accession
python main.py "beta-2 adrenergic receptor"   # protein name
python main.py OPRM1 --verbose                # also print per-PMID links/evidence in the terminal
python main.py ADRB2 -o my_result.json        # custom JSON path
python main.py --refresh                      # (optional) download UniProt Swiss-Prot human flat file (~120 MB)
```

---

## How it works

```
Input (gene / accession / name)
  → resolve to UniProt entry
  → fetch UniProt + iPTMnet
  → build membrane topology (TM / ECL / ICL / N-term / C-tail)
  → extract verified PTMs (L1/L2) with evidence grading
      → per-PMID PubMed abstract verification (direct / indirect / unsupported)
  → rule-based prediction (L3): motif + topology filter
  → cross-species conservation (Needleman-Wunsch global alignment)
  → scoring & ranking
  → HTML / JSON / Markdown output
```

Key design decisions:

- **Topology-aware**: N-glycosylation is only considered in extracellular domains; phosphorylation / palmitoylation only in intracellular domains — this filters out the most common false positives.
- **Evidence grading**: UniProt experimental evidence (`ECO:0000269`) is ranked above inferred or high-throughput hits; iPTMnet sites without PMIDs are flagged "single high-throughput hit, needs experimental confirmation".
- **Automated literature verification**: for each PMID, the abstract is fetched and checked for the specific residue (e.g. `Ser355`, `Ser(345,346)`, `tyrosine-141`) and the PTM type. This catches mis-attributions such as a **Tyr364** paper being attached to a **Ser364** site, or a vasopressin-receptor paper being attached to a β2-AR site.
- **Motif specificity**: no more "all S/T in the intracellular tail" noise — each prediction must match a real consensus motif.

---

## Data sources & caching

- **UniProt** REST API (sequence, topology, features, evidence codes)
- **iPTMnet** REST API (literature-curated PTM sites, with PMIDs)
- **dbPTM / iPTMnet-bulk** (optional offline flat files, see `update_ptm_db.py`)
- **PubMed** (abstracts for verification)
- **Offline UniProt Swiss-Prot human flat file** (optional, `python main.py --refresh`)

Caching:

- `data/pmid_abstract_cache.json` — PubMed title/abstract cache (repeated queries skip re-fetching)
- In-memory result cache in the web app — the same protein is served instantly on repeat queries

---

## Output files

For a query on `P07550`, three files are written to the working directory:

| File | Contents |
|------|----------|
| `P07550_report.html` | Interactive collapsible report (open in any browser) |
| `P07550_ptm.json` | Structured results incl. `pmid_verification` per site |
| `P07550_verification.md` | Markdown verification table with links + evidence |

---

## Project layout

```
main.py               CLI entry point + analysis pipeline (analyze())
webapp.py             Flask web interface (async jobs + progress polling)
output.py             console report, JSON, Markdown, HTML (build_html)
verify_pmids.py       per-PMID PubMed verification (verdict + evidence)
verified.py           L1/L2 extraction & evidence grading
predict.py            rule-based L3 prediction (motifs)
conservation.py       cross-species conservation (Needleman-Wunsch)
topology.py           membrane topology map
scoring.py            prediction confidence scoring
resolve.py            input → UniProt resolution
fetch.py              UniProt / iPTMnet fetching
ptm_sources.py        iPTMnet / dbPTM parsers
local_db.py           offline Swiss-Prot flat-file parser
download_data.py      offline UniProt flat-file download
update_ptm_db.py      optional dbPTM / iPTMnet-bulk download
run.sh / run.bat      one-command setup + start
requirements.txt      dependencies
data/                 caches & user-extensions (data/ptm_extensions.tsv)
```

---

## Notes & limitations

- **Consensus motifs are necessary but not sufficient.** Only ~30–40% of N-glycosylation sequons are actually occupied; phospho motifs are weak predictors. Scores are for **ranking only**, never probabilities.
- **Conservation** is computed against orthologs found by exact gene name (mostly mammals), so conserved regions saturate near 1.0.
- The **"unsupported" verdict** means the cited abstract does not support *that residue + that PTM type* — it does not necessarily mean the site is wrong, just that the cited evidence is questionable.
- Predictions are hypotheses and should be validated experimentally (e.g. by mass spectrometry).
- For research/educational use only; the data is drawn from public databases.

## License

Released under the [MIT License](LICENSE). Please cite the underlying data sources (UniProt, iPTMnet, dbPTM, PubMed) when using results.

---

---

# GPCR-PTM

预测并核验 **G 蛋白偶联受体（GPCR）的翻译后修饰（PTM）位点**。

输入一个 GPCR（基因名 / UniProt 编号 / 蛋白名），输出：

- **已查证的 PTM 位点** —— 来自 UniProt / iPTMnet / dbPTM，并带三级证据分级
- **预测的 PTM 位点** —— 基于共识 motif + 膜拓扑 + 跨物种保守性的规则预测
- **逐条文献核验** —— 对引用的每一篇文献抓取 PubMed 摘要，自动判断它是否真的支持"该残基 + 该 PTM 类型"（判定：**直接 / 间接 / 不支持**），附可点击链接与摘要证据

同时提供 **Flask 网页界面**（带异步进度条）与**命令行（CLI）**两种用法。

---

## 功能特性

| 层级 | 含义 |
|------|------|
| **L1 已查证** | UniProt 实验证据（`ECO:0000269`，带 PMID） |
| **L2 有支持** | 推断 / 相似性（`ECO:0000255/0305/7744`）、iPTMnet / dbPTM 文献或高通量单一检出 |
| **L3 预测** | 规则预测候选（共识 motif + 拓扑 + 保守性），**未见于任何数据库** |

覆盖的 PTM 类型：

- 磷酸化（S/T：PKA、PKC、CK1、CK2、CDK/MAPK、proline-directed、GRK 位置提示）
- N-糖基化（N-X-S/T sequon，X≠P，仅胞外）
- 棕榈酰化（C 尾膜旁 / K/R 富集的半胱氨酸）
- 泛素化**仅从数据库报告**（无可靠线性 consensus，不做规则预测）

输出格式：

- **交互式 HTML 报告** —— 可折叠位点卡片、判定徽章配色、可点击 PubMed 链接、序列高亮标注
- **JSON** —— 机器可读，含每个位点的逐条 `pmid_verification` 数据
- **Markdown** —— 带链接与摘要证据的核验表

---

## 快速开始

### Linux / macOS

```bash
bash run.sh            # 首次运行自动装好环境并启动网页服务
# 浏览器打开 http://127.0.0.1:8000
```

`run.sh` 会自动：创建 `venv`（用 `--without-pip` 规避 Debian/Ubuntu 缺 `ensurepip` 导致失败）、用 `get-pip.py` 引导安装 `pip`、安装 `requirements.txt`、启动服务。

### Windows

直接双击 **`run.bat`**（原生脚本，无需 Git-Bash / WSL），然后打开 http://127.0.0.1:8000。

### 手动安装

```bash
python3 -m venv --without-pip venv
curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
venv/bin/python get-pip.py
venv/bin/pip install -r requirements.txt
venv/bin/python webapp.py            # 网页界面
# 或
venv/bin/python main.py ADRB2        # 命令行
```

依赖：**Python 3.8+**、`requests`、`flask`（见 `requirements.txt`）。

---

## 使用方法

### 网页界面

1. 启动服务：`bash run.sh`（或 `./venv/bin/python webapp.py --host 127.0.0.1 --port 8000`）
2. 打开 http://127.0.0.1:8000
3. 输入 GPCR，如 `ADRB2`、`P07550`、`OPRM1` 或 `beta-2 adrenergic receptor`
4. 查看进度步骤（首次会逐条核验文献摘要，约 15–30 秒；之后有缓存秒回）
5. 展示可折叠报告；顶部搜索框可继续查询

### 命令行

```bash
python main.py ADRB2                          # 基因名
python main.py P07550                         # UniProt 编号
python main.py "beta-2 adrenergic receptor"   # 蛋白名
python main.py OPRM1 --verbose                # 同时在终端打印每条 PMID 的链接与证据
python main.py ADRB2 -o my_result.json        # 自定义 JSON 输出路径
python main.py --refresh                      # （可选）下载 UniProt 人源 Swiss-Prot 离线文件（约 120 MB）
```

---

## 工作原理

```
输入（基因名 / 编号 / 蛋白名）
  → 解析为 UniProt 条目
  → 获取 UniProt + iPTMnet
  → 构建膜拓扑（TM / ECL / ICL / N端 / C尾）
  → 提取已查证 PTM（L1/L2）并进行证据分级
      → 逐条 PMID 抓取 PubMed 摘要核验（直接 / 间接 / 不支持）
  → 规则预测（L3）：motif + 拓扑过滤
  → 跨物种保守性（Needleman-Wunsch 全局比对）
  → 评分排序
  → HTML / JSON / Markdown 输出
```

关键设计：

- **拓扑感知**：N-糖基化只考虑胞外结构域；磷酸化 / 棕榈酰化只考虑胞内结构域——过滤掉最常见的假阳性。
- **证据分级**：UniProt 实验证据（`ECO:0000269`）高于推断或高通量检出；iPTMnet 无 PMID 的位点标注为"高通量单一检出，需实验确认"。
- **自动文献核验**：对每条 PMID 抓取摘要，检查是否点名该残基（如 `Ser355`、`Ser(345,346)`、`tyrosine-141`）与该 PTM 类型，能自动抓出误归属——例如把 **Tyr364** 的论文挂在 **Ser364** 位点上、或把加压素受体论文挂在 β2-AR 位点上。
- **motif 特异度**：不再输出"胞内所有 S/T"这类噪声，每个预测都必须命中真实共识 motif。

---

## 数据来源与缓存

- **UniProt** REST API（序列、拓扑、feature、证据代码）
- **iPTMnet** REST API（带 PMID 的文献位点）
- **dbPTM / iPTMnet-bulk**（可选的离线文件，见 `update_ptm_db.py`）
- **PubMed**（核验用摘要）
- **UniProt 人源 Swiss-Prot 离线文件**（可选，`python main.py --refresh`）

缓存：

- `data/pmid_abstract_cache.json` —— PubMed 标题/摘要缓存（重复查询不重复抓取）
- 网页端内存结果缓存 —— 同一蛋白重复查询直接命中

---

## 输出文件

以查询 `P07550` 为例，会在当前目录生成三个文件：

| 文件 | 内容 |
|------|------|
| `P07550_report.html` | 交互式可折叠报告（浏览器双击打开） |
| `P07550_ptm.json` | 结构化结果，含每位的 `pmid_verification` |
| `P07550_verification.md` | 带链接与摘要证据的核验表 |

---

## 项目结构

```
main.py               CLI 入口 + 分析流水线（analyze()）
webapp.py             Flask 网页界面（异步任务 + 进度轮询）
output.py             控制台报告、JSON、Markdown、HTML（build_html）
verify_pmids.py       逐条 PMID 的 PubMed 核验（判定 + 证据）
verified.py           L1/L2 提取与证据分级
predict.py            规则预测 L3（motif）
conservation.py       跨物种保守性（Needleman-Wunsch）
topology.py           膜拓扑图谱
scoring.py            预测置信度评分
resolve.py            输入 → UniProt 解析
fetch.py              UniProt / iPTMnet 获取
ptm_sources.py        iPTMnet / dbPTM 解析器
local_db.py           离线 Swiss-Prot 文件解析
download_data.py      离线 UniProt 文件下载
update_ptm_db.py      可选的 dbPTM / iPTMnet-bulk 下载
run.sh / run.bat      一键安装并启动
requirements.txt      依赖清单
data/                 缓存与用户扩展（data/ptm_extensions.tsv）
```

---

## 注意事项与局限

- **共识 motif 只是必要条件而非充分条件。** 大约只有 30–40% 的 N-糖基化 sequon 实际被占用；磷酸化 motif 预测能力较弱。分数**仅用于排序，不代表概率**。
- **保守性**基于按基因名匹配到的直系同源（主要是哺乳动物），保守区域会饱和到接近 1.0。
- **"不支持"判定**指所引摘要不支持"该残基 + 该 PTM"，并不代表位点一定错误，只说明所引证据存疑。
- 预测结果属于假设，建议通过实验（如质谱）验证。
- 数据来自公开数据库，供研究/教学使用。

## License

本项目采用 [MIT License](LICENSE) 开源。引用结果时请注明底层数据来源（UniProt、iPTMnet、dbPTM、PubMed）。
