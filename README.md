# GPCR-PTM

输入GPCR（基因名/UniProt编号/蛋白名），即可在线查询PTM位点，包括已查证有实验报道的位点，有支持但尚未定论的位点（推断或相似性注释、数据库或文献命中、高通量单一检出），以及基于规则提出的预测位点。

---

## 使用方法

### 安装依赖

```text
requests>=2.28
```

### 网页界面

网页服务基于 Python 标准库，无需额外安装 Web 框架。

```bash
python webapp.py                 # 默认 http://127.0.0.1:8000，启动后自动打开浏览器
python webapp.py --port 9000     # 自定义端口
python webapp.py --host 0.0.0.0  # 允许局域网访问
python webapp.py --no-browser    # 启动后不自动打开浏览器
```

### 命令行界面

```bash
python main.py ADRB2
python main.py P07550
python main.py OPRM1 --verbose
```

---

## 结果说明

### 预测流程

**基于规则的预测，不是机器学习。** 程序不训练数据、不学习参数，而是套用明确、成熟的生物学规则，因此每个结果都有可追溯的理由、完全可复现。候选位点分三步产生：

1. **motif扫描** —— 在序列中搜索公认的共识motif，每类 PTM 背后都有明确的生化机制：

   - **磷酸化** 扫描丝氨酸/苏氨酸上的激酶共识（PKA、PKC、CK1、CK2、CDK/MAPK 等经典 motif [1,2]）。GRK（GPCR 激酶）没有可靠的线性 motif [3,4]，仅以"邻近酸性残基 + S/T 簇集"作弱位置提示 [5,6]。此外还扫描 **GPCR 专属的 pXpp 三磷酸簇 `[S/T]-X-[S/T]-[S/T]`**——三个磷酸化 S/T 残基（p 代表磷酸化的 S/T，不是脯氨酸），是 arrestin 募集的关键模块 [7]。pXpp 在不同 GPCR 中的分布位置不同：趋化因子/肽类受体（class B）的 pXpp 富集于 C 端尾，胺类受体如多巴胺/肾上腺素能（class A）的 pXpp 则富集于较长的第三胞内环（ICL3）[7]。因此程序会先按受体的 arrestin 行为判组（依据 Isaikina 2023 的 Table S2 实验分类 [7]，26 个受体；未命中者按 ICL3 长度回退），**只在 pXpp 典型分布的区域报告，过滤原文标明 pXpp 缺失的区域**——例如 class A 受体的 C 端尾 pXpp 命中会被过滤，反之 class B 受体的 ICL3 pXpp 会被过滤。
   - **N-糖基化** 扫描 N-X-S/T sequon（X≠P）。糖基化的本质是膜表面 trafficking 的质量控——缺少糖基化的受体会滞留胞内、难以到达细胞膜，但糖基化位点在 N 端还是胞外环并不关键（任一处都能补救表面表达）[8]。
   - **棕榈酰化** 扫描胞内半胱氨酸（Cys），重点是 TM7 后膜旁区的 Cys——硫酯键把 Cys 锚定到质膜内侧，形成"第四胞内环"，稳定受体构象 [9]。
   - **泛素化**（保守版）扫描胞内赖氨酸（Lys）+ 邻近 PPxY 或 Pro-rich 区。泛素化决定受体内吞后的去向（降解或循环）[10]；Nedd4 家族 E3 连接酶通过其 WW 结构域识别底物的 PPxY (P-P-x-Y) motif 来完成泛素化"挂载" [10,11]。但 PPxY 多出现在 adaptor 蛋白（如 ARRDC3）上而非 GPCR 自身，因此此规则只能作保守弱提示，恒低置信。

2. **膜拓扑过滤** —— 只保留生化上正确的位置：N-糖基化必须在胞外（酶在内质网腔）；磷酸化/棕榈酰化/泛素化必须在胞内（激酶/E3 连接酶在胞浆）。这一步滤掉绝大多数假阳性。对磷酸化位点还会按其相对 TM7 的位置标注**近端**（TM7 后 ≤15aa，主要驱动脱敏）或**远端**（>15aa，主要驱动内吞与 arrestin 介导的信号）——对应 arrestin 的近端/远端 micro-locks 机制 [12]。

3. **保守性与评分** —— 功能性 PTM 位点通常跨物种保守，因此与其他物种的直系同源比对，保守的位点权重更高，最后统一排序。得分为0–1的排序值，**只用于排序，不代表概率**；置信度High/Medium/Low由motif强度+保守性决定。共识motif只是必要而非充分条件，请把排序当作假设参考。

每条预测位点都附带支撑该规则的文献依据（HTML 报告中以可点击 DOI 形式给出，如 pXpp 位点对应 Isaikina 2023 [7]、C 端尾位点对应 Sente 2018 [12] 等），便于直接核验规则来源。

### 文献核验

除了上述预测规则的来源文献，程序对**数据库已有点位**的文献引用也会做核验：抓取所引用 PMID 的 PubMed 摘要，确认是否真正支持"**该残基+该PTM类型**"，给出三种判定：

- **直接** —— 摘要点名了该残基与该修饰。
- **间接** —— 论文确实在讲该蛋白的该修饰，但未点名该残基。
- **不支持** —— 摘要指向的是别的残基或别的蛋白（例如把Tyr364的论文挂在Ser364位点上）。这类情况往往暴露了数据库的误归属。

### 局限性

这套规则本质是序列 motif + 膜拓扑 + 跨物种保守性的启发式组合，不是机器学习模型（如 NetPhos 或 PhosphoSitePlus 的 SVM/DNN），也没有做溶剂可及性计算（用"胞内环 + C 端尾"作膜可及性代理）。具体到各类 PTM：

- **磷酸化**没有强线性 motif，假阳性偏高——这也是程序设置保守度权重而非单纯 motif 匹配的原因。pXpp 三磷酸簇的判组过滤依赖 Table S2 的 26 个受体（一手实验数据），未命中的受体走 ICL3 长度回退：只有 ICL3 明确很短（<~5aa）或很长（>100aa）才能判组，中间地带（5–100aa）无法判组、两区 pXpp 均保留（不过滤）。"多磷酸化簇"标注是基于局部 S/T 密度的代理，真正的 barcode 取决于磷酸化位点的组合与排列而非单纯密度 [13]，所以这些标注只能定位候选区，不能判定确切的磷酸化模式。CK1 规则同样只能作弱提示：其共识为 pS/pT-X-X-S/T，要求候选位点**上游第 3 位**存在已磷酸化的引物残基 [14]，而纯序列无法确认引物是否真的被磷酸化，因此 CK1 命中在所有激酶 motif 中权重最低、恒为低置信。
- **泛素化**无可靠线性共识，PPxY 又常位于 adaptor（ARRDC3 等）而非 GPCR 自身，故 adaptor 介导的泛素化无法检出；该规则仅作保守弱提示，恒低置信、需实验确认。
- **保守性**通过精确基因名查直系同源，主要命中哺乳动物，保守区域易饱和到接近 1.0。

总体而言，所有预测位点都应视为可检验的假设，而非确证结果；排序用于优先级参考。

---

## 数据来源

- **UniProt** —— 序列、拓扑、PTM注释与证据代码
- **iPTMnet** —— 带PMID的文献PTM位点
- **dbPTM** —— 可选的离线PTM数据集
- **PubMed** —— 逐条核验所用的摘要

## License

本项目采用 [MIT License](LICENSE) 开源。引用结果时请注明底层数据来源（UniProt、iPTMnet、dbPTM、PubMed）。

---

## 参考文献

正文中 [n] 标注对应的文献清单如下（引文信息均经 PubMed / Europe PMC / AMiner 核验）。

1. Kennelly PJ, Krebs EG. Consensus sequences as substrate specificity determinants for protein kinases and protein phosphatases. *J Biol Chem*. 1991;266(24):15555–15558. doi: [10.1016/S0021-9258(18)98436-X](https://doi.org/10.1016/S0021-9258(18)98436-X)
2. Pinna LA, Ruzzene M. How do protein kinases recognize their substrates? *Biochim Biophys Acta*. 1996;1314(3):191–225. doi: [10.1016/S0167-4889(96)00083-3](https://doi.org/10.1016/S0167-4889(96)00083-3)
3. Tobin AB. G-protein-coupled receptor phosphorylation: where, when and by whom. *Br J Pharmacol*. 2008;153(Suppl 1):S167–S176. doi: [10.1038/sj.bjp.0707662](https://doi.org/10.1038/sj.bjp.0707662)
4. Komolov KE, Benovic JL. G protein-coupled receptor kinases: past, present and future. *Cell Signal*. 2018;41:17–24. doi: [10.1016/j.cellsig.2017.07.004](https://doi.org/10.1016/j.cellsig.2017.07.004)
5. Onorato JJ, Palczewski K, Regan JW, Caron MG, Lefkowitz RJ, Benovic JL. Role of acidic amino acids in peptide substrates of the beta-adrenergic receptor kinase and rhodopsin kinase. *Biochemistry*. 1991;30(21):5118–5125. doi: [10.1021/bi00235a002](https://doi.org/10.1021/bi00235a002)
6. Asai D, Toita R, Murata M, Katayama Y, Nakashima H, Kang JH. Peptide substrates for G protein-coupled receptor kinase 2. *FEBS Lett*. 2014;588(13):2129–2132. doi: [10.1016/j.febslet.2014.04.038](https://doi.org/10.1016/j.febslet.2014.04.038)
7. Isaikina P, Petrovic I, Jakob RP, et al. A key GPCR phosphorylation motif discovered in arrestin2⋅CCR5 phosphopeptide complexes. *Mol Cell*. 2023;83(12):2108–2121.e7. doi: [10.1016/j.molcel.2023.05.002](https://doi.org/10.1016/j.molcel.2023.05.002)
8. García Rodríguez C, Cundell DR, Tuomanen EI, Kolakowski LF Jr, Gerard C, Gerard NP. The role of N-glycosylation for functional expression of the human platelet-activating factor receptor. *J Biol Chem*. 1995;270(42):25178–25184. doi: [10.1074/jbc.270.42.25178](https://doi.org/10.1074/jbc.270.42.25178)
9. Hussain W, Khan YD, Rasool N, Khan SA, Chou KC. SPalmitoylC-PseAAC: a sequence-based model developed via Chou's 5-steps rule and general PseAAC for identifying S-palmitoylation sites in proteins. *Anal Biochem*. 2019;568:14–23. doi: [10.1016/j.ab.2018.12.019](https://doi.org/10.1016/j.ab.2018.12.019)
10. Kennedy JE, Marchese A. Regulation of GPCR trafficking by ubiquitin. *Prog Mol Biol Transl Sci*. 2015;132:15–38. doi: [10.1016/bs.pmbts.2015.02.005](https://doi.org/10.1016/bs.pmbts.2015.02.005)
11. Min B. In silico identification and in vitro validation of NEDD4-mediated GPCR ubiquitination. 2012. 见 [AMiner 条目](https://www.aminer.cn/pub/56d87a47dabfae2eee362f64)
12. Sente A, Peer R, Srivastava A, et al. Molecular mechanism of modulating arrestin conformation by GPCR phosphorylation. *Nat Struct Mol Biol*. 2018;25(6):538–545. doi: [10.1038/s41594-018-0071-3](https://doi.org/10.1038/s41594-018-0071-3)
13. Latorraca NR, Masureel M, Hollingsworth SA, et al. How GPCR phosphorylation patterns orchestrate arrestin-mediated signaling. *Cell*. 2020;183(7):1813–1825.e18. doi: [10.1016/j.cell.2020.11.014](https://doi.org/10.1016/j.cell.2020.11.014)
14. Venerando A, Ruzzene M, Pinna LA. Casein kinase: the triple meaning of a misnomer. *Biochem J*. 2014;460(2):141–156. doi: [10.1042/BJ20140178](https://doi.org/10.1042/BJ20140178)
