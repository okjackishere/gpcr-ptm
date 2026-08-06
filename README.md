# GPCR-PTM

Predict and verify **post-translational modification (PTM) sites** on G-protein-coupled receptors (GPCRs). Given a GPCR (gene name, UniProt accession, or protein name), it outputs both **verified** PTM sites (curated from public databases) and **predicted** PTM sites (rule-based candidates), and it checks every cited publication to confirm it really supports each site.

Two ways to use it: a **web interface** (Flask, with progress bar) or a **command line**.

---

## Method

### What kind of prediction is this?

**Rule-based, not machine learning.** The program does not train on data or learn parameters. It applies explicit, well-established biological rules, so every result has a traceable reason and is fully reproducible.

### How verified sites are obtained (L1 / L2)

Verified sites come from curated databases — UniProt, iPTMnet and dbPTM — and are graded by **evidence strength**:

- **Experimental evidence** (reported in the literature and curated in UniProt) is ranked highest.
- **Inferred / by-similarity annotations** and **single high-throughput hits** are ranked lower and clearly labelled (the latter are flagged "needs experimental confirmation").

### How prediction works (L3)

For modifications not yet recorded in any database, candidates are produced in three steps:

1. **Motif scan** — search the receptor sequence for well-established consensus motifs (phosphorylation, N-glycosylation sequons, palmitoylation).
2. **Topology filter** — keep only sites in the biochemically correct location: N-glycosylation must be extracellular; phosphorylation / palmitoylation must be intracellular. This removes most false positives.
3. **Conservation & scoring** — compare the receptor with orthologs from other species; sites conserved across species get more weight, then all candidates are scored and ranked.

### Literature verification (every cited paper is checked)

For every PMID cited by a database, the program fetches the paper's PubMed abstract and judges whether it actually supports **that residue + that PTM type**, returning one of:

- **Direct** — the abstract names this residue and this modification.
- **Indirect** — the paper is about this protein and this modification, but does not name the residue.
- **Unsupported** — the abstract points to a different residue or a different protein (e.g. a Tyr364 paper attached to a Ser364 site). This usually flags a database mis-annotation.

### Why rule-based instead of machine learning?

- **Transparent** — each site carries the motif, region and evidence behind it.
- **Reproducible** — same input, same output, no randomness.
- The PTM consensus motifs of GPCRs are already well characterised, so explicit rules are appropriate and interpretable.

---

## Usage

### Quick start

**Linux / macOS**

```bash
bash run.sh          # installs everything on first run, then starts the web server
# open http://127.0.0.1:8000
```

**Windows** — double-click `run.bat`, then open http://127.0.0.1:8000

**Manual install**

```bash
python3 -m venv --without-pip venv
curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
venv/bin/python get-pip.py
venv/bin/pip install -r requirements.txt
venv/bin/python webapp.py
```

Requires Python 3.8+.

### Web interface

Open http://127.0.0.1:8000, type a GPCR (e.g. `ADRB2`, `P07550`, `OPRM1`), and wait for the progress bar. Results are cached, so repeat queries are instant.

### Command line

```bash
python main.py ADRB2
python main.py P07550
python main.py OPRM1 --verbose      # also print per-PMID links/evidence in the terminal
```

---

## Reading the Output

### The three layers

Results are organised into three layers; read them in this order:

- **L1 · Verified** — experimentally reported sites, curated with PMIDs. Highest confidence; treat these as known modifications.
- **L2 · Supported** — support exists but is not fully settled: inferred / by-similarity annotations, database or literature hits, or single high-throughput detections (flagged "needs experimental confirmation"). Likely but unconfirmed.
- **L3 · Predicted** — rule-based candidates not found in any database. Ranked suggestions for follow-up, not established facts.

### The literature verdicts

Each PMID is marked **Direct** (green) / **Indirect** (amber) / **Unsupported** (red). Use them to judge how solid a site's citation list really is — "Unsupported" usually points to a database mis-annotation.

### Score & confidence

- **Score** — a ranking value (0–1); it only orders candidates, it is **not** a probability.
- **Confidence** — High / Medium / Low, derived from motif strength + conservation.

### Output files

| File | What it is |
|------|-----------|
| `P07550_report.html` | Interactive report: collapsible site cards, clickable PubMed links, sequence overview — open in any browser |
| `P07550_ptm.json` | Structured data for further analysis |
| `P07550_verification.md` | Literature-verification table with links and abstract evidence |

---

## Scientific Basis & Caveats

- Consensus motifs are **necessary but not sufficient** — e.g. only ~30–40% of N-glycosylation sequons are actually used, and phospho motifs are weak predictors. Scores are for ranking only.
- Conservation is computed against orthologs found by exact gene name (mostly mammals), so conserved regions saturate near 1.0.
- An **"Unsupported"** verdict means the cited paper does not support that residue + PTM; it does not necessarily mean the site is wrong.
- Predictions are hypotheses and should be confirmed experimentally (e.g. by mass spectrometry) before drawing conclusions.

## Data Sources

- **UniProt** — sequence, topology, PTM annotations and evidence codes
- **iPTMnet** — literature-curated PTM sites (with PMIDs)
- **dbPTM** — optional offline PTM dataset
- **PubMed** — abstracts used for per-PMID verification

## License

Released under the [MIT License](LICENSE). Please cite the underlying data sources (UniProt, iPTMnet, dbPTM, PubMed) when using results.
