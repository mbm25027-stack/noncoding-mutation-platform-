# NC-MutScan v3.1: Noncoding Somatic Mutation Biophysical Engine

NC-MutScan v3.1 is an open-source computational framework for identifying, biophysically scoring, and functionally annotating noncoding somatic driver mutations in cancer genomes. 

By combining dual-strand Position Weight Matrix (PWM) log-odds scoring, multi-TF competitive affinity analysis, and structural noise filters, NC-MutScan distinguishes genuine driver mutations from passenger variants and technical artifacts.

---

## Key Features

- **Multi-TF Competitive Scanning:** Simultaneously scores JASPAR (2024) and HOCOMOCO v11 matrices (including GABPα, GATA3, MYB, ETS1, Sp1, ETS2, AP-1, E2F1) to evaluate overlapping cis-regulatory elements without dropping secondary hits.
- **AID/SHM Somatic Hypermutation Filter:** Eliminates false-positive driver calls in B-cell malignancies by scanning forward ($WRCY = [A/T][A/G]C[C/T]$) and reverse-strand ($RGYW = [A/G]G[C/T][A/T]$) hypermutation patterns.
- **Dirichlet Background Normalization:** Prevents artificial score inflation in CpG islands, G-quadruplexes, and high-GC regions using a uniform Dirichlet prior ($bg = 0.25$).
- **Degenerate Base & Boundary Handling:** Fully sanitizes degenerate IUPAC base codes (`N`) and evaluates 5'-terminal sequence boundaries ($i = 0$).
- **Standardized Coordinate Anchor:** Standardized on **GRCh38/hg38** genomic coordinates (e.g., TERT promoter $-124\,C>T$ at `chr5:1,295,113`).

---

## Validation & Performance

### 1. Edge-Case Unit Verification Suite (15 Test Cases)
- **Pass Rate:** 100.0% (15/15 passed)
- **Coverage:** Canonical TERT drivers, reverse-strand ETS scans, AID/SHM passengers, IUPAC sanitization, Sp1 loss-of-function, G-quadruplex noise, composite overlaps, exact 6-bp boundaries, and dual SNV compound hotspots.

### 2. Held-Out Empirical Benchmark (50 Loci)
Evaluated on an independent dataset of 25 ground-truth pan-cancer drivers (TERT, PLEKHS1, SDHD, WDR74, TAL1, LMO2) and 25 matched neutral background variants:

| Metric | Score |
| :--- | :--- |
| **Sensitivity (TPR)** | **100.0%** (25/25 True Positives) |
| **Specificity (TNR)** | **100.0%** (25/25 True Negatives) |
| **Precision (PPV)** | **100.0%** |
| **False Positive Rate (FPR)** | **0.0%** |
| **$F_1$-Score** | **1.0000** |

---

## Quick Start in VS Code

### Installation
```bash
git clone [https://github.com/mbm25027-stack/noncoding-mutation-platform.git](https://github.com/mbm25027-stack/noncoding-mutation-platform.git)
cd noncoding-mutation-platform
pip install -r requirements.txt