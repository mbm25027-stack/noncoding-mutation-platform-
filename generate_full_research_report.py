import csv
import datetime
import os
import re

ROOT = os.path.abspath(os.getcwd())
OUT_TXT = os.path.join(ROOT, "FULL_RESEARCH_REPORT_DAY1_TO_NOW.txt")
SKIP_DIRS = {"venv", ".venv", ".git", "__pycache__", "node_modules", ".idea", ".vscode"}
SQ3 = chr(39) * 3
DQ3 = chr(34) * 3


def fmt_ts(ts):
    return datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S")


def fmt_date(ts):
    return datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d")


def count_lines_fast(path):
    count = 0
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                count += chunk.count(b"\n")
    except OSError:
        return 0
    return count


def count_vcf_variants(path):
    n = 0
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    n += 1
    except OSError:
        return 0
    return n


def extract_py_summary(path):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read(4000)
        m = re.search(DQ3 + r"(.*?)" + DQ3, text, re.DOTALL)
        if not m:
            m = re.search(SQ3 + r"(.*?)" + SQ3, text, re.DOTALL)
        if m:
            doc = " ".join(m.group(1).strip().split())
            return doc[:180] + ("..." if len(doc) > 180 else "")
        comments = []
        for line in text.splitlines()[:25]:
            s = line.strip()
            if s.startswith("#") and not s.startswith("#!"):
                comments.append(s.lstrip("#").strip())
        if comments:
            joined = " ".join(comments)
            return joined[:180] + ("..." if len(joined) > 180 else "")
    except OSError:
        pass
    return "Python analysis / pipeline script"


def collect_inventory():
    files = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".")]
        for fn in filenames:
            if fn.startswith("~$"):
                continue
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, ROOT)
            try:
                st = os.stat(full)
                earliest = min(st.st_ctime, st.st_mtime)
                files.append({
                    "full": full,
                    "rel": rel,
                    "name": fn,
                    "ext": os.path.splitext(fn)[1].lower(),
                    "size": st.st_size,
                    "ctime": st.st_ctime,
                    "mtime": st.st_mtime,
                    "earliest": earliest,
                })
            except OSError:
                pass
    files.sort(key=lambda x: (x["earliest"], x["rel"]))
    return files


def read_table_rows(path, max_rows=500):
    delim = "\t" if path.lower().endswith(".tsv") else ","
    rows = []
    headers = []
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f, delimiter=delim)
            headers = reader.fieldnames or []
            for i, r in enumerate(reader):
                if i >= max_rows:
                    break
                rows.append(r)
    except Exception:
        pass
    return headers, rows


def extract_loci_and_genes(files):
    gene_keys = ["Gene", "Hugo_Symbol", "Target_Gene", "Element", "Locus", "gene", "target_gene"]
    chr_keys = ["Chromosome", "Chrom", "chrom", "chr"]
    pos_keys = [
        "Start_Position", "Position", "Peak_Center", "Center_Position",
        "GRCh38_Coord", "hg19_Coord", "Locus_GRCh38", "Locus_hg19", "pos"
    ]
    end_keys = ["End_Position", "Window_End", "Cluster_End"]
    extra_keys = [
        "Cancer_Type", "Cancer_Types", "Projects", "Dataset", "Element_Class", "Biotype",
        "Region_Group", "Variant_Classification", "Mutated_Donors", "Recurrent_Donors",
        "Observed_Mutations", "Observed_SNVs", "Fold_Enrichment", "Fold_Enrichment_Repaired",
        "Obs_Over_Exp", "q_value", "Q_NB_ZeroInc", "Q_Global_Cluster_NB", "Q_Local_Cluster_NB",
        "Q_eQTL_BH", "CovAdj_Beta_Mut", "Log2_Fold_Change", "Precision_PPV_%", "Sensitivity_%",
        "Specificity_%", "Classification", "Mechanism", "Predicted_TF_Impact"
    ]

    loci_records = []
    seen = set()

    for f in files:
        if f["ext"] not in (".csv", ".tsv") or f["size"] > 15 * 1024 * 1024:
            continue
        headers, rows = read_table_rows(f["full"], max_rows=400)
        if not headers or not rows:
            continue

        g_col = next((k for k in gene_keys if k in headers), None)
        c_col = next((k for k in chr_keys if k in headers), None)
        p_col = next((k for k in pos_keys if k in headers), None)
        e_col = next((k for k in end_keys if k in headers), None)

        if not g_col and not (c_col and p_col):
            continue

        for r in rows:
            gene = str(r.get(g_col, "")).strip() if g_col else ""
            chrom = str(r.get(c_col, "")).strip() if c_col else ""
            pos = str(r.get(p_col, "")).strip() if p_col else ""
            endp = str(r.get(e_col, "")).strip() if e_col else ""

            if not gene and not pos:
                continue
            if gene.startswith("Pooled") or gene.lower() in ("", "na", "none", "nan"):
                if not pos:
                    continue

            coord = ""
            if chrom and pos and ":" not in pos:
                coord = f"{chrom}:{pos}-{endp}" if endp and endp != pos else f"{chrom}:{pos}"
            elif pos:
                coord = pos
            elif chrom:
                coord = chrom

            details = []
            for ek in extra_keys:
                if ek in headers:
                    val = str(r.get(ek, "")).strip()
                    if val and val.lower() not in ("na", "nan", "none", "null", ""):
                        details.append(f"{ek}={val}")

            sig_key = (f["rel"], gene, coord)
            if sig_key in seen:
                continue
            seen.add(sig_key)

            loci_records.append({
                "source": f["rel"],
                "gene": gene or "(intergenic/window)",
                "coord": coord or "(locus-level)",
                "details": "; ".join(details[:8]),
            })

    return loci_records


def extract_key_metric_tables(files):
    priority_names = {
        "41_manuscript_final_numbers.csv",
        "01_dataset_summary.csv",
        "04_region_groups.csv",
        "32_melanoma_nb_bmr_calibration.csv",
        "34_breast_nb_bmr_calibration.csv",
        "39_liver_nb_bmr_calibration.csv",
        "36_ml_loocv_benchmark_metrics.csv",
        "43_loocv_bootstrap_and_permutation.csv",
        "Table2_Blind_Ablations_SpliceAI.csv",
        "Table4_PCAWG_Melanoma_WGS_Hotspots.csv",
        "pcawg_histology_adjusted_eqtl.csv",
        "Table16_blind_FN_mechanisms.csv",
        "kircher_mpra_multi_element_benchmark.csv",
    }
    extracted = []
    for f in files:
        if f["name"] in priority_names:
            headers, rows = read_table_rows(f["full"], max_rows=100)
            if headers and rows:
                extracted.append((f["rel"], headers, rows))
    return extracted


def build_report():
    files = collect_inventory()
    if not files:
        print("[ERROR] No files found in current directory.")
        return

    first_ts = files[0]["earliest"]
    last_ts = max(f["mtime"] for f in files)

    py_files = [f for f in files if f["ext"] == ".py"]
    data_files = [
        f for f in files
        if f["ext"] in (".csv", ".tsv", ".vcf", ".json", ".bed", ".maf", ".txt")
        and f["name"] != os.path.basename(OUT_TXT)
    ]
    docx_files = [f for f in files if f["ext"] == ".docx"]
    fig_files = [f for f in files if f["ext"] in (".png", ".pdf", ".svg")]

    loci_records = extract_loci_and_genes(files)
    metric_tables = extract_key_metric_tables(files)

    lines = []
    w = lines.append

    w("=" * 98)
    w("COMPREHENSIVE RESEARCH REPORT: FROM DAY 1 TO PRESENT")
    w(f"Workspace Directory : {ROOT}")
    w(f"Report Generated    : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    w(f"Earliest File Stamp : {fmt_ts(first_ts)}")
    w(f"Latest File Stamp   : {fmt_ts(last_ts)}")
    w(
        f"Total Project Files : {len(files)} ({len(py_files)} Python scripts, "
        f"{len(data_files)} data/table files, {len(fig_files)} figures, {len(docx_files)} manuscripts)"
    )
    w("=" * 98)
    w("")

    w("1. EXECUTIVE OVERVIEW OF THE RESEARCH JOURNEY (PHASE BY PHASE)")
    w("-" * 98)
    w("This workspace documents a two-stage computational genomics research program investigating")
    w("noncoding somatic mutations and regulatory variant pathogenicity:")
    w("")
    w("PHASE 1 (Early Exploration -> Statistical Calibration -> ML Driver Classification):")
    w("  - Objective: Discover recurrent noncoding driver hotspots across pan-cancer cohorts")
    w("    (Melanoma/TCGA-SKCM/MELA-AU/PCAWG, Breast, Liver, Blood, Brain, Colon, Lung).")
    w("  - Core Methodology:")
    w("      1. Stratified noncoding mutations across regulatory regions (Promoters, 5'/3'-UTRs,")
    w("         Splice sites, Introns, lncRNAs, miRNAs, Enhancers, Intergenic regions).")
    w("      2. Identified severe overdispersion in raw somatic mutation counts (Fano factor >> 1),")
    w("         showing that uncalibrated Poisson models produce massive false-positive inflation.")
    w("      3. Built Negative Binomial (NB) and Zero-Inflated Negative Binomial (ZINB) background")
    w("         mutation rate (BMR) models accounting for local sequence context and mutational")
    w("         signatures (UV C>T / G>A pyrimidine dimer enrichment in melanoma; AID/SHM WRCY/RGYW")
    w("         hypermutation in lymphoid/blood cancers).")
    w("      4. Validated candidate hotspots using cis-eQTL transcriptional association (unadjusted")
    w("         and histology/covariate-adjusted linear models with Benjamini-Hochberg FDR).")
    w("      5. Trained supervised Machine Learning classifiers (Gradient Boosted Decision Trees")
    w("         and L2-regularized Logistic Regression) evaluated via Leave-One-Out Cross-Validation")
    w("         (LOOCV) with bootstrap 95% CIs and permutation testing.")
    w("")
    w("PHASE 2 (Deterministic Rule-In Scanner Engineering -> Locked Validation -> Multi-Element MPRA Audit):")
    w("  - Objective: Build NC-MutScan v3.3 (formerly prototyped as NC-DriverEngine v3.1-v3.3),")
    w("    a zero-dependency, CPU-only deterministic sequence scanner for high-precision rule-in")
    w("    prioritization of noncoding SNVs without relying on precomputed genome-wide databases.")
    w("  - Core Architecture (Frozen NC-MutScan v3.3):")
    w("      1. Module 1 - Human U1/U2 5' Splice Donor (9-mer, -3..+6) & 3' Splice Acceptor")
    w("         (23-mer, -20..+3) PWMs (Burge et al. 1999 / Yeo & Burge 2004 consensus):")
    w("         Splice loss: |delta_PWM| >= 2.8 bits (WT >= 5.6 bits); Cryptic gain: delta_PWM >= +3.8 bits")
    w("         (Mut >= 6.4 bits; canonical anchor >= 7.5 bits).")
    w("      2. Module 2 - 3'-UTR Polyadenylation Signal (PAS) Loss Scanner:")
    w("         Canonical AATAAA / ATTAAA disruption in >= 65% A+T 21-bp window (4.2-bit floor).")
    w("      3. Module 3 - 9-TF Competitive PWM Arbitration + AID/SHM Passenger Masking:")
    w("         Sp1, GABPA, MYB, ETS1, E2F1, AP-1, GATA3 (MA0037.4), FOXA1 (MA0148.4), CTCF (MA0139.1)")
    w("         with |delta_PWM_TF| >= 4.2 bits (active site >= 7.8 bits) and WRCY/RGYW passenger filter.")
    w("  - Benchmarks & Audits Completed:")
    w("      * Cryptographic SHA-256 code freeze on dev_chr1_11.vcf (N=150 drawn from N=955 pool)")
    w("        before unlocking locked_blind_chr12_22.vcf (N=150: 75 Pathogenic, 75 Benign).")
    w("      * Locked Blind ClinVar Evaluation (N=150): 98.67% Specificity (95% CI: 92.8-99.8%),")
    w("        97.37% PPV, 49.33% Sensitivity (37/75 TP, 1/75 FP), 495.7 us/variant latency.")
    w("      * 23-Chromosome Leave-One-Chromosome-Out (LOCO, chr1-22,X; N=300): 94.00% Specificity")
    w("        (95% CI: 89.0-96.8%), 44.67% Sensitivity, pooled LR+ = 7.44 (95% CI: 3.85-14.39).")
    w("      * Matched ClinVar Comparator Benchmark (N=139): Pre-threshold continuous ROC-AUC = 0.785")
    w("        (thresholded S_unified AUC = 0.733) vs. SpliceAI (0.946), CADD (0.871), MaxEntScan (0.832),")
    w("        and exon-junction proximity (0.819).")
    w("      * Quantitative False-Negative Audit (all 38 FNs on chr12-22): 57.9% (22/38) in extended")
    w("        +/-3-17 bp splice regions (15 sub-threshold delta_PWM, 7 outside PWM window); 28.9% (11/38)")
    w("        at canonical +/-1-2 bp sites below the 5.6-bit WT floor (including GC-AG junctions).")
    w("      * Saturation-Mutagenesis MPRA Audit across all 29 Kircher et al. (2019) GRCh38 libraries")
    w("        (38,556 library SNV observations; 31,482 unique genomic SNVs across 20 elements):")
    w("        95.25% pooled specificity against neutral controls (|log2| < 0.10); precision enriched")
    w("        in in-registry TERT promoter controls (53.3-66.7% PPV) while suppressing false positives")
    w("        across unmodeled enhancers (pooled LR+ = 0.78).")
    w("")

    w("2. CHRONOLOGICAL RESEARCH TIMELINE (DAY 1 TO PRESENT)")
    w("-" * 98)
    by_date = {}
    for f in files:
        d = fmt_date(f["earliest"])
        by_date.setdefault(d, []).append(f)

    for d in sorted(by_date.keys()):
        day_files = by_date[d]
        w(f"[{d}]  ({len(day_files)} files created/updated)")
        for f in day_files[:25]:
            w(f"    - {f['rel']}  ({f['size']:,} bytes)")
        if len(day_files) > 25:
            w(f"    ... plus {len(day_files) - 25} additional files on this date")
    w("")

    w("3. COMPLETE INVENTORY OF DATASETS, COHORTS, AND INPUT/OUTPUT DATA FILES")
    w("-" * 98)
    w(f"{'Relative Path':64s} | {'Format':6s} | {'Rows/Vars':10s} | {'Columns / Description'}")
    w("-" * 98)
    for f in data_files:
        ext = f["ext"]
        if ext == ".vcf":
            n_rows = count_vcf_variants(f["full"])
            desc = "VCF variant dataset"
        elif ext in (".csv", ".tsv"):
            n_rows = max(0, count_lines_fast(f["full"]) - 1)
            hdrs, _ = read_table_rows(f["full"], max_rows=1)
            desc = ", ".join(hdrs[:6]) + (" ..." if len(hdrs) > 6 else "")
        elif ext == ".json":
            n_rows = count_lines_fast(f["full"])
            desc = "JSON cache / annotations / scores"
        else:
            n_rows = count_lines_fast(f["full"])
            desc = "Text / manifest / auxiliary data"
        w(f"{f['rel'][:64]:64s} | {ext:6s} | {n_rows:10,d} | {desc[:60]}")
    w("")

    w("4. COMPLETE COMPENDIUM OF GENOMIC LOCI, GENES, HOTSPOTS, eQTLs, AND MPRA ELEMENTS")
    w("-" * 98)
    w(f"Total distinct (source, gene/element, locus) entries extracted: {len(loci_records)}")
    w("")
    unique_genes = sorted(set(r["gene"] for r in loci_records if r["gene"] and not r["gene"].startswith("(")))
    w(f"Unique Genes / Regulatory Elements Studied ({len(unique_genes)}):")
    for i in range(0, len(unique_genes), 8):
        w("  " + ", ".join(unique_genes[i:i + 8]))
    w("")
    w("Detailed Locus-by-Locus Breakdown Across Result Tables:")
    w("-" * 98)
    current_src = None
    for rec in loci_records:
        if rec["source"] != current_src:
            current_src = rec["source"]
            w(f"\n[Source Table: {current_src}]")
        w(f"  * {rec['gene']:24s} | {rec['coord']:28s} | {rec['details']}")
    w("")

    w("5. VERIFIED QUANTITATIVE BENCHMARKS & STATISTICAL METRICS")
    w("-" * 98)
    for src_rel, headers, rows in metric_tables:
        w(f"\n>>> Table: {src_rel}  (Columns: {', '.join(headers)})")
        for r in rows[:45]:
            row_str = " | ".join(f"{h}={r.get(h, '')}" for h in headers if str(r.get(h, "")).strip() != "")
            w(f"    {row_str}")
        if len(rows) > 45:
            w(f"    ... ({len(rows) - 45} more rows in {src_rel})")
    w("")

    w("6. COMPLETE CATALOG OF PYTHON SCRIPTS, MANUSCRIPTS, AND FIGURES")
    w("-" * 98)
    w("A. Python Scripts (Chronological Order):")
    for f in py_files:
        desc = extract_py_summary(f["full"])
        w(f"  - [{fmt_date(f['earliest'])}] {f['rel']:48s} | {desc}")
    w("")
    w("B. Manuscript Documents (.docx):")
    for f in docx_files:
        w(f"  - [{fmt_date(f['mtime'])}] {f['rel']:58s} ({f['size']/1024:.1f} KB)")
    w("")
    w("C. Generated Figures (.png / .pdf):")
    for f in fig_files:
        w(f"  - [{fmt_date(f['mtime'])}] {f['rel']:58s} ({f['size']/1024:.1f} KB)")
    w("")
    w("=" * 98)
    w("END OF FULL RESEARCH REPORT")
    w("=" * 98)

    with open(OUT_TXT, "w", encoding="utf-8") as out_f:
        out_f.write("\n".join(lines) + "\n")

    print(f"[SUCCESS] Full research report written to: {OUT_TXT}")
    print(
        f"[STATS] Files scanned: {len(files)} | Datasets: {len(data_files)} | "
        f"Unique Genes/Elements: {len(unique_genes)} | Loci records: {len(loci_records)}"
    )


if __name__ == "__main__":
    build_report()
