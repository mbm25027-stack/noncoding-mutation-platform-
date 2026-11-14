import datetime
import os
import shutil

ROOT = os.path.abspath(os.getcwd())
BUNDLE_DIR = os.path.join(ROOT, "B1_B5_SELF_RESEARCH_BUNDLE")
FILES_DIR = os.path.join(BUNDLE_DIR, "copied_source_files")
PACK_TXT = os.path.join(BUNDLE_DIR, "B1_B5_MASTER_EVIDENCE_PACK.txt")
SKIP_DIRS = {"venv", ".venv", ".git", "__pycache__", "B1_B5_SELF_RESEARCH_BUNDLE"}

KEYWORDS = (
    "spliceai", "cadd", "exon_dist", "blind", "loco", "mpra", "kircher",
    "clinvar", "nc_mutscan", "table2", "sha256", "flank", "dev_chr", "locked"
)


def count_vcf_by_chr(path):
    counts = {}
    total = 0
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t", 2)
                c = parts[0].strip()
                c_norm = c[3:] if c.lower().startswith("chr") else c
                counts[c_norm] = counts.get(c_norm, 0) + 1
                total += 1
    except OSError:
        return 0, {}
    return total, counts


def collect_files():
    os.makedirs(FILES_DIR, exist_ok=True)
    matched = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".")]
        for fn in filenames:
            if fn.startswith("~$"):
                continue
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, ROOT)
            low = rel.lower()
            ext = os.path.splitext(fn)[1].lower()
            if ext in (".csv", ".tsv", ".json", ".vcf", ".txt", ".py", ".docx"):
                if any(k in low for k in KEYWORDS) or ext == ".docx":
                    matched.append((rel, full, ext))
    matched.sort(key=lambda x: x[0])

    copied_map = []
    for rel, full, ext in matched:
        safe_name = rel.replace("\\", "__").replace("/", "__")
        dst = os.path.join(FILES_DIR, safe_name)
        try:
            shutil.copy2(full, dst)
            copied_map.append((rel, safe_name, full, ext))
        except OSError:
            pass
    return copied_map


def dump_docx_sections(docx_path, out_lines):
    try:
        import docx
    except ImportError:
        out_lines.append("[INFO] python-docx not installed; skipping DOCX text dump.")
        return
    if not os.path.isfile(docx_path):
        return
    doc = docx.Document(docx_path)
    out_lines.append(f"=== CURRENT MANUSCRIPT TEXT & TABLES: {os.path.basename(docx_path)} ===")
    for p in doc.paragraphs:
        t = p.text.strip()
        if t:
            out_lines.append(t)
    out_lines.append("")
    for idx, tbl in enumerate(doc.tables, 1):
        out_lines.append(f"--- DOCX Table #{idx} ---")
        for row in tbl.rows:
            out_lines.append(" | ".join(c.text.strip().replace("\n", " ") for c in row.cells))
        out_lines.append("")


def main():
    copied = collect_files()
    out = []
    w = out.append

    w("=" * 96)
    w("B.1 - B.5 SELF-RESEARCH EVIDENCE PACK & FILE MANIFEST")
    w(f"Generated : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    w(f"SourceDir : {ROOT}")
    w(f"BundleDir : {BUNDLE_DIR}")
    w("=" * 96)
    w("")

    w("1. COPIED SOURCE FILES INDEX (All files placed in copied_source_files/)")
    w("-" * 96)
    for rel, safe_name, full, ext in copied:
        sz = os.path.getsize(full)
        w(f"  * {rel:60s} -> {safe_name} ({sz:,} bytes)")
    w("")

    w("2. ITEM B.5 - VCF & CLINVAR CHROMOSOME VARIANT COUNTS (chr1-11,X vs chr12-22)")
    w("-" * 96)
    chr_1_11_x = {str(i) for i in range(1, 12)} | {"X"}
    chr_12_22 = {str(i) for i in range(12, 23)}
    for rel, safe_name, full, ext in copied:
        if ext == ".vcf":
            tot, by_c = count_vcf_by_chr(full)
            n_1_11_x = sum(v for k, v in by_c.items() if k in chr_1_11_x)
            n_12_22 = sum(v for k, v in by_c.items() if k in chr_12_22)
            w(f"  VCF: {rel}")
            w(f"       Total={tot:,} | chr1-11,X={n_1_11_x:,} | chr12-22={n_12_22:,} | Per-chr={dict(sorted(by_c.items()))}")
    w("")

    w("3. ITEM B.4 - FROZEN SCANNER SPLICE PARAMETERS IN PYTHON CODE")
    w("-" * 96)
    for rel, safe_name, full, ext in copied:
        if ext == ".py" and "nc_mutscan" in rel.lower():
            w(f"  Script: {rel}")
            try:
                with open(full, "r", encoding="utf-8", errors="ignore") as f:
                    for lineno, line in enumerate(f, 1):
                        if any(k in line.lower() for k in ("acceptor", "donor", "cryptic", "6.4", "7.5", "2.8", "3.8", "5.6", "s_raw", "pwm")):
                            w(f"    L{lineno:04d}: {line.rstrip()}")
            except OSError:
                pass
            w("")

    w("4. ITEMS B.1, B.2, B.3 - FULL CONTENTS OF CSV / TSV / TXT BENCHMARK TABLES")
    w("-" * 96)
    for rel, safe_name, full, ext in copied:
        if ext in (".csv", ".tsv", ".txt") and os.path.getsize(full) < 512 * 1024:
            w("")
            w(f">>> FILE: {rel} <<<")
            try:
                with open(full, "r", encoding="utf-8", errors="ignore") as f:
                    for i, line in enumerate(f):
                        if i >= 250:
                            w("    ... [truncated at 250 lines; open copied CSV directly for full file]")
                            break
                        w("    " + line.rstrip())
            except OSError:
                pass
    w("")

    w("5. CURRENT MANUSCRIPT SECTIONS & TABLES (FOR SIDE-BY-SIDE COMPARISON)")
    w("-" * 96)
    target_docx = os.path.join(ROOT, "NC-MutScan_v3.3_PreFilter_Revised_Manuscript.docx")
    if not os.path.isfile(target_docx):
        for rel, _, full, ext in copied:
            if ext == ".docx":
                target_docx = full
                break
    dump_docx_sections(target_docx, out)

    with open(PACK_TXT, "w", encoding="utf-8") as out_f:
        out_f.write("\n".join(out) + "\n")

    print(f"[DONE] Copied {len(copied)} files into: {FILES_DIR}")
    print(f"[DONE] Master evidence pack written to: {PACK_TXT}")

if __name__ == "__main__":
    main()
