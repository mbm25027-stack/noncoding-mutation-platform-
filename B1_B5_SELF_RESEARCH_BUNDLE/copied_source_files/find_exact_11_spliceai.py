"""
find_exact_11_spliceai.py
Locates the per-variant SpliceAI table in the workspace, identifies the exact
11 locked chr12-22 ClinVar variants (6 pathogenic, 5 benign) without SpliceAI
scores, prints them, and replaces '[variant identifiers to be supplied]'
(removing yellow highlight) in the active manuscript .docx.
"""

import csv
import glob
import json
import os
import re
import docx

# 1. Inspect and parse locked_blind_chr12_22.vcf
with open("locked_blind_chr12_22.vcf", "r", encoding="utf-8", errors="replace") as fh:
    raw_lines = [ln.rstrip("\r\n") for ln in fh if ln.strip()]

print("Sample lines from locked_blind_chr12_22.vcf:")
for ln in raw_lines[:5]:
    print("  ", ln[:140])

blind_variants = []
for ln in raw_lines:
    if ln.startswith("#"):
        continue
    parts = ln.split("\t")
    # Handle either standard VCF (CHROM POS ID REF ALT ...) or prefixed label VCF (LABEL CHROM POS ID REF ALT ...)
    if parts[0].upper() in ("PATHOGENIC", "BENIGN", "LIKELY_PATHOGENIC", "LIKELY_BENIGN"):
        label = "pathogenic" if "PATHOGENIC" in parts[0].upper() else "benign"
        chrom, pos, vid, ref, alt = parts[1], parts[2], parts[3], parts[4], parts[5]
    else:
        chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        info = parts[7] if len(parts) > 7 else ln
        label = "pathogenic" if "pathogenic" in info.lower() else "benign"

    c_clean = chrom.replace("chr", "")
    clinvar_id = vid if vid != "." else f"chr{c_clean}:{pos}:{ref}>{alt}"
    disp = f"ClinVar ID {clinvar_id} (chr{c_clean}:{pos} {ref}>{alt})" if vid != "." else f"chr{c_clean}:{pos} {ref}>{alt}"
    blind_variants.append({
        "vid": vid,
        "coord_key": f"{c_clean}:{pos}:{ref}:{alt}",
        "pos_key": f"{c_clean}:{pos}",
        "disp": disp,
        "label": label,
    })

print(f"\nParsed {len(blind_variants)} blind variants ({sum(v['label']=='pathogenic' for v in blind_variants)} pathogenic, {sum(v['label']=='benign' for v in blind_variants)} benign).")

# 2. Find all workspace files containing SpliceAI per-variant scores (search for ClinVar ID 31336918 or DS_max/SpliceAI)
print("\nScanning workspace files for per-variant SpliceAI scores...")
candidate_files = []
for root, _, files in os.walk("."):
    if "venv" in root or ".git" in root:
        continue
    for fn in files:
        if fn.endswith((".csv", ".tsv", ".json", ".txt", ".vcf", ".py")):
            fp = os.path.join(root, fn)
            if os.path.abspath(fp) == os.path.abspath("locked_blind_chr12_22.vcf"):
                continue
            try:
                with open(fp, "r", encoding="utf-8", errors="replace") as f:
                    txt = f.read()
                if "31336918" in txt or ("spliceai" in txt.lower() and ("ds_max" in txt.lower() or "ds_ag" in txt.lower() or "139" in txt)):
                    candidate_files.append((fp, txt))
                    print("  Matched file:", fp)
            except Exception:
                pass

# 3. Extract scored variant IDs / coordinates from matched data files
scored_vids = set()
scored_coords = set()

for fp, txt in candidate_files:
    if fp.endswith(".py"):
        continue
    lines = [l for l in txt.splitlines() if l.strip()]
    if len(lines) < 50:
        continue
    delim = "\t" if fp.endswith((".tsv", ".vcf")) else ","
    reader = csv.DictReader(lines, delimiter=delim)
    if not reader.fieldnames:
        continue
    fmap = {c.lower().strip(): c for c in reader.fieldnames}
    splice_cols = [orig for low, orig in fmap.items() if any(k in low for k in ("spliceai", "ds_max", "ds_ag", "ds_al", "ds_dg", "ds_dl"))]
    if not splice_cols:
        continue
    print(f"  Using SpliceAI columns {splice_cols} from {fp}")
    for row in reader:
        # Check if SpliceAI has a valid numeric score (not empty, NA, None, nan, '.')
        has_score = False
        for sc in splice_cols:
            val = str(row.get(sc, "")).strip()
            if val and val.lower() not in ("na", "nan", "none", "null", ".", "-"):
                try:
                    float(val.split("|")[0])
                    has_score = True
                    break
                except ValueError:
                    if re.search(r"\d+\.\d+", val):
                        has_score = True
                        break
        if has_score:
            for col_name, val in row.items():
                v_str = str(val).strip()
                if v_str.isdigit():
                    scored_vids.add(v_str)
            # Also capture chrom:pos if present
            c_col = next((orig for low, orig in fmap.items() if low in ("chrom", "chr", "#chrom")), None)
            p_col = next((orig for low, orig in fmap.items() if low in ("pos", "position")), None)
            if c_col and p_col:
                scored_coords.add(f"{str(row[c_col]).replace('chr','')}:{str(row[p_col]).strip()}")

unscored = [
    v for v in blind_variants
    if (v["vid"] not in scored_vids) and (v["pos_key"] not in scored_coords)
]

print(f"\nUnscored count: {len(unscored)} (Pathogenic: {sum(v['label']=='pathogenic' for v in unscored)}, Benign: {sum(v['label']=='benign' for v in unscored)})")
for i, v in enumerate(unscored, 1):
    print(f"  {i:2d}. [{v['label']}] {v['disp']}")

if len(unscored) == 11:
    path_list = [v["disp"] for v in unscored if v["label"] == "pathogenic"]
    ben_list = [v["disp"] for v in unscored if v["label"] == "benign"]
    replacement_str = (
        f"pathogenic: {', '.join(path_list)}; "
        f"benign: {', '.join(ben_list)}"
    )
    docx_files = sorted(glob.glob("NC-DriverEngine_v3.3*.docx"), key=os.path.getmtime, reverse=True)
    for doc_fp in docx_files:
        doc = docx.Document(doc_fp)
        updated = False
        for p in doc.paragraphs:
            if "[variant identifiers to be supplied]" in p.text:
                new_text = p.text.replace("[variant identifiers to be supplied]", replacement_str)
                fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
                fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else None
                for r in list(p.runs):
                    p._p.remove(r._r)
                nr = p.add_run(new_text)
                nr.font.name = fn
                if fs:
                    nr.font.size = fs
                nr.font.highlight_color = None
                updated = True
        if updated:
            doc.save(doc_fp)
            print(f"\nSuccessfully inserted all 11 variant IDs and cleared yellow highlight in: {doc_fp}")
            break
