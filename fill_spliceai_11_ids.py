"""
fill_spliceai_11_ids.py
Finds the 11 locked chr12-22 blind ClinVar variants (6 pathogenic, 5 benign)
unscored by SpliceAI from local workspace files, prints their ClinVar IDs
and coordinates, and replaces '[variant identifiers to be supplied]' (removing
yellow highlight) in the active manuscript .docx.
"""

import csv
import glob
import json
import os
import re
import docx

# 1. Load the 150 locked blind variants from locked_blind_chr12_22.vcf (or matching VCF/CSV)
vcf_candidates = glob.glob("**/locked_blind_chr12_22*.vcf", recursive=True) + glob.glob("**/locked_blind_chr12_22*.tsv", recursive=True)
if not vcf_candidates:
    raise FileNotFoundError("Could not find locked_blind_chr12_22.vcf in workspace.")

blind_variants = []  # list of dicts: id, chrom, pos, ref, alt, label
vcf_path = vcf_candidates[0]
print("Reading blind cohort:", vcf_path)

with open(vcf_path, "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\r\n").split("\t")
        chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        info = parts[7] if len(parts) > 7 else ""
        # Extract ClinVar ID / AlleleID if present in ID or INFO
        vcv_m = re.search(r"(VCV\d+|rs\d+)", vid + ";" + info)
        aid_m = re.search(r"ALLELEID=(\d+)", info)
        disp_id = vcv_m.group(1) if vcv_m else (f"ClinVar:{aid_m.group(1)}" if aid_m else vid)
        if disp_id == ".":
            disp_id = f"{chrom}:{pos}:{ref}>{alt}"
        else:
            disp_id = f"{disp_id} ({chrom}:{pos}:{ref}>{alt})"

        clnsig_m = re.search(r"CLNSIG=([^;]+)", info, re.I)
        sig_txt = clnsig_m.group(1).lower() if clnsig_m else line.lower()
        label = "pathogenic" if "pathogenic" in sig_txt else "benign"
        blind_variants.append({
            "key": f"{chrom.replace('chr','')}:{pos}:{ref}:{alt}",
            "disp": disp_id,
            "label": label,
            "scored": False,
        })

# 2. Scan workspace CSV/TSV/JSON files for SpliceAI scored variants
scored_keys = set()
for ext in ("*.csv", "*.tsv", "*.json", "*.vcf"):
    for fp in glob.glob(f"**/{ext}", recursive=True):
        if os.path.abspath(fp) == os.path.abspath(vcf_path):
            continue
        try:
            with open(fp, "r", encoding="utf-8", errors="replace") as fh:
                head = fh.read(4096)
                if "spliceai" not in head.lower() and "ds_max" not in head.lower() and "ds_ag" not in head.lower():
                    continue
                fh.seek(0)
                print("Found SpliceAI results file:", fp)
                if fp.endswith(".json"):
                    data = json.load(fh)
                    items = data if isinstance(data, list) else data.values()
                    for it in items:
                        s = json.dumps(it)
                        m = re.search(r"(?:chr)?(1[2-9]|2[0-2])[:_-](\d+)[:_-]([ACGT]+)[:_>-]([ACGT]+)", s, re.I)
                        if m and re.search(r"(?:ds_max|DS_AG|spliceai)", s, re.I):
                            if not re.search(r'"ds_max"\s*:\s*(?:null|"NA"|"")', s, re.I):
                                scored_keys.add(f"{m.group(1)}:{m.group(2)}:{m.group(3).upper()}:{m.group(4).upper()}")
                else:
                    delim = "\t" if fp.endswith((".tsv", ".vcf")) else ","
                    for row in csv.reader(fh, delimiter=delim):
                        if not row or row[0].startswith("##"):
                            continue
                        row_s = "\t".join(row)
                        m = re.search(r"(?:chr)?(1[2-9]|2[0-2])[:\t_-](\d+)[:\t_-]([ACGT]+)[:\t_>-]([ACGT]+)", row_s, re.I)
                        if m:
                            # Check that a numeric SpliceAI delta score is present on the line
                            if re.search(r"(?:^|[\t,|])0?\.\d+", row_s):
                                scored_keys.add(f"{m.group(1)}:{m.group(2)}:{m.group(3).upper()}:{m.group(4).upper()}")
        except Exception:
            pass

for v in blind_variants:
    if v["key"] in scored_keys:
        v["scored"] = True

unscored = [v for v in blind_variants if not v["scored"]]
print(f"\nTotal blind variants: {len(blind_variants)} | Scored: {len(blind_variants)-len(unscored)} | Unscored: {len(unscored)}")
for idx, v in enumerate(unscored, 1):
    print(f"  {idx:2d}. [{v['label']}] {v['disp']}")

if len(unscored) == 11:
    id_list_str = "; ".join(v["disp"] for v in unscored)
    # 3. Replace placeholder and clear yellow highlight in the manuscript .docx
    docx_files = sorted(glob.glob("NC-DriverEngine_v3.3*.docx"), key=os.path.getmtime, reverse=True)
    for doc_fp in docx_files:
        doc = docx.Document(doc_fp)
        replaced = False
        for p in doc.paragraphs:
            if "[variant identifiers to be supplied]" in p.text:
                full_t = p.text.replace("[variant identifiers to be supplied]", id_list_str)
                fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
                fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else None
                for r in list(p.runs):
                    p._p.remove(r._r)
                nr = p.add_run(full_t)
                nr.font.name = fn
                if fs:
                    nr.font.size = fs
                nr.font.highlight_color = None
                replaced = True
        if replaced:
            doc.save(doc_fp)
            print(f"\nReplaced placeholder and cleared yellow highlight in: {doc_fp}")
            break
else:
    print("\nNote: Expected 11 unscored variants; paste the SpliceAI output filename if stored elsewhere.")
