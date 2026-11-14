"""
extract_11_from_spliceai_json.py
Reads results/clean_submission/blind_150_spliceai.json and locked_blind_chr12_22.vcf,
identifies the exact 11 variants (6 pathogenic, 5 benign) unscored by SpliceAI,
prints them, and replaces '[variant identifiers to be supplied]' (clearing yellow
highlight) in the manuscript .docx.
"""

import glob
import json
import os
import re
import docx

# 1. Load locked_blind_chr12_22.vcf
blind_variants = []
with open("locked_blind_chr12_22.vcf", "r", encoding="utf-8", errors="replace") as fh:
    for ln in fh:
        ln = ln.rstrip("\r\n")
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if parts[0].upper() in ("PATHOGENIC", "BENIGN", "LIKELY_PATHOGENIC", "LIKELY_BENIGN"):
            label = "pathogenic" if "PATHOGENIC" in parts[0].upper() else "benign"
            chrom, pos, vid, ref, alt = parts[1], parts[2], parts[3], parts[4], parts[5]
        else:
            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parts[7] if len(parts) > 7 else ln
            label = "pathogenic" if "pathogenic" in info.lower() else "benign"

        c_clean = chrom.replace("chr", "")
        disp = f"ClinVar ID {vid} (chr{c_clean}:{pos} {ref}>{alt})"
        blind_variants.append({
            "vid": str(vid).strip(),
            "c_clean": c_clean,
            "pos": str(pos).strip(),
            "ref": ref.strip().upper(),
            "alt": alt.strip().upper(),
            "key1": f"{c_clean}:{pos}:{ref}:{alt}",
            "key2": f"chr{c_clean}:{pos}:{ref}:{alt}",
            "key3": f"{c_clean}-{pos}-{ref}-{alt}",
            "disp": disp,
            "label": label,
        })

# 2. Load results\clean_submission\blind_150_spliceai.json
json_path = os.path.join("results", "clean_submission", "blind_150_spliceai.json")
with open(json_path, "r", encoding="utf-8") as jf:
    sdata = json.load(jf)

print("Loaded JSON type:", type(sdata).__name__, "| length:", len(sdata))

scored_keys = set()
scored_vids = set()

def is_valid_entry(val):
    if val is None:
        return False
    if isinstance(val, (int, float)):
        return True
    if isinstance(val, dict):
        if val.get("error") or val.get("status") == "error":
            return False
        for k in ("ds_max", "DS_MAX", "max_ds", "score", "DS_AG", "DS_AL", "DS_DG", "DS_DL"):
            if k in val and val[k] is not None and str(val[k]).lower() not in ("na", "nan", "none", "null", ""):
                return True
        # Check if any numeric score exists anywhere in dict
        for v in val.values():
            if isinstance(v, (int, float)):
                return True
    return False

if isinstance(sdata, dict):
    for k, v in sdata.items():
        if is_valid_entry(v):
            k_norm = str(k).replace("chr", "").replace("_", ":").replace("-", ":").replace(">", ":")
            scored_keys.add(k_norm)
            if isinstance(v, dict) and "clinvar_id" in v:
                scored_vids.add(str(v["clinvar_id"]))
elif isinstance(sdata, list):
    for item in sdata:
        if is_valid_entry(item):
            s = json.dumps(item)
            for field in ("id", "vid", "clinvar_id", "variation_id"):
                if field in item and item[field] is not None:
                    scored_vids.add(str(item[field]))
            c = str(item.get("chrom", item.get("chr", ""))).replace("chr", "")
            p = str(item.get("pos", item.get("position", "")))
            r = str(item.get("ref", ""))
            a = str(item.get("alt", ""))
            if c and p:
                scored_keys.add(f"{c}:{p}:{r}:{a}")
                scored_keys.add(f"{c}:{p}")
            for k_field in ("variant", "key", "locus"):
                if k_field in item and item[k_field]:
                    scored_keys.add(str(item[k_field]).replace("chr", "").replace("_", ":").replace("-", ":").replace(">", ":"))

unscored = []
for v in blind_variants:
    k_colon = f"{v['c_clean']}:{v['pos']}:{v['ref']}:{v['alt']}"
    k_cp = f"{v['c_clean']}:{v['pos']}"
    if (v["vid"] in scored_vids) or (k_colon in scored_keys) or (k_cp in scored_keys):
        continue
    unscored.append(v)

print(f"\nScored: {150 - len(unscored)} | Unscored: {len(unscored)} (Pathogenic: {sum(x['label']=='pathogenic' for x in unscored)}, Benign: {sum(x['label']=='benign' for x in unscored)})")
for i, v in enumerate(unscored, 1):
    print(f"  {i:2d}. [{v['label']}] {v['disp']}")

if len(unscored) == 11:
    path_list = [v["disp"] for v in unscored if v["label"] == "pathogenic"]
    ben_list = [v["disp"] for v in unscored if v["label"] == "benign"]
    replacement_str = (
        f"pathogenic: {', '.join(path_list)}; "
        f"benign: {', '.join(ben_list)}"
    )
    docx_files = sorted(glob.glob("*.docx"), key=os.path.getmtime, reverse=True)
    for doc_fp in docx_files:
        if doc_fp.startswith("~$"):
            continue
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
