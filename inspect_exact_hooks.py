import importlib.util, json, re

spec = importlib.util.spec_from_file_location("nc33", "nc_mutscan_research_suite.py")
nc33 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nc33)

wt = "CCCCTCCCGGGTCCCCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCC"
mut = "CCCCTCCCGGGTCCCCGGCCCAGCCCCTTCCGGGCCCTCCCAGCCCCTCCC"
res = nc33.analyze_variant_delta(wt, mut, 25, "chr5:1295113", "TERT_-124C>T")
print("=== 1. EXACT RETURN DICT OF analyze_variant_delta ===")
print(json.dumps(res, indent=2, default=str))

html = open("index.html", "r", encoding="utf-8-sig", errors="ignore").read()
print("\n=== 2. ALL '3.1' OCCURRENCES IN index.html ===")
for i, line in enumerate(html.splitlines(), 1):
    if "3.1" in line or "analyzeSequence(" in line:
        print(f"Line {i}: {line.strip()[:130]}")