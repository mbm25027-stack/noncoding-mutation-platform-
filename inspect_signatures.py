src = open("nc_mutscan_research_suite.py", "r", encoding="utf-8-sig", errors="ignore").read().splitlines()
for i, line in enumerate(src):
    if line.startswith("def analyze_variant_delta") or line.startswith("def analyze_fasta_entry"):
        print("=" * 70)
        print(f"nc_mutscan_research_suite.py lines {i+1}-{min(i+45, len(src))}:")
        print("\n".join(src[i:i+45]))

print("=" * 70)
print("index.html lines 1740-1795:")
html = open("index.html", "r", encoding="utf-8-sig", errors="ignore").read().splitlines()
for idx in range(1739, min(1795, len(html))):
    print(f"{idx+1}: {html[idx]}")