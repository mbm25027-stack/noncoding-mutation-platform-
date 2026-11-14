src = open("nc_mutscan_research_suite.py", "r", encoding="utf-8-sig", errors="ignore").read().splitlines()
print("=== nc_mutscan_research_suite.py lines 1-115 ===")
print("\n".join(src[:115]))

html = open("index.html", "r", encoding="utf-8-sig", errors="ignore").read().splitlines()
print("\n=== index.html lines 1800-1875 ===")
print("\n".join(html[1799:1875]))