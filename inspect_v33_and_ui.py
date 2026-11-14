import ast, re

print("=== 1. PUBLIC FUNCTIONS IN FROZEN nc_mutscan_research_suite.py ===")
src = open("nc_mutscan_research_suite.py", "r", encoding="utf-8-sig", errors="ignore").read()
tree = ast.parse(src)
for node in tree.body:
    if isinstance(node, ast.FunctionDef):
        args = [a.arg for a in node.args.args]
        print(f"  def {node.name}({', '.join(args)})")

print("\n=== 2. CHECKING INDEX.HTML ENGINE / VERSION / FUNCTIONS ===")
html = open("index.html", "r", encoding="utf-8-sig", errors="ignore").read()
for pat in (r"NC-MutScan[^\n<\"']{0,60}", r"function\s+[a-zA-Z0-9_]+\s*\([^)]*\)", r"http[s]?://[^\s\"'<>]+"):
    matches = list(dict.fromkeys(re.findall(pat, html)))
    print(f"Pattern {pat}: {len(matches)} matches -> {matches[:12]}")