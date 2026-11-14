import re

for fn in ("Procfile", "render.yaml", "app.py"):
    print("=" * 70)
    print("FILE:", fn)
    print("=" * 70)
    print(open(fn, "r", encoding="utf-8", errors="ignore").read())

print("=" * 70)
print("FILE: server.py")
print("=" * 70)
print(open("server.py", "r", encoding="utf-8", errors="ignore").read())

print("=" * 70)
print("INDEX.HTML API / ENDPOINT REFERENCES")
print("=" * 70)
html = open("index.html", "r", encoding="utf-8", errors="ignore").read()
for m in re.finditer(r".{0,120}(onrender\.com|/predict|/api/analyze|nc_mutscan|Mutated_Donors).{0,120}", html, re.I):
    print("->", m.group(0).replace("\n", " ").strip())