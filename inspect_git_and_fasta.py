import subprocess

print("=== GIT TRACKED FILES (git ls-files) ===")
print(subprocess.check_output(["git", "ls-files"], stderr=subprocess.STDOUT, text=True))

html = open("index.html", "r", encoding="utf-8-sig", errors="ignore").read().splitlines()
print("=== index.html lines 1240-1379 ===")
print("\n".join(html[1239:1379]))