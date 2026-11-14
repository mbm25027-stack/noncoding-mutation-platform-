import subprocess

src = open("nc_mutscan_research_suite.py", "r", encoding="utf-8-sig", errors="ignore").read().splitlines()
for i, line in enumerate(src):
    if "TC-01" in line or "run_regression" in line or "def main" in line:
        start = max(0, i - 5)
        print(f"=== nc_mutscan_research_suite.py lines {start+1}-{min(start+85, len(src))} ===")
        print("\n".join(src[start:start+85]))
        break

print("\n=== GIT STATUS & BRANCHES ===")
print(subprocess.check_output(["git", "status", "-s"], stderr=subprocess.STDOUT, text=True))
print(subprocess.check_output(["git", "branch", "-a"], stderr=subprocess.STDOUT, text=True))