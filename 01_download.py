import glob
import os
import subprocess
import time
from pathlib import Path

BASE = Path(r"D:\noncoding mutation finding surver")
manifest = glob.glob(str(BASE / "gdc_manifest*"))[0]
out_dir = BASE / "maf_files"
tmp_dir = BASE / "tmp"
out_dir.mkdir(exist_ok=True)
tmp_dir.mkdir(exist_ok=True)

# Keep temp files on the same drive (D:) as the output folder
env = os.environ.copy()
env["TEMP"] = str(tmp_dir)
env["TMP"] = str(tmp_dir)
env["TMPDIR"] = str(tmp_dir)

TARGET = 472
MAX_ROUNDS = 15

def n_done():
    return len(list(out_dir.rglob("*.maf.gz")))

for attempt in range(1, MAX_ROUNDS + 1):
    have = n_done()
    print(f"Attempt {attempt}: {have}/{TARGET} files present")
    if have >= TARGET:
        break
    subprocess.run(
        [str(BASE / "gdc-client.exe"), "download",
         "-m", manifest, "-d", str(out_dir)],
        env=env,
    )
    time.sleep(30)

print(f"Finished: {n_done()}/{TARGET} files in {out_dir}")