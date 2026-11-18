#!/usr/bin/env python
"""
inspect_frozen_build.py  (read-only: changes nothing)

Finds your "frozen" build and collects evidence for WHY it is frozen:
  - candidate build folders (names like frozen / deploy / trained / server / holdout / bundle)
  - model files, manifests / checksums, lock files, read-only flags
  - text mentioning freezing, "do not modify", pre-registration, hold-out, no-retrain
  - git state (last commit, tags, uncommitted changes)
  - sha256 verification against any manifest it can parse
  - which local ports already have a server listening

It cannot read minds: if no file states a reason, it says so and shows the evidence it did find.

Run from the project folder:  python inspect_frozen_build.py
Optional:                      python inspect_frozen_build.py --root "D:\\noncoding mutation finding surver"
Output:                        printed + results/frozen_build_inspection.txt
"""
import argparse
import hashlib
import os
import re
import socket
import subprocess
import sys
import time
from pathlib import Path

NAME_RE = re.compile(r"frozen|freeze|deploy|trained|server|holdout|hold_out|release|bundle|final|model", re.I)
MODEL_EXT = {".pt", ".pth", ".pkl", ".pickle", ".joblib", ".onnx", ".h5", ".keras", ".ckpt", ".safetensors",
             ".bin", ".npz", ".json", ".yaml", ".yml", ".cfg", ".toml", ".ini"}
TEXT_EXT = {".md", ".txt", ".json", ".yaml", ".yml", ".py", ".cfg", ".toml", ".ini", ".sh", ".ps1", ".csv", ".log"}
KEYWORDS = [
    (r"\bfrozen\b|\bfreeze\b|\bfreezing\b", "freeze"),
    (r"do not (modify|edit|change|retrain|touch)|don't (modify|edit|change|retrain)", "do-not-modify"),
    (r"\block(ed)?\b|\bimmutable\b|read[- ]only", "locked"),
    (r"pre-?regist|preregist|prespecif|pre-specif", "pre-registered"),
    (r"hold-?out|held[- ]out|blind(ed)?\b|\bunseen\b", "holdout/blind"),
    (r"no (further )?(retrain|tuning|re-?fit)|without (retraining|tuning)", "no-retrain"),
    (r"sha-?256|checksum|md5", "checksum"),
    (r"leakage|leak\b|contaminat", "leakage"),
    (r"\bv\d+(\.\d+)*\b.*(final|release|frozen)|(final|release|frozen).*\bv\d+", "version-tag"),
]
KW = [(re.compile(p, re.I), lab) for p, lab in KEYWORDS]
PORTS = [5000, 5001, 7860, 8000, 8001, 8008, 8080, 8081, 8501, 8888, 9000, 3000]
MAX_TEXT_BYTES = 2_000_000
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def human(n):
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or u == "TB":
            return f"{n:.1f} {u}" if u != "B" else f"{int(n)} B"
        n /= 1024


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def walk_files(d, limit=20000):
    n = 0
    for p in d.rglob("*"):
        if p.is_file() and "venv" not in p.parts and "__pycache__" not in p.parts and ".git" not in p.parts:
            yield p
            n += 1
            if n >= limit:
                return


def run(cmd, cwd):
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=20)
        return r.stdout.strip()
    except Exception:
        return None


def find_candidates(root):
    out = []
    for p in sorted(root.iterdir()):
        if p.is_dir() and NAME_RE.search(p.name) and p.name.lower() not in {"venv", "__pycache__", "results"}:
            out.append(p)
    return out


def scan_text(path, hits):
    try:
        if path.suffix.lower() not in TEXT_EXT or path.stat().st_size > MAX_TEXT_BYTES:
            return
        with open(path, "r", errors="ignore") as fh:
            for i, line in enumerate(fh, 1):
                for rx, lab in KW:
                    if rx.search(line):
                        hits.append((lab, path, i, line.strip()[:200]))
                if len(hits) > 4000:
                    return
    except Exception:
        pass


def parse_manifest(path):
    """Lines like '<64 hex>  relative/path' or 'relative/path <64 hex>' or 'path: hex'."""
    pairs = []
    try:
        for line in open(path, "r", errors="ignore"):
            m = re.search(r"\b([0-9a-fA-F]{64})\b", line)
            if not m:
                continue
            rest = (line[:m.start()] + line[m.end():]).strip(" \t\r\n:*,\"'=")
            if rest:
                pairs.append((rest.replace("\\", "/"), m.group(1).lower()))
    except Exception:
        pass
    return pairs


def check_ports():
    open_ports = []
    for p in PORTS:
        s = socket.socket()
        s.settimeout(0.25)
        try:
            if s.connect_ex(("127.0.0.1", p)) == 0:
                open_ports.append(p)
        except Exception:
            pass
        finally:
            s.close()
    return open_ports


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent))
    ap.add_argument("--hash-limit-mb", type=int, default=300, help="skip hashing files larger than this")
    a = ap.parse_args()
    root = Path(a.root)
    say("=" * 78)
    say("FROZEN BUILD INSPECTION (read-only)")
    say(f"Root: {root}    Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say("=" * 78)

    cands = find_candidates(root)
    say(f"\n1. Candidate build folders at the top level ({len(cands)}):")
    for c in cands:
        files = list(walk_files(c))
        total = sum(f.stat().st_size for f in files)
        newest = max((f.stat().st_mtime for f in files), default=0)
        say(f"  - {c.name:35s} files={len(files):6d} size={human(total):>10s} "
            f"newest={time.strftime('%Y-%m-%d', time.localtime(newest)) if newest else '-'}")
    if not cands:
        say("  none found. Use --root to point at the folder that holds your build.")

    say("\n2. Model files, manifests, locks, README per candidate:")
    manifests = []
    for c in cands:
        files = list(walk_files(c))
        models = [f for f in files if f.suffix.lower() in MODEL_EXT and f.stat().st_size > 50_000]
        marks = [f for f in files if re.search(r"manifest|checksum|sha256|\.lock$|frozen|readme|license|changelog|version",
                                               f.name, re.I)]
        ro = [f for f in files if not os.access(f, os.W_OK)]
        say(f"\n  [{c.name}]")
        for f in sorted(models, key=lambda x: -x.stat().st_size)[:8]:
            say(f"    model-like : {f.relative_to(root)}  ({human(f.stat().st_size)}, "
                f"{time.strftime('%Y-%m-%d', time.localtime(f.stat().st_mtime))})")
        for f in marks[:10]:
            say(f"    marker     : {f.relative_to(root)}")
            if re.search(r"manifest|checksum|sha256", f.name, re.I):
                manifests.append(f)
        say(f"    read-only files: {len(ro)} of {len(files)}"
            + ("  (the build is write-protected)" if files and len(ro) == len(files) else ""))

    say("\n3. Manifest / checksum verification:")
    if not manifests:
        say("  no manifest or checksum file found in the candidate folders.")
    for m in manifests:
        pairs = parse_manifest(m)
        if not pairs:
            say(f"  {m.relative_to(root)}: no sha256 entries parsed")
            continue
        ok = bad = miss = skipped = 0
        for rel, digest in pairs:
            tgt = (m.parent / rel)
            if not tgt.exists():
                tgt = root / rel
            if not tgt.exists():
                miss += 1
                continue
            if tgt.stat().st_size > a.hash_limit_mb * 1024 * 1024:
                skipped += 1
                continue
            if sha256(tgt) == digest:
                ok += 1
            else:
                bad += 1
                say(f"    CHANGED since freeze: {tgt.relative_to(root)}")
        say(f"  {m.relative_to(root)}: entries={len(pairs)} match={ok} changed={bad} missing={miss} "
            f"skipped(large)={skipped}")

    say("\n4. Text that explains or implies the freeze (first 12 per label, candidate folders + root docs):")
    hits = []
    scan_roots = list(cands) + [root]
    seen = set()
    for r in scan_roots:
        if r == root:
            files = [f for f in root.glob("*") if f.is_file()]
        else:
            files = list(walk_files(r))
        for f in files:
            if f in seen or f.name in {"inspect_frozen_build.py", "serve_frozen_build.py"}:
                continue
            seen.add(f)
            scan_text(f, hits)
    by = {}
    for lab, p, i, t in hits:
        by.setdefault(lab, []).append((p, i, t))
    if not by:
        say("  no text found that mentions freezing, locking, hold-out or checksums.")
    for lab in ("freeze", "do-not-modify", "locked", "pre-registered", "holdout/blind", "no-retrain", "leakage",
                "checksum", "version-tag"):
        for p, i, t in by.get(lab, [])[:12]:
            try:
                rel = p.relative_to(root)
            except ValueError:
                rel = p
            say(f"  [{lab}] {rel}:{i}: {t}")

    say("\n5. Git state:")
    if (root / ".git").exists() and run(["git", "--version"], root):
        say("  last commit : " + str(run(["git", "log", "-1", "--format=%h %ad %s", "--date=short"], root)))
        say("  tags        : " + str(run(["git", "tag", "--sort=-creatordate"], root) or "none").replace("\n", ", ")[:300])
        st = run(["git", "status", "--porcelain"], root)
        say(f"  uncommitted : {len(st.splitlines()) if st else 0} changed files")
    else:
        say("  not a git repository (or git not installed): no commit or tag evidence of a freeze.")

    say("\n6. Servers already listening on common local ports:")
    op = check_ports()
    say("  " + (", ".join(str(p) for p in op) if op else "none. No server is running."))

    say("\n7. WHY IS IT FROZEN? (evidence only; this script cannot read intent)")
    reasons = []
    if by.get("holdout/blind") or by.get("pre-registered"):
        reasons.append("files mention hold-out / blind / pre-registered evaluation, so the freeze most likely protects "
                       "an unbiased test of an unchanged model")
    if by.get("no-retrain") or by.get("do-not-modify"):
        reasons.append("files state not to retrain or modify, so changes are explicitly disallowed")
    if manifests:
        reasons.append("a checksum manifest exists, so the freeze is verifiable by hash")
    if by.get("leakage"):
        reasons.append("files discuss leakage, so the freeze may guard against train/test contamination")
    if by.get("freeze") and not reasons:
        reasons.append("files use the word 'frozen' but give no stated reason")
    if reasons:
        for r in reasons:
            say("  - " + r)
    else:
        say("  no file states a reason. The usual reason for freezing a model before external validation is to stop\n"
            "  any tuning against the validation data, so the reported performance is unbiased. Confirm that with\n"
            "  whoever froze it, or add a FROZEN.md that records the date, commit, metrics and reason.")

    (root / "results").mkdir(exist_ok=True)
    outp = root / "results" / "frozen_build_inspection.txt"
    outp.write_text("\n".join(OUT), encoding="utf-8")
    say(f"\nSaved: {outp}")


if __name__ == "__main__":
    main()
