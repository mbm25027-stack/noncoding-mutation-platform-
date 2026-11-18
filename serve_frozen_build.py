#!/usr/bin/env python
"""
serve_frozen_build.py

A small, always-on, READ-ONLY scoring server around your frozen build. It keeps the build
frozen while making it "active": it only scores, it never trains or writes into the build.

  - Verifies the build against a sha256 manifest at start-up and refuses to serve if a file changed
    (override with --allow-drift, which is logged and shown on /health).
  - Endpoints:  GET /health      status, manifest hash, uptime, request count
                POST /score      {"variants":[{"chrom":"1","pos":123,"ref":"A","alt":"G"}, ...]}
  - Every request is appended to logs/served_requests.jsonl (what was scored, when).
  - Uses only the Python standard library (no installs). Binds to 127.0.0.1 by default.

You must tell it how to call your model: --entry package.module:function
  The function receives a list of dicts (chrom, pos, ref, alt, plus any extra fields you send) and
  returns a list of numbers or dicts of the same length. Example wrapper file:  my_entry.py

      def score(variants):
          return [{"score": run_my_frozen_model(v)} for v in variants]

Examples
  Create the manifest once (hashes every file in the build folder):
      python serve_frozen_build.py --model-dir trained_server_model --make-manifest
  Start the server:
      python serve_frozen_build.py --model-dir trained_server_model --entry my_entry:score --port 8000
  Test the plumbing with a dummy scorer (not your model):
      python serve_frozen_build.py --demo --port 8000
  Query it:
      curl -X POST http://127.0.0.1:8000/score -H "Content-Type: application/json" ^
           -d "{\\"variants\\":[{\\"chrom\\":\\"1\\",\\"pos\\":12345,\\"ref\\":\\"A\\",\\"alt\\":\\"G\\"}]}"
"""
import argparse
import hashlib
import importlib
import json
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

STATE = dict(started=time.time(), requests=0, variants=0, drift=[], manifest_sha=None, model_dir=None, demo=False)
LOCK = threading.Lock()
SCORER = None
MAX_BATCH = 5000


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def build_files(model_dir):
    return sorted(p for p in model_dir.rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts and p.name != "FROZEN_MANIFEST.sha256")


def make_manifest(model_dir):
    lines = [f"{sha256(p)}  {p.relative_to(model_dir).as_posix()}" for p in build_files(model_dir)]
    out = model_dir / "FROZEN_MANIFEST.sha256"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out} with {len(lines)} entries. Keep a copy of it elsewhere; it is the proof of the freeze.")


def verify_manifest(model_dir):
    mf = model_dir / "FROZEN_MANIFEST.sha256"
    if not mf.exists():
        return None, ["no FROZEN_MANIFEST.sha256 (run with --make-manifest first)"]
    expected = {}
    for line in mf.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2:
            expected[parts[1].strip()] = parts[0].lower()
    problems = []
    actual = {p.relative_to(model_dir).as_posix(): p for p in build_files(model_dir)}
    for rel, digest in expected.items():
        p = actual.get(rel)
        if p is None:
            problems.append(f"missing: {rel}")
        elif sha256(p) != digest:
            problems.append(f"changed: {rel}")
    for rel in actual:
        if rel not in expected:
            problems.append(f"new file not in manifest: {rel}")
    return hashlib.sha256(mf.read_bytes()).hexdigest(), problems


def load_scorer(entry, demo):
    if demo:
        def demo_score(variants):
            return [{"score": round((abs(hash((v.get('chrom'), v.get('pos'), v.get('ref'), v.get('alt')))) % 1000) / 1000, 3),
                     "note": "DEMO scorer, not a real model"} for v in variants]
        return demo_score
    if not entry or ":" not in entry:
        sys.exit("Give --entry module:function (or --demo).")
    mod, fn = entry.split(":", 1)
    sys.path.insert(0, str(Path.cwd()))
    return getattr(importlib.import_module(mod), fn)


class Handler(BaseHTTPRequestHandler):
    server_version = "FrozenBuild/1.0"

    def _send(self, code, obj):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *a):
        pass

    def do_GET(self):
        if self.path.rstrip("/") in ("", "/health", "health"):
            with LOCK:
                self._send(200, dict(status="ok" if not STATE["drift"] or STATE.get("allow_drift") else "drifted",
                                     frozen=not STATE["drift"], drift=STATE["drift"][:20],
                                     manifest_sha256=STATE["manifest_sha"], model_dir=STATE["model_dir"],
                                     demo=STATE["demo"], uptime_s=round(time.time() - STATE["started"]),
                                     requests=STATE["requests"], variants_scored=STATE["variants"]))
        else:
            self._send(404, dict(error="use GET /health or POST /score"))

    def do_POST(self):
        if self.path.rstrip("/") != "/score":
            return self._send(404, dict(error="use POST /score"))
        try:
            n = int(self.headers.get("Content-Length", 0))
            data = json.loads(self.rfile.read(n) or b"{}")
            vs = data.get("variants")
            if not isinstance(vs, list) or not vs:
                return self._send(400, dict(error="send {'variants': [ {chrom,pos,ref,alt}, ... ]}"))
            if len(vs) > MAX_BATCH:
                return self._send(413, dict(error=f"max {MAX_BATCH} variants per request"))
            for v in vs:
                if not all(k in v for k in ("chrom", "pos", "ref", "alt")):
                    return self._send(400, dict(error="each variant needs chrom, pos, ref, alt"))
            res = SCORER(vs)
            if len(res) != len(vs):
                return self._send(500, dict(error="scorer returned a different number of results"))
            with LOCK:
                STATE["requests"] += 1
                STATE["variants"] += len(vs)
            Path("logs").mkdir(exist_ok=True)
            with open("logs/served_requests.jsonl", "a", encoding="utf-8") as fh:
                fh.write(json.dumps(dict(t=time.strftime("%Y-%m-%dT%H:%M:%S"), n=len(vs), client=self.client_address[0],
                                         manifest=STATE["manifest_sha"])) + "\n")
            self._send(200, dict(results=res, n=len(res), manifest_sha256=STATE["manifest_sha"]))
        except Exception as e:  # keep the server alive
            self._send(500, dict(error=f"{type(e).__name__}: {e}"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", help="folder holding the frozen build")
    ap.add_argument("--entry", help="module:function that scores a list of variant dicts")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--make-manifest", action="store_true", help="hash the build folder and exit")
    ap.add_argument("--allow-drift", action="store_true", help="serve even if files differ from the manifest")
    ap.add_argument("--demo", action="store_true", help="dummy scorer to test the server only")
    a = ap.parse_args()

    if a.make_manifest:
        if not a.model_dir:
            sys.exit("--make-manifest needs --model-dir")
        return make_manifest(Path(a.model_dir))

    STATE["demo"] = a.demo
    STATE["allow_drift"] = a.allow_drift
    if not a.demo:
        if not a.model_dir:
            sys.exit("Give --model-dir (or --demo).")
        md = Path(a.model_dir)
        STATE["model_dir"] = str(md.resolve())
        sha, problems = verify_manifest(md)
        STATE["manifest_sha"] = sha
        STATE["drift"] = problems
        if problems:
            print("BUILD DOES NOT MATCH ITS FREEZE MANIFEST:")
            for p in problems[:20]:
                print("  -", p)
            if not a.allow_drift:
                sys.exit("Refusing to serve a build that differs from the frozen manifest. "
                         "Restore the files, or pass --allow-drift to serve anyway (it will be flagged on /health).")
        else:
            print(f"Build verified against manifest {sha[:12]}...")
    global SCORER
    SCORER = load_scorer(a.entry, a.demo)
    srv = ThreadingHTTPServer((a.host, a.port), Handler)
    print(f"Serving on http://{a.host}:{a.port}  (GET /health, POST /score).  Ctrl+C to stop.")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("stopped")


if __name__ == "__main__":
    main()
