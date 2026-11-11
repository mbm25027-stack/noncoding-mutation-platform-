#!/usr/bin/env python3
"""
audit_claims.py - checks marketing claims in an HTML page (NC-MutScan style)
against what the page's own code actually does.

Usage:
    python audit_claims.py index.html
    python audit_claims.py index.html --json report.json

Verdicts:
    SUPPORTED     - code evidence found
    UNSUPPORTED   - claim present, no evidence in code
    CONTRADICTED  - code evidence points the other way
    MANUAL        - cannot be decided statically (needs a real run / data)

Stdlib only. Static analysis, so it is heuristic: read the evidence lines.
"""
import re
import sys
import json
from pathlib import Path

# ------------------------------------------------------------------ helpers
def load(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    return text, text.splitlines()

def find_lines(lines, pattern, flags=re.I, limit=8):
    rx = re.compile(pattern, flags)
    out = []
    for i, ln in enumerate(lines, 1):
        if rx.search(ln):
            out.append(f"L{i}: {ln.strip()[:140]}")
            if len(out) >= limit:
                break
    return out

def revcomp(s):
    return s.upper()[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def scripts_only(text):
    return "\n".join(re.findall(r"<script[^>]*>(.*?)</script>", text, re.S | re.I))

def external_scripts(text, base_dir):
    """Read <script src="..."> files that exist next to the HTML."""
    out, found, missing = [], [], []
    for src in re.findall(r"<script[^>]*\bsrc\s*=\s*[\"']([^\"']+)[\"']", text, re.I):
        if re.match(r"https?://|//", src):
            missing.append(f"{src} (remote CDN, not read)")
            continue
        p = (base_dir / src.split("?")[0]).resolve()
        if p.exists():
            out.append(p.read_text(encoding="utf-8", errors="replace"))
            found.append(str(p))
        else:
            missing.append(f"{src} (not found at {p})")
    return "\n".join(out), found, missing

def blank_blocks(text, tags=("style", "script")):
    """Blank out <style>/<script> bodies but keep line numbers intact."""
    for t in tags:
        text = re.sub(rf"(<{t}[^>]*>)(.*?)(</{t}>)",
                      lambda m: m.group(1) + "\n" * m.group(2).count("\n") + m.group(3),
                      text, flags=re.S | re.I)
    return text

# ------------------------------------------------------------------ checks
def check_claim_text(lines):
    """Which headline claims actually appear in the page text."""
    pats = {
        "100% precision":        r"100\s*%\s*(targeted|hotspot|precision)|targeted hotspot precision",
        "Empirically validated": r"empirically\s+validated|6\s*/\s*6",
        "100% sensitivity":      r"100(\.0)?\s*%\s*sensitivity",
        "CADD/GWAVA/FunSeq2":    r"CADD|GWAVA|FunSeq2?",
        "Live chromatin":        r"live\s+chromatin|ATAC-?seq|H3K27ac",
        "JASPAR/HOCOMOCO/800":   r"JASPAR|HOCOMOCO|800\s*\+?\s*TF",
        "Multithreaded/<15ms":   r"multi-?thread|sub-?15\s*ms|under\s*15\s*ms|<\s*15\s*ms",
        "REST API/benchmarks":   r"REST API|benchmark scripts|source repositor",
        "Luciferase/ChIP-seq":   r"luciferase|ChIP-?seq",
    }
    return {k: find_lines(lines, v) for k, v in pats.items()}

def check_live_data(js_lines):
    hits = find_lines(js_lines, r"\bfetch\s*\(|XMLHttpRequest|axios|WebSocket|\$\.ajax", limit=15)
    urls = find_lines(js_lines, r"https?://[^\s\"'`]*(jaspar|hocomoco|encode|roadmap|ucsc|ensembl|ncbi)[^\s\"'`]*", limit=10)
    if hits:
        return "SUPPORTED", hits + urls
    return "UNSUPPORTED", ["No fetch/XHR/WebSocket calls in any <script>. Chromatin/PWM data cannot be 'live'."]

def check_tissue_prior(js_lines):
    ev = find_lines(js_lines, r"(tissue|chromatin|accessib|atac|h3k27)\w*\s*[:=]\s*[\{\[]|\b(weight|prior)s?\b\s*[:=]", limit=10)
    if not ev:
        return "MANUAL", ["Could not locate a weight table. Check the external JS or search for the tissue names by hand."]
    return "CONTRADICTED", ev + ["=> weights are stored in code, so 'tissue prior' is accurate and 'live' is not."]

def check_threading(js_lines):
    pos = find_lines(js_lines, r"new\s+(Shared)?Worker\s*\(|SharedArrayBuffer|hardwareConcurrency|Atomics\.")
    if pos:
        return "SUPPORTED", pos
    delays = find_lines(js_lines, r"setTimeout\s*\([^)]*,\s*\d{2,}")
    return "CONTRADICTED", ["No Worker/SharedArrayBuffer: single-threaded JS."] + delays

def check_latency(js_lines):
    perf = find_lines(js_lines, r"performance\.now|console\.time")
    delays = find_lines(js_lines, r"setTimeout\s*\([^)]*,\s*\d{2,}")
    ev = []
    if not perf:
        ev.append("No performance.now()/console.time(): latency is never measured.")
    ev += perf + delays
    verdict = "CONTRADICTED" if any(re.search(r",\s*(1[5-9]|[2-9])\d{1,2}\b", d) for d in delays) else "MANUAL"
    return verdict, ev

def check_matrices(text, js_lines):
    """Count PWM-like objects and flag placeholders."""
    # JASPAR-ish IDs and matrix objects with A/C/G/T arrays
    ids = set(re.findall(r"\bMA\d{4}\.\d\b", text))
    acgt_objs = re.findall(r"\{\s*[\"']?A[\"']?\s*:\s*\[[^\]]+\]\s*,\s*[\"']?C[\"']?\s*:", text)
    named = find_lines(js_lines, r"\b(pwm|matrix|matrices|motif)s?\b\s*[:=]\s*[\{\[]", limit=10)
    placeholder = find_lines(js_lines, r"placeholder|todo|stub|neutral|0\.25\s*,\s*0\.25\s*,\s*0\.25", limit=10)
    tf_hits = {tf: bool(re.search(tf, text, re.I)) for tf in
               ["CTCF", r"NF-?[kκ]B|NFKB", "E2F1", "ETS", "GABPA", "TERT"]}
    n = max(len(ids), len(acgt_objs))
    ev = [f"JASPAR IDs referenced: {len(ids)}", f"A/C/G/T matrix literals: {len(acgt_objs)}",
          f"TF names present: {tf_hits}"] + named + ["Placeholder-ish lines:"] + placeholder
    if n == 0:
        return "MANUAL", ["Found 0 matrices in the scripts read. If the engine is in a separate .js file, "
                          "pass it: python audit_claims.py index.html engine.js"] + ev
    verdict = "SUPPORTED" if n >= 800 else "CONTRADICTED"
    ev.insert(0, f"Counted ~{n} matrices vs claim of 800+ (verify the count by hand).")
    return verdict, ev

def check_samples(js):
    """Do curated sequences contain the motifs they are labelled with?"""
    seqs = re.findall(r"[\"'`]([ACGTacgtNn]{20,})[\"'`]", js)
    motifs = {"CCGGAA (ETS, TERT)": "CCGGAA", "GGAA core": "GGAA"}
    rows = []
    for s in seqs:
        u = s.upper()
        row = {"len": len(u), "head": u[:24] + "..."}
        for name, m in motifs.items():
            row[name] = (m in u) or (revcomp(m) in u)
        rows.append(row)
    if not rows:
        return "MANUAL", ["No literal ACGT strings >=20 nt found; samples may be built dynamically."]
    missing = [r for r in rows if not r["CCGGAA (ETS, TERT)"]]
    ev = [json.dumps(r) for r in rows]
    ev.append(f"{len(missing)}/{len(rows)} sequences lack CCGGAA/TTCCGG (either strand).")
    return ("CONTRADICTED" if missing else "SUPPORTED"), ev

def check_hardcoded_numbers(js_lines, all_lines):
    """Headline numbers that appear as literals, i.e. stored not computed."""
    nums = {"ΔPWM +6.42": r"6\.42", "CADD 64.2": r"64\.2", "GWAVA 58.5": r"58\.5",
            "FunSeq2 61.0": r"61\.0", "100% sens/prec": r"100(\.0)?\s*%"}
    out = {}
    for k, p in nums.items():
        hits = find_lines(all_lines, rf"(?<![\d.]){p}(?![\d])", limit=5)
        out[k] = hits or ["not found as a literal"]
    return out

def check_benchmark_data(text):
    ev = find_lines(text.splitlines(), r"auc|auroc|roc|confusion|true[_ ]?positive|fn\b|fp\b|ground[_ ]?truth|test[_ ]?set", limit=8)
    if ev:
        return "MANUAL", ev + ["Computation code exists; run it and compare to the displayed numbers."]
    return "UNSUPPORTED", ["No ground-truth labels, ROC/AUC, or confusion-matrix code. "
                           "Competitor percentages cannot be derived from this file."]

def check_links(text):
    hrefs = re.findall(r"<a\b[^>]*href\s*=\s*[\"']([^\"']*)[\"'][^>]*>(.*?)</a>", text, re.S | re.I)
    rows = []
    for h, label in hrefs:
        label = re.sub(r"<[^>]+>", "", label).strip()[:40]
        dead = h.strip() in ("", "#") or h.startswith("javascript:")
        rows.append(f"{'DEAD ' if dead else 'ok   '} {h[:60]!r}  [{label}]")
    wanted = bool(re.search(r"github\.com|/api|swagger|openapi|docs", text, re.I))
    ev = rows[:20] or ["No <a href> links at all."]
    return ("MANUAL" if wanted else "UNSUPPORTED"), ev

# ------------------------------------------------------------------ main
def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    path = sys.argv[1]
    text, lines = load(path)
    base = Path(path).resolve().parent
    ext_js, ext_found, ext_missing = external_scripts(text, base)
    # also pick up sibling .js files (e.g. engine.js) that are imported as modules
    for extra in sys.argv[2:]:
        if extra.endswith(".js") and Path(extra).exists():
            ext_js += "\n" + Path(extra).read_text(encoding="utf-8", errors="replace")
            ext_found.append(extra)
    js = scripts_only(text) + "\n" + ext_js
    js_lines = js.splitlines()
    page_lines = blank_blocks(text).splitlines()   # visible HTML only, line numbers kept
    print("External scripts read:", ext_found or "none")
    for m in ext_missing:
        print("  skipped:", m)

    report = {"file": path, "claims_in_text": check_claim_text(page_lines), "checks": {}}
    lines = page_lines + js_lines
    checks = {
        "Live chromatin (ATAC/H3K27ac) fetched":      check_live_data(js_lines),
        "Chromatin = hardcoded tissue weights":       check_tissue_prior(js_lines),
        "Multithreaded":                              check_threading(js_lines),
        "Sub-15ms latency":                           check_latency(js_lines),
        "800+ TF profiles (JASPAR/HOCOMOCO)":         check_matrices(text, js_lines),
        "Sample seqs contain labelled motifs":        check_samples(js),
        "CADD/GWAVA/FunSeq2 benchmark computed":      check_benchmark_data(text),
        "Docs / API / repo links exist":              check_links(text),
    }
    for k, (v, ev) in checks.items():
        report["checks"][k] = {"verdict": v, "evidence": ev}
    report["hardcoded_numbers"] = check_hardcoded_numbers(js_lines, lines)

    # ----- print
    print(f"\n=== AUDIT: {path} ===\n")
    print("Claims present in page text:")
    for k, v in report["claims_in_text"].items():
        print(f"  [{'FOUND' if v else 'absent'}] {k}")
        for ln in v[:2]:
            print(f"        {ln}")
    print("\nCode checks:")
    for k, d in report["checks"].items():
        print(f"\n  {d['verdict']:<13} {k}")
        for e in d["evidence"][:8]:
            print(f"      {e}")
    print("\nHeadline numbers as literals (stored, not computed):")
    for k, v in report["hardcoded_numbers"].items():
        print(f"  {k}")
        for e in v[:3]:
            print(f"      {e}")
    print("\nNot checkable statically: ΔPWM for TERT -124C>T. Run the real scorer on the")
    print("sequence and diff against the stored 6.42 (wire the engine's function name in).")

    if "--json" in sys.argv:
        out = sys.argv[sys.argv.index("--json") + 1]
        Path(out).write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"\nJSON written to {out}")

if __name__ == "__main__":
    main()