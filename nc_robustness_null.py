#!/usr/bin/env python3
"""
nc_robustness_null.py - robustness sweeps and an empirical null for the frozen NC-MutScan engine.
Needs numpy + the standard library (matplotlib is used for PNG heatmaps only if it is importable).
NOTHING is written into server.py: parameters are patched in memory for the sweep and restored afterwards.

Three sub-commands
  discover  list numeric constants / threshold-looking lines in server.py (and other modules) and show the keys
            of the engine's result dict, so you can name the parameters to sweep and the "module" field.
  sweep     heatmaps of sensitivity / FPR / PPV over tau_splice x tau_TF, tau_splice x window, tau_TF x window
            (the third axis held at its frozen value). Frozen value = the heatmap cell that is outlined.
  null      dinucleotide-shuffled flanks (Altschul-Erickson): false-call rate per module from chance motif
            matches, compared with the real flanks, plus an empirical per-variant p-value for S_raw.

Typical use (from the project folder, venv active, $fa set):
  python nc_robustness_null.py discover --server-dir . --pool check_dev_v2\\pool_all.vcf --fasta $fa
  python nc_robustness_null.py sweep --server-dir . --pool check_dev_v2\\pool_all.vcf --fasta $fa --exclude-mc none ^
       --scores eval_dev_v2\\scores.tsv --out eval_sweep_dev
  python nc_robustness_null.py null  --server-dir . --pool check_dev_v2\\pool_all.vcf --fasta $fa --exclude-mc none ^
       --out eval_null_dev
Targets are  module:NAME  or  module:DICT[key]  (several targets for one parameter: comma separated).

IMPORTANT (honesty rule): the sweep is a SENSITIVITY ANALYSIS of the frozen engine. Never pick the best cell and
call it the new setting. Explore on the development pool (check_dev_v2); show the blind pool only with the frozen
cell outlined.
"""
import argparse, importlib, math, os, random, re, sys, time
from collections import Counter, defaultdict
import numpy as np

PATH_SIG = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}
BEN_SIG = {"benign", "likely_benign", "benign/likely_benign"}


# --------------------------------------------------------------------------- frozen call (verbatim copy)
def is_frozen_driver_call(res):   # verbatim from build_dev_flanks_and_loco.py
    cls_txt = str(res.get("classification", "")).lower()
    if any(neg in cls_txt for neg in ("benign", "neutral", "polymorphism", "passenger", "sub-threshold", "subthreshold", "non-driver", "compet")):
        return 0
    if any(pos in cls_txt for pos in ("driver", "pathogenic", "loss", "gain", "splice", "polya", "disrupt")):
        return 1
    return 0


# --------------------------------------------------------------------------- pool / fasta
def load_pool(path, exclude_mc):
    rows, skipped, excl = [], Counter(), set(exclude_mc)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
            sig = info.get("CLNSIG", "").lower()
            if sig in PATH_SIG:
                y = 1
            elif sig in BEN_SIG:
                y = 0
            else:
                skipped["not_clean_label"] += 1
                continue
            terms = {x.split("|")[-1] for x in info.get("MC", "").split(",") if x}
            if terms & excl:
                skipped["excluded_mc"] += 1
                continue
            rows.append(dict(chrom=t[0].replace("chr", ""), pos=int(t[1]), ref=t[3].upper(), alt=t[4].upper(), y=y))
    return rows, skipped


def sample_pool(rows, n_benign, seed):
    P = [r for r in rows if r["y"] == 1]
    B = [r for r in rows if r["y"] == 0]
    rng = random.Random(seed)
    if n_benign and n_benign < len(B):
        B = rng.sample(B, n_benign)
    return P + B, len([r for r in rows if r["y"] == 1]), len([r for r in rows if r["y"] == 0])


class Fasta:
    def __init__(self, fa):
        self.f = open(fa, "rb")
        self.idx = {}
        with open(fa + ".fai") as g:
            for line in g:
                t = line.split("\t")
                self.idx[t[0]] = (int(t[1]), int(t[2]), int(t[3]), int(t[4]))

    def name(self, chrom):
        for c in (chrom, "chr" + chrom, "chrM" if chrom == "MT" else chrom):
            if c in self.idx:
                return c
        return None

    def fetch(self, chrom, start1, end1):
        """1-based inclusive; returns upper-case string or None when out of range."""
        nm = self.name(chrom)
        if nm is None:
            return None
        L, off, lb, lw = self.idx[nm]
        if start1 < 1 or end1 > L:
            return None
        s0, e0 = start1 - 1, end1
        a = off + (s0 // lb) * lw + s0 % lb
        b = off + ((e0 - 1) // lb) * lw + (e0 - 1) % lb + 1
        self.f.seek(a)
        return self.f.read(b - a).decode().replace("\n", "").replace("\r", "").upper()


def build_contexts(rows, fa, maxw):
    """left/right flanks of length maxw (reference) per variant; drops variants near contig ends / ref mismatches."""
    keep, drop = [], Counter()
    for r in rows:
        ref = fa.fetch(r["chrom"], r["pos"], r["pos"] + len(r["ref"]) - 1)
        if ref is None:
            drop["no_contig_or_range"] += 1
            continue
        if ref != r["ref"]:
            drop["ref_mismatch"] += 1
            continue
        L = fa.fetch(r["chrom"], r["pos"] - maxw, r["pos"] - 1)
        R = fa.fetch(r["chrom"], r["pos"] + len(r["ref"]), r["pos"] + len(r["ref"]) + maxw - 1)
        if L is None or R is None or len(L) != maxw or len(R) != maxw:
            drop["flank_out_of_range"] += 1
            continue
        r["L"], r["R"] = L, R
        keep.append(r)
    return keep, drop


def seqs_for_window(r, w, L=None, R=None):
    L = r["L"] if L is None else L
    R = r["R"] if R is None else R
    lw, rw = L[len(L) - w:], R[:w]
    return lw + r["ref"] + rw, lw + r["alt"] + rw


# --------------------------------------------------------------------------- engine
def load_server(server_dir):
    sys.path.insert(0, os.path.abspath(server_dir))
    return importlib.import_module("server")


def run_engine(server, wt, mut, w):
    try:
        res = server.run_v33_prediction({"wt_seq": wt, "mut_seq": mut, "flank_bp": w})
    except Exception:
        return None
    return res if isinstance(res, dict) else None


def clear_caches(mods):
    for m in mods:
        for v in list(vars(m).values()):
            cc = getattr(v, "cache_clear", None)
            if callable(cc):
                try:
                    cc()
                except Exception:
                    pass


class Param:
    """One tunable engine parameter, possibly stored in several places. Patched in memory only."""

    def __init__(self, spec, label):
        self.label, self.targets = label, []
        for s in spec.split(","):
            mod, _, rest = s.strip().partition(":")
            m = re.match(r"^([A-Za-z_]\w*)(?:\[(.+)\])?$", rest)
            if not mod or not m:
                raise SystemExit(f"Bad target '{s}'. Use module:NAME or module:DICT[key].")
            module = importlib.import_module(mod)
            name, key = m.group(1), m.group(2)
            if not hasattr(module, name):
                raise SystemExit(f"{mod} has no attribute {name}. Run 'discover' to list candidates.")
            if key is not None:
                key = key.strip("'\"")
            self.targets.append((module, name, key))
        self.orig = [self._get(t) for t in self.targets]
        d = self.orig[0]
        if not isinstance(d, (int, float)) or isinstance(d, bool):
            raise SystemExit(f"{label}: target {spec} is {type(d).__name__} ({d!r}); only a single number can be swept.")
        self.default = d

    def _get(self, t):
        module, name, key = t
        v = getattr(module, name)
        return v[key] if key is not None else v

    def set(self, value):
        for module, name, key in self.targets:
            if key is not None:
                getattr(module, name)[key] = value
            else:
                setattr(module, name, value)

    def restore(self):
        for t, v in zip(self.targets, self.orig):
            module, name, key = t
            if key is not None:
                getattr(module, name)[key] = v
            else:
                setattr(module, name, v)

    def modules(self):
        return [t[0] for t in self.targets]


# --------------------------------------------------------------------------- statistics
def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


# --------------------------------------------------------------------------- discover
def cmd_discover(a):
    pat = re.compile(a.pattern, re.I)
    server = load_server(a.server_dir)
    mods = [server] + [importlib.import_module(m) for m in a.extra_modules.split(",") if m]
    import inspect
    for m in mods:
        print(f"\n=== module {m.__name__}  ({getattr(m, '__file__', '?')})")
        shown = 0
        for n, v in vars(m).items():
            if n.startswith("__"):
                continue
            numeric = (isinstance(v, (int, float)) and not isinstance(v, bool)) or \
                      (isinstance(v, (tuple, list)) and 0 < len(v) <= 8 and all(isinstance(x, (int, float)) for x in v)) or \
                      (isinstance(v, dict) and 0 < len(v) <= 30 and all(isinstance(x, (int, float)) for x in v.values()))
            if numeric and (a.all or pat.search(n)):
                s = repr(v)
                print(f"  {n} = {s[:160]}")
                shown += 1
        print(f"  ({shown} numeric constants shown; --all lists every numeric constant)")
        src = getattr(m, "__file__", None)
        if src and src.endswith(".py"):
            print("  source lines that look like thresholds / windows (first 60):")
            c = 0
            with open(src, encoding="utf-8", errors="replace") as f:
                for i, line in enumerate(f, 1):
                    if pat.search(line) and re.search(r"\d", line) and not line.lstrip().startswith("#"):
                        print(f"    {i}: {line.rstrip()[:150]}")
                        c += 1
                        if c >= 60:
                            break
    if a.pool and a.fasta:
        rows, _ = load_pool(a.pool, [x for x in a.exclude_mc.split(",") if x and x != "none"])
        rows = rows[:3] + [r for r in rows if r["y"] == 0][:2]
        fa = Fasta(a.fasta)
        keep, _ = build_contexts(rows, fa, 25)
        print("\n=== result dict of the engine on a few pool variants (keys and values, truncated)")
        for r in keep[:4]:
            wt, mut = seqs_for_window(r, 25)
            res = run_engine(server, wt, mut, 25)
            print(f"  {r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']} label={r['y']}")
            if res:
                for k, v in res.items():
                    print(f"     {k}: {repr(v)[:110]}")
    if a.tally and a.pool and a.fasta:
        rows, _ = load_pool(a.pool, [x for x in a.exclude_mc.split(",") if x and x != "none"])
        samp, _, _ = sample_pool(rows, a.tally, a.seed)
        P = [r for r in samp if r["y"] == 1][:a.tally]
        samp = P + [r for r in samp if r["y"] == 0]
        keep, _ = build_contexts(samp, Fasta(a.fasta), 25)
        allc, calledc = defaultdict(Counter), defaultdict(Counter)
        ncall = 0
        print(f"\n=== tally of engine output strings over {len(keep):,} variants")
        for r in keep:
            wt, mut = seqs_for_window(r, 25)
            res = run_engine(server, wt, mut, 25)
            if not res:
                continue
            c = is_frozen_driver_call(res)
            ncall += c
            for f in ("classification", "mechanism", "primary_tf", "direction"):
                v = str(res.get(f, ""))[:110]
                allc[f][v] += 1
                if c:
                    calledc[f][v] += 1
        print(f"  frozen calls: {ncall:,} of {len(keep):,}")
        for f in ("classification", "mechanism", "primary_tf", "direction"):
            print(f"  -- {f}: values among FROZEN CALLS (top 25)")
            for v, n in calledc[f].most_common(25):
                print(f"       {n:>6}  {v}")
            print(f"  -- {f}: values among ALL variants (top 8)")
            for v, n in allc[f].most_common(8):
                print(f"       {n:>6}  {v}")
    print("\nNext: the sweep patches the threshold literals of the scoring function automatically (see 'sweep').")
    print("and the result-dict field that names the module, e.g.  --module-fields mechanism")


# --------------------------------------------------------------------------- sweep
# (key, regex with the threshold literal as group 2, description, group, required)
RULES = [
    ("sp_loss", r'(best_splice_loss\["delta"\]\s*<=\s*)(-?\d+(?:\.\d+)?)', "splice LOSS delta-PWM threshold", "splice_delta", True),
    ("sp_gain", r'(best_splice_gain\["delta"\]\s*>=\s*)(-?\d+(?:\.\d+)?)', "splice GAIN delta-PWM threshold", "splice_delta", True),
    ("sp_gain_score", r'(m_res\["score"\]\s*>=\s*)(\d+(?:\.\d+)?)', "splice gain: minimum mutant motif score", "splice_min", True),
    ("tf_gain_score", r'(mut_res\["score"\]\s*>=\s*)(\d+(?:\.\d+)?)', "TF gain: minimum mutant motif score", "tf_score", True),
    ("tf_loss_score", r'(wt_res\["score"\]\s*>=\s*)(\d+(?:\.\d+)?)', "TF loss: minimum wild-type motif score", "tf_score", True),
    ("tf_gain_delta", r'(best_gain\["delta"\]\s*>=\s*)(-?\d+(?:\.\d+)?)', "TF GAIN delta-PWM threshold", "tf_delta", False),
    ("tf_loss_delta", r'(best_loss\["delta"\]\s*(?:<=|<)\s*)(-?\d+(?:\.\d+)?)', "TF LOSS delta-PWM threshold", "tf_delta", False),
]
MINWT_RE = r'(splice_min_wt\s*=\s*\{\s*"SPLICE_5SS"\s*:\s*)(\d+(?:\.\d+)?)(\s*,\s*"SPLICE_3SS"\s*:\s*)(\d+(?:\.\d+)?)(\s*\})'


class SuitePatcher:
    """Runs a REWRITTEN IN-MEMORY COPY of analyze_variant_delta whose threshold literals are read from a dict (SWEEP).
    The files on disk are never changed (so the engine's build hash is untouched). Every rewrite is checked:
    patterns must match exactly where expected, and at the frozen setting the copy must give results identical to the original."""

    def __init__(self, module, funcname, include_min, include_tf_delta):
        self.orig_mod = module
        modname = module.__name__
        self.modname, self.funcname = modname, funcname
        self.origf = getattr(self.orig_mod, funcname)
        path = self.orig_mod.__file__
        src = open(path, encoding="utf-8-sig").read()
        self.defaults, self.groups, self.desc, self.report = {}, {}, {}, []
        new = src

        def one(pattern, what, required):
            ms = list(re.finditer(pattern, new))
            if len(ms) == 1:
                return ms[0]
            ctx = "\n".join("      " + new[max(0, new.rfind("\n", 0, m.start()) + 1): new.find("\n", m.end())].strip() for m in ms)
            if len(ms) == 0 and not required:
                self.report.append(f"  (not present) {what}")
                return None
            if len(ms) > 1 and not required:
                self.report.append(f"  (AMBIGUOUS, {len(ms)} matches, so it stays FIXED and is NOT swept) {what}")
                for ln in ctx.splitlines():
                    self.report.append("        " + ln.strip())
                return None
            raise SystemExit(f"Cannot patch '{what}': pattern matched {len(ms)} times in {path} (expected 1).\n{ctx}\n"
                             "Send this message and the lines of analyze_variant_delta around the thresholds so the rule can be fixed.")

        m = one(MINWT_RE, "splice_min_wt dictionary (5' / 3' minimum wild-type scores)", True)
        new = new[:m.start()] + m.group(1) + 'SWEEP["sp_min5"]' + m.group(3) + 'SWEEP["sp_min3"]' + m.group(5) + new[m.end():]
        for k, v, d in (("sp_min5", m.group(2), "splice loss: minimum wild-type score, 5' donor"), ("sp_min3", m.group(4), "splice loss: minimum wild-type score, 3' acceptor")):
            self.defaults[k], self.groups[k], self.desc[k] = float(v), "splice_min", d
            self.report.append(f"  {k:<14} = {v:<6} {d}")
        for key, pat, desc, grp, req in RULES:
            m = one(pat, desc, req)
            if m is None:
                continue
            new = new[:m.start()] + m.group(1) + f'SWEEP["{key}"]' + new[m.end():]
            self.defaults[key], self.groups[key], self.desc[key] = float(m.group(2)), grp, desc
            self.report.append(f"  {key:<14} = {m.group(2):<6} {desc}")
        self.sweep = dict(self.defaults)
        self.splice_groups = {"splice_delta"} | ({"splice_min"} if include_min else set())
        self.tf_groups = {"tf_score"} | ({"tf_delta"} if include_tf_delta else set())
        ns = {"__name__": modname + "__sweep", "__file__": path, "SWEEP": self.sweep}
        exec(compile(new, path + " [sweep copy]", "exec"), ns)
        if funcname not in ns:
            raise SystemExit(f"{funcname} not defined by the rewritten copy.")
        self.newf = ns[funcname]
        self.replaced = []
        # numeric comparisons on the TF decision variables that were NOT parameterised (so you can see what stays fixed)
        self.fixed_lines = []
        for i, line in enumerate(src.splitlines(), 1):
            if re.search(r'(best_gain|best_loss|best_splice_gain|best_splice_loss)\[.*?\]\s*(?:>=|<=|>|<)\s*-?\d', line) and "SWEEP[" not in line:
                if not any(re.search(pat, line) for _, pat, _, _, _ in RULES):
                    self.fixed_lines.append(f"    {i}: {line.strip()[:120]}")

    def install(self):
        self.replaced = []
        self._slots = []
        dicts = [(f"{self.modname}(engine module object)", self.orig_mod.__dict__)]
        for name, mod in list(sys.modules.items()):
            d = getattr(mod, "__dict__", None)
            if isinstance(d, dict) and d is not self.orig_mod.__dict__:
                dicts.append((name, d))
        for name, d in dicts:
            try:
                for k, v in list(d.items()):
                    if v is self.origf:
                        d[k] = self.newf
                        self._slots.append((d, k))
                        self.replaced.append((name, k))
            except Exception:
                pass

    def uninstall(self):
        for d, k in getattr(self, "_slots", []):
            d[k] = self.origf
        self.replaced, self._slots = [], []

    def apply(self, ms, mt):
        for k, g in self.groups.items():
            m = ms if g in self.splice_groups else (mt if g in self.tf_groups else 1.0)
            self.sweep[k] = self.defaults[k] * m

    def describe(self, ms, mt):
        out = []
        for k, g in self.groups.items():
            m = ms if g in self.splice_groups else (mt if g in self.tf_groups else 1.0)
            out.append(f"{k}={self.defaults[k] * m:.4g}")
        return ", ".join(out)

    def self_check(self, server, sample, w=25, n=600):
        sub = sample[:n // 2] + sample[-(n // 2):]
        self.apply(1.0, 1.0)
        keys = ("S_raw", "classification", "delta_pwm", "mechanism", "primary_tf", "direction")

        def run():
            out = []
            for r in sub:
                wt, mut = seqs_for_window(r, w)
                res = run_engine(server, wt, mut, w)
                out.append(None if res is None else tuple(str(res.get(k)) for k in keys))
            return out
        new_out = run()
        self.uninstall()
        old_out = run()
        self.install()
        n_none = sum(1 for x in old_out if x is None)
        if n_none > 0.05 * len(old_out):
            try:
                server.run_v33_prediction({"wt_seq": "ACGT" * 20, "mut_seq": "ACGA" * 20, "flank_bp": w})
                why = "no error on a test call"
            except Exception as e:
                why = f"{type(e).__name__}: {e}"
            raise SystemExit(f"The engine returned nothing for {n_none} of {len(old_out)} self-check variants ({why}); "
                             "a comparison of empty results proves nothing. Aborting.")
        diff = sum(1 for x, y in zip(new_out, old_out) if x != y)
        print(f"  SELF-CHECK: the rewritten copy at the frozen setting vs the original engine on {len(sub)} variants: {diff} differ "
              f"({'OK, identical' if diff == 0 else 'MISMATCH'})")
        if diff:
            raise SystemExit("The rewritten copy does not reproduce the original engine at the frozen setting; aborting (no sweep is trustworthy).")
        # the patch must really reach the engine: extreme thresholds have to change at least one classification
        reach = 0
        for ms, mt in ((0.1, 0.1), (3.0, 3.0)):
            self.apply(ms, mt)
            ext = run()
            reach += sum(1 for x, y in zip(ext, old_out) if x is not None and y is not None and x[1] != y[1])
        self.apply(1.0, 1.0)
        print(f"  REACH CHECK: thresholds x0.1 and x3 changed the classification of {reach} variant evaluations "
              f"({'OK, the patched function is the one the engine runs' if reach else 'NO EFFECT'})")
        if not reach:
            raise SystemExit("Changing the thresholds changed nothing: the patched function is NOT the one the engine runs "
                             "(it is probably referenced from an object this script cannot see). Aborting.")


def evaluate_point(server, rows, w, mods_to_clear):
    clear_caches(mods_to_clear)
    n = len(rows)
    calls = np.zeros(n, dtype=int)
    grp = np.zeros(n, dtype=int)   # 0 none, 1 splice, 2 polyA, 3 TF, 4 other
    bad = 0
    for i, r in enumerate(rows):
        wt, mut = seqs_for_window(r, w)
        res = run_engine(server, wt, mut, w)
        if res is None:
            bad += 1
            continue
        c = is_frozen_driver_call(res)
        calls[i] = c
        if c:
            lab = derived_module(res)
            grp[i] = 1 if lab.startswith("splice") else 2 if lab.startswith("polyA") else 3 if lab.startswith("TF") else 4
    return calls, grp, bad


def cmd_sweep(a):
    t0 = time.time()
    server = load_server(a.server_dir)
    rows, skipped = load_pool(a.pool, [x for x in a.exclude_mc.split(",") if x and x != "none"])
    sample, P_tot, B_tot = sample_pool(rows, a.n_benign, a.seed)
    fa = Fasta(a.fasta)
    windows = sorted({int(x) for x in a.windows.split(",")} | {a.base_window})
    sample, drop = build_contexts(sample, fa, max(windows))
    y = np.array([r["y"] for r in sample])
    nP, nB = int((y == 1).sum()), int((y == 0).sum())
    print(f"Pool {len(rows):,} (P {P_tot:,} / B {B_tot:,}); sweep sample: P {nP:,} + B {nB:,}; dropped {dict(drop)}")
    eng = getattr(server, a.engine_attr, None)
    if eng is not None and hasattr(eng, a.func):
        print(f"Engine module taken from server.{a.engine_attr} ({getattr(eng, '__file__', '?')}); it is patched in place in memory.")
        target = eng
    else:
        print(f"server.{a.engine_attr} not found: importing {a.suite_module} by name (works only if server.py uses that import).")
        target = importlib.import_module(a.suite_module)
    patcher = SuitePatcher(target, a.func, a.splice_include_min, a.tf_include_delta)
    print(f"Parameterised thresholds of {a.suite_module}.{a.func} (frozen values):")
    print("\n".join(patcher.report))
    if patcher.fixed_lines:
        print("  numeric comparisons that stay FIXED (not swept):\n" + "\n".join(patcher.fixed_lines))
    sv = sorted({float(x) for x in a.splice_factors.split(",")} | {1.0})
    tv = sorted({float(x) for x in a.tf_factors.split(",")} | {1.0})
    print(f"tau_splice = multiplier on groups {sorted(patcher.splice_groups)}: {sv}\ntau_TF = multiplier on groups {sorted(patcher.tf_groups)}: {tv}\n"
          f"windows (flank_bp) {windows} (base {a.base_window})")
    s0, t0v, w0 = 1.0, 1.0, a.base_window
    points = {}
    for s_ in sv:
        for t_ in tv:
            points[(s_, t_, w0)] = None
    for s_ in sv:
        for w_ in windows:
            points[(s_, t0v, w_)] = None
    for t_ in tv:
        for w_ in windows:
            points[(s0, t_, w_)] = None
    if a.full3d:
        for s_ in sv:
            for t_ in tv:
                for w_ in windows:
                    points[(s_, t_, w_)] = None
    print(f"{len(points)} grid points x {len(sample):,} variants; est. {len(points) * len(sample) / 1600 / 60:,.1f} min at ~1,600 variants/s")
    cached = {}
    if a.scores and os.path.exists(a.scores):
        with open(a.scores) as f:
            next(f, None)
            for line in f:
                k, _, c, st = line.rstrip("\n").split("\t")
                if st == "ok":
                    cached[k] = int(c)
    res_rows, done, nbad_total = [], 0, 0
    mods = [server, patcher.orig_mod]
    patcher.install()
    print(f"  patched {len(patcher.replaced)} reference(s) to {a.func} in memory: {[f'{m}.{k}' for m, k in patcher.replaced]}")
    if not patcher.replaced:
        raise SystemExit("No reference to the engine function was found to replace.")
    try:
        patcher.self_check(server, sample)
        for (s_, t_, w_) in sorted(points):
            patcher.apply(s_, t_)
            calls, grp, bad = evaluate_point(server, sample, w_, mods)
            nbad_total += bad
            kP, kB = int(calls[y == 1].sum()), int(calls[y == 0].sum())
            sens, fpr = kP / max(nP, 1), kB / max(nB, 1)
            ppv = sens * P_tot / (sens * P_tot + fpr * B_tot) if (sens * P_tot + fpr * B_tot) > 0 else float("nan")
            lo_s, hi_s = wilson(kP, nP)
            lo_f, hi_f = wilson(kB, nB)
            row = dict(tau_splice=s_, tau_tf=t_, window=w_, sens=sens, sens_lo=lo_s, sens_hi=hi_s, fpr=fpr, fpr_lo=lo_f,
                       fpr_hi=hi_f, ppv_pool=ppv, calls_P=kP, calls_B=kB, engine_errors=bad)
            for gi, gn in ((1, "splice"), (2, "polyA"), (3, "tf"), (4, "other")):
                row[f"calls_P_{gn}"] = int(((grp == gi) & (y == 1)).sum())
                row[f"calls_B_{gn}"] = int(((grp == gi) & (y == 0)).sum())
            res_rows.append(row)
            if (s_, t_, w_) == (s0, t0v, w0) and cached and w0 == 25:
                keyf = lambda r: f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}"
                mism = sum(1 for r, c in zip(sample, calls) if keyf(r) in cached and cached[keyf(r)] != c)
                nchk = sum(1 for r in sample if keyf(r) in cached)
                print(f"  REPRODUCTION CHECK at the frozen cell: {mism} of {nchk:,} calls differ from {a.scores} "
                      f"({'OK' if mism == 0 else 'MISMATCH - the frozen cell does not reproduce the cached run'})")
                if mism:
                    raise SystemExit("The frozen cell does not reproduce the cached calls; aborting.")
            done += 1
            if done % 5 == 0 or done == len(points):
                el = time.time() - t0
                print(f"  {done}/{len(points)} points | {el / 60:,.1f} min elapsed, ~{el / done * (len(points) - done) / 60:,.0f} min left")
    finally:
        patcher.apply(1.0, 1.0)
        patcher.uninstall()
    os.makedirs(a.out, exist_ok=True)
    cols = ["tau_splice", "tau_tf", "window", "sens", "sens_lo", "sens_hi", "fpr", "fpr_lo", "fpr_hi", "ppv_pool", "calls_P", "calls_B", "engine_errors",
            "calls_P_splice", "calls_P_polyA", "calls_P_tf", "calls_P_other", "calls_B_splice", "calls_B_polyA", "calls_B_tf", "calls_B_other"]
    with open(os.path.join(a.out, "sweep_grid.tsv"), "w", encoding="utf-8") as o:
        o.write("\t".join(cols) + "\n")
        for r in res_rows:
            o.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    lookup = {(r["tau_splice"], r["tau_tf"], r["window"]): r for r in res_rows}
    txt = [f"ROBUSTNESS SWEEP  pool={a.pool}  sample P {nP:,} + B {nB:,} (pool P {P_tot:,} / B {B_tot:,}); seed {a.seed}",
           f"frozen cell: tau_splice x{s0:g}, tau_TF x{t0v:g}, window={w0}   ( [ ] marks it ); tau values are MULTIPLIERS of the frozen thresholds",
           "tau_splice scales: " + ", ".join(k for k, g in patcher.groups.items() if g in patcher.splice_groups) + "   | tau_TF scales: " + ", ".join(k for k, g in patcher.groups.items() if g in patcher.tf_groups),
           "frozen values: " + patcher.describe(1.0, 1.0),
           "absolute values per level: " + " ; ".join(f"splice x{v:g}: " + patcher.describe(v, 1.0) for v in sv if v != 1.0) + " ; " + " ; ".join(f"TF x{v:g}: " + patcher.describe(1.0, v) for v in tv if v != 1.0),
           "sens = calls/P in the sample; FPR = calls/B in the (sub-sampled) benign; PPV_pool = at the pool's own prevalence using the full-pool P and B counts", ""]
    slices = [("tau_splice (rows) x tau_TF (cols), window=%d" % w0, "tau_splice", sv, "tau_tf", tv, lambda r, c: (r, c, w0)),
              ("tau_splice (rows) x window (cols), tau_TF=%g" % t0v, "tau_splice", sv, "window", windows, lambda r, c: (r, t0v, c)),
              ("tau_TF (rows) x window (cols), tau_splice=%g" % s0, "tau_tf", tv, "window", windows, lambda r, c: (s0, r, c))]
    figdata = []
    for title, rn, rv, cn, cv, fn in slices:
        for metric in ("sens", "fpr", "ppv_pool"):
            txt.append(f"--- {metric}: {title}")
            txt.append(f"{'':>12}" + "".join(f"{c:>11g}" for c in cv))
            grid = np.full((len(rv), len(cv)), np.nan)
            for i, rr in enumerate(rv):
                cells = []
                for j, cc in enumerate(cv):
                    r = lookup.get(fn(rr, cc))
                    v = r[metric] if r else float("nan")
                    grid[i, j] = v
                    mark = fn(rr, cc) == (s0, t0v, w0)
                    cells.append(f"{('[%.4f]' if mark else '%.4f') % v:>11}")
                txt.append(f"{rr:>12g}" + "".join(cells))
            txt.append("")
            figdata.append((title, metric, rn, rv, cn, cv, grid, fn))
    flat = [(r["tau_splice"], r["tau_tf"], r["window"], r["calls_P"], r["calls_B"]) for r in res_rows]
    for nm, axis_idx, vals, other in (("tau_splice", 0, sv, (None, t0v, w0)), ("tau_TF", 1, tv, (s0, None, w0)), ("window", 2, windows, (s0, t0v, None))):
        row = [x for x in flat if all(o is None or x[i] == o for i, o in enumerate(other))]
        ncalls = {(x[3], x[4]) for x in row}
        if len(row) > 1 and len(ncalls) == 1:
            txt.append(f"WARNING: varying {nm} (others frozen) changed NOTHING in the number of calls - either the sample has no variant near that threshold or the "
                       "parameter group does not reach any decision. Check the per-module call columns in sweep_grid.tsv.")
    txt += ["", "NOTES",
            "  - Explore on the development pool; on the blind pool report the grid with the frozen cell marked. Do not choose a 'best' cell.",
            "  - A robust engine shows slowly varying sens/FPR around the frozen cell. Sharp cliffs next to it mean the frozen setting sits on a knife edge.",
            "  - Intervals in sweep_grid.tsv are Wilson 95% intervals on the SAMPLE (variants treated as independent; gene clustering not applied).",
            f"  - Engine errors across the sweep: {nbad_total}. Dropped before scoring: {dict(drop)}.",
            "  - A rewritten in-memory copy of the engine function was used; no file on disk was modified (build hash untouched); at the frozen setting it was verified identical to the original.",
            "  - Thresholds scaled by a factor > 1 become STRICTER (fewer calls) for both signs: loss deltas (negative) move further from 0, gain deltas and scores move up. If motif scores have an upper limit (the engine caps the REPORTED score at 10), a large TF score multiplier can switch TF calls off completely: check calls_P_tf / calls_B_tf in sweep_grid.tsv.",
            "  - Thresholds that are not parameterised (listed at the start of the run under 'stay FIXED') were not swept."]
    open(os.path.join(a.out, "report.txt"), "w", encoding="utf-8").write("\n".join(txt) + "\n")
    print("\n".join(txt[:6]))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(3, 3, figsize=(15, 12))
        for k, (title, metric, rn, rv, cn, cv, grid, fn) in enumerate(figdata):
            ax = axs[k // 3][k % 3]
            im = ax.imshow(grid, aspect="auto", origin="lower")
            ax.set_xticks(range(len(cv)))
            ax.set_xticklabels([f"{c:g}" for c in cv], rotation=45)
            ax.set_yticks(range(len(rv)))
            ax.set_yticklabels([f"{r:g}" for r in rv])
            ax.set_xlabel(cn)
            ax.set_ylabel(rn)
            ax.set_title(f"{metric}\n{title}", fontsize=8)
            for i in range(len(rv)):
                for j in range(len(cv)):
                    if fn(rv[i], cv[j]) == (s0, t0v, w0):
                        ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, lw=2.5, ec="red"))
                    ax.text(j, i, f"{grid[i, j]:.3f}", ha="center", va="center", fontsize=6, color="white")
            fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "heatmaps.png"), dpi=150)
        print(f"PNG heatmaps -> {a.out}\\heatmaps.png")
    except ImportError:
        print("matplotlib not importable here: text heatmaps are in report.txt and numbers in sweep_grid.tsv (plot them elsewhere).")
    print(f"\nOutputs in {a.out}\\  (report.txt, sweep_grid.tsv[, heatmaps.png])  [{(time.time() - t0) / 60:,.1f} min]")


# --------------------------------------------------------------------------- dinucleotide shuffle null
def dinuc_shuffle(seq, rng):
    """Altschul & Erickson (1985): random sequence with exactly the same dinucleotide counts and the same first/last letter."""
    n = len(seq)
    if n <= 3:
        return seq
    last = seq[-1]
    edges = defaultdict(list)
    for x, y in zip(seq, seq[1:]):
        edges[x].append(y)
    while True:
        lastedge = {v: rng.choice(nx) for v, nx in edges.items() if v != last}
        ok = True
        for v in lastedge:
            seen, u = set(), v
            while u != last:
                if u in seen or u not in lastedge:
                    ok = False
                    break
                seen.add(u)
                u = lastedge[u]
            if not ok:
                break
        if ok:
            break
    for v, nx in edges.items():
        lst = list(nx)
        if v != last:
            lst.remove(lastedge[v])
            rng.shuffle(lst)
            lst.append(lastedge[v])
        else:
            rng.shuffle(lst)
        edges[v] = lst
    ptr = defaultdict(int)
    out, v = [seq[0]], seq[0]
    for _ in range(n - 1):
        w = edges[v][ptr[v]]
        ptr[v] += 1
        out.append(w)
        v = w
    return "".join(out)


def derived_module(res):
    """Clean module label from primary_tf + direction (the mechanism text embeds the PWM value, so it fragments)."""
    tf = str(res.get("primary_tf", ""))
    dr = str(res.get("direction", ""))
    sense = "gain" if "\u2191" in dr else ("loss" if "\u2193" in dr else "")
    if "SPLICE_5SS" in tf:
        base = "splice_5SS (U1 donor)"
    elif "SPLICE_3SS" in tf:
        base = "splice_3SS (U2 acceptor)"
    elif "PolyA" in tf:
        base = "polyA signal"
    elif "AID" in tf:
        return "AID/SHM"
    elif tf in ("", "None"):
        return "none"
    else:
        base = "TF motif"
    return f"{base} {sense}".strip()


def modlabel(res, field):
    if field == "module":
        return derived_module(res)
    v = str(res.get(field, "")).strip().lower()
    return (v[:60] or "(empty)")


def cmd_null(a):
    t0 = time.time()
    server = load_server(a.server_dir)
    rows, skipped = load_pool(a.pool, [x for x in a.exclude_mc.split(",") if x and x != "none"])
    sample, P_tot, B_tot = sample_pool(rows, a.n_benign, a.seed)
    fa = Fasta(a.fasta)
    w = a.window
    sample, drop = build_contexts(sample, fa, w)
    y = np.array([r["y"] for r in sample])
    n, K = len(sample), a.k
    fields = [x for x in a.module_fields.split(",") if x]
    print(f"Pool {len(rows):,}; null sample P {int((y == 1).sum()):,} + B {int((y == 0).sum()):,}; window {w}; {K} shuffles each; "
          f"~{n * (K + 1):,} engine calls (~{n * (K + 1) / 1600 / 60:,.1f} min); dropped {dict(drop)}")
    rng = random.Random(a.seed)
    real_call, real_S = np.zeros(n, int), np.full(n, np.nan)
    sh_call, sh_S = np.zeros((n, K), int), np.full((n, K), np.nan)
    real_mod = {f: [None] * n for f in fields}
    sh_mod = {f: [[None] * K for _ in range(n)] for f in fields}
    bad = 0
    for i, r in enumerate(sample):
        wt, mut = seqs_for_window(r, w)
        res = run_engine(server, wt, mut, w)
        if res is None:
            bad += 1
        else:
            real_call[i] = is_frozen_driver_call(res)
            s = res.get("S_raw")
            real_S[i] = float(s) if isinstance(s, (int, float)) else np.nan
            for f in fields:
                real_mod[f][i] = modlabel(res, f) if real_call[i] else None
        for k in range(K):
            Ls, Rs = dinuc_shuffle(r["L"][len(r["L"]) - w:], rng), dinuc_shuffle(r["R"][:w], rng)
            wt, mut = seqs_for_window(r, w, Ls, Rs)
            res = run_engine(server, wt, mut, w)
            if res is None:
                bad += 1
                continue
            sh_call[i, k] = is_frozen_driver_call(res)
            s = res.get("S_raw")
            sh_S[i, k] = float(s) if isinstance(s, (int, float)) else np.nan
            for f in fields:
                sh_mod[f][i][k] = modlabel(res, f) if sh_call[i, k] else None
        if (i + 1) % 2000 == 0 or i + 1 == n:
            el = time.time() - t0
            print(f"  {i + 1:,}/{n:,} variants | {el / 60:,.1f} min elapsed, ~{el / (i + 1) * (n - i - 1) / 60:,.0f} min left")
    brng = np.random.default_rng(a.seed)
    nb = a.n_boot

    def rate_ci(per_variant_counts, mask, denom_per_variant):
        idx = np.nonzero(mask)[0]
        x = per_variant_counts[idx].astype(float)
        est = x.sum() / (len(idx) * denom_per_variant)
        bs = np.empty(nb)
        for b in range(nb):
            s = brng.integers(0, len(idx), len(idx))
            bs[b] = x[s].sum() / (len(idx) * denom_per_variant)
        return est, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))

    txt = [f"DINUCLEOTIDE-SHUFFLE EMPIRICAL NULL  pool={a.pool}  window +/-{w} bp  K={K} shuffles/variant  n={n:,}  boot {nb}  seed {a.seed}",
           "Each flank (left and right of the variant, separately) is replaced by a random sequence with the SAME dinucleotide counts; "
           "the variant allele and the flank lengths are unchanged. Calls on shuffled flanks arise from chance motif matches / composition only.", ""]
    lines = ["group\tfield\tmodule\treal_rate\treal_lo\treal_hi\tshuf_rate\tshuf_lo\tshuf_hi\tshuf_over_real"]
    for lab, m in (("benign", y == 0), ("pathogenic", y == 1)):
        rr = rate_ci(real_call, m, 1)
        ss = rate_ci(sh_call.sum(1), m, K)
        txt.append(f"OVERALL frozen-call rate, {lab} variants: real flanks {rr[0]:.2%} [{rr[1]:.2%}, {rr[2]:.2%}]  | shuffled flanks {ss[0]:.2%} [{ss[1]:.2%}, {ss[2]:.2%}]"
                   f"  | shuffled/real = {ss[0] / rr[0] if rr[0] > 0 else float('nan'):.2f}")
        lines.append(f"{lab}\tALL\tALL\t{rr[0]:.6g}\t{rr[1]:.6g}\t{rr[2]:.6g}\t{ss[0]:.6g}\t{ss[1]:.6g}\t{ss[2]:.6g}\t{ss[0] / rr[0] if rr[0] > 0 else float('nan'):.4g}")
    txt.append("")
    for f in fields:
        cnt = Counter()
        for i in range(n):
            if real_mod[f][i]:
                cnt[real_mod[f][i]] += 1
            for lab_ in sh_mod[f][i]:
                if lab_:
                    cnt[lab_] += 1
        mods = [m for m, _ in cnt.most_common(a.top)]
        txt.append(f"PER MODULE (field '{f}', top {len(mods)} labels among frozen calls): rate = calls with that label / evaluations")
        txt.append(f"  {'module':<44}{'benign real':>14}{'benign shuffled [95% CI]':>34}{'ratio':>8}{'path real':>12}{'path shuffled':>15}")
        for m in mods:
            rc = np.array([1 if real_mod[f][i] == m else 0 for i in range(n)])
            sc = np.array([sum(1 for l in sh_mod[f][i] if l == m) for i in range(n)])
            rb, rp = rate_ci(rc, y == 0, 1), rate_ci(rc, y == 1, 1)
            sb, sp = rate_ci(sc, y == 0, K), rate_ci(sc, y == 1, K)
            ratio = sb[0] / rb[0] if rb[0] > 0 else float("nan")
            txt.append(f"  {m:<44}{rb[0]:>13.3%} {sb[0]:>12.3%} [{sb[1]:.3%}, {sb[2]:.3%}]{ratio:>8.2f}{rp[0]:>12.3%}{sp[0]:>15.3%}")
            lines.append(f"benign\t{f}\t{m}\t{rb[0]:.6g}\t{rb[1]:.6g}\t{rb[2]:.6g}\t{sb[0]:.6g}\t{sb[1]:.6g}\t{sb[2]:.6g}\t{ratio:.4g}")
            lines.append(f"pathogenic\t{f}\t{m}\t{rp[0]:.6g}\t{rp[1]:.6g}\t{rp[2]:.6g}\t{sp[0]:.6g}\t{sp[1]:.6g}\t{sp[2]:.6g}\t{sp[0] / rp[0] if rp[0] > 0 else float('nan'):.4g}")
        txt.append("")
    ge = np.array([[(1 + np.sum(sh_S[i] >= real_S[i])) / (K + 1) if not np.isnan(real_S[i]) else np.nan for _ in [0]][0] for i in range(n)])
    for lab, m in (("pathogenic", y == 1), ("benign", y == 0)):
        sel = m & (real_call == 1) & ~np.isnan(ge)
        if sel.sum():
            parts = []
            for th in (0.05, 0.10):
                parts.append(f"p <= {th:.2f}: " + (f"{np.mean(ge[sel] <= th):.1%}" if 1.0 / (K + 1) <= th + 1e-12 else f"not resolvable with K={K} (needs K >= {int(math.ceil(1 / th)) - 1})"))
            txt.append(f"EMPIRICAL P (S_raw real vs shuffled flanks), {lab} variants with a frozen call (n={int(sel.sum()):,}): "
                       f"median p = {np.median(ge[sel]):.3f}; fraction with " + "; ".join(parts) +
                       f" (smallest possible p = 1/{K + 1} = {1 / (K + 1):.3f})")
    for lab, m in (("benign", y == 0), ("pathogenic", y == 1)):
        txt.append(f"S_raw > 0: {lab} real {np.nanmean(real_S[m] > 0):.1%} vs shuffled {np.nanmean(sh_S[m] > 0):.1%}")
    txt += ["", "HOW TO READ IT",
            "  - shuffled/real ratio near 1 for a module: its calls on real benign variants are explained by chance motif matches (composition), so the module's",
            "    benign false calls are not evidence of biology; ratio well below 1: real flanks trigger the module more than chance would.",
            "  - Pathogenic calls with a small empirical p (<= 0.05) exceed what shuffled flanks produce; calls with large p are indistinguishable from chance.",
            "  - Shuffling destroys dinucleotide-context beyond order-1, repeats and positional structure; it is a null for motif content, not for all biology.",
            f"  - Engine errors: {bad}. Dropped before scoring: {dict(drop)}. Rate CIs: variant-level bootstrap (all {K} shuffles of a variant move together).",
            "  - Module labels are whatever the engine returns in the field(s) above; map them to your modules (splice, TF, polyA, ...) in the manuscript."]
    os.makedirs(a.out, exist_ok=True)
    open(os.path.join(a.out, "report.txt"), "w", encoding="utf-8").write("\n".join(txt) + "\n")
    open(os.path.join(a.out, "null_summary.tsv"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    with open(os.path.join(a.out, "null_per_variant.tsv"), "w", encoding="utf-8") as o:
        o.write("key\tlabel\treal_S_raw\treal_call\tshuffled_call_frac\tshuffled_S_raw_mean\tempirical_p\n")
        for i, r in enumerate(sample):
            o.write(f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}\t{r['y']}\t{real_S[i]:.6g}\t{real_call[i]}\t{sh_call[i].mean():.4f}\t{np.nanmean(sh_S[i]):.6g}\t{ge[i]:.4g}\n")
    print("\n".join(txt))
    print(f"\nOutputs in {a.out}\\  (report.txt, null_summary.tsv, null_per_variant.tsv)  [{(time.time() - t0) / 60:,.1f} min]")


# --------------------------------------------------------------------------- main
def main():
    for st in (sys.stdout, sys.stderr):   # engine text contains non-ASCII characters (e.g. a Greek delta); Windows consoles default to cp1252
        try:
            st.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--server-dir", default=".", help="folder containing server.py (the frozen engine)")
        p.add_argument("--pool", help="pool VCF (e.g. check_dev_v2\\pool_all.vcf)")
        p.add_argument("--fasta", help="GRCh38 FASTA (with .fai)")
        p.add_argument("--exclude-mc", default="none", help="comma list of MC terms to drop; 'none' for the rebuilt _v2 pools")
        p.add_argument("--seed", type=int, default=12345)

    d = sp.add_parser("discover")
    common(d)
    d.add_argument("--extra-modules", default="", help="more modules to scan, comma separated (e.g. nc_mutscan_research_suite)")
    d.add_argument("--pattern", default=r"tau|thresh|cutoff|window|flank|splice|tf|pwm|polya|min_|max_|score")
    d.add_argument("--all", action="store_true")
    d.add_argument("--tally", type=int, default=0, help="score N pathogenic + N benign variants and list the distinct classification/mechanism strings (e.g. 3000)")
    s = sp.add_parser("sweep")
    common(s)
    s.add_argument("--suite-module", default="nc_mutscan_research_suite", help="module that holds the scoring function with the threshold literals")
    s.add_argument("--func", default="analyze_variant_delta")
    s.add_argument("--engine-attr", default="NC_ENGINE", help="attribute of server.py that holds the loaded engine module (your server.py: NC_ENGINE)")
    s.add_argument("--splice-factors", default="0.6,0.8,1.0,1.2,1.4", help="multipliers for tau_splice (1.0 always included)")
    s.add_argument("--tf-factors", default="0.6,0.8,1.0,1.2,1.4", help="multipliers for tau_TF (1.0 always included)")
    s.add_argument("--splice-include-min", action="store_true", help="tau_splice also scales the splice minimum motif scores (5.6 / 6.4 / 7.5)")
    s.add_argument("--tf-include-delta", action="store_true", help="tau_TF also scales the TF delta-PWM thresholds (if found)")
    s.add_argument("--windows", default="15,20,25,35,50", help="flank_bp values (bp each side); base window must be included")
    s.add_argument("--base-window", type=int, default=25)
    s.add_argument("--n-benign", type=int, default=20000, help="benign variants sampled (all pathogenic kept); 0 = all")
    s.add_argument("--full3d", action="store_true", help="evaluate the whole 3-D grid, not just the three 2-D slices")
    s.add_argument("--scores", help="cached scores.tsv of the same pool: checks the frozen cell reproduces its calls")
    s.add_argument("--out", default="eval_sweep")
    nl = sp.add_parser("null")
    common(nl)
    nl.add_argument("--window", type=int, default=25)
    nl.add_argument("--k", type=int, default=10, help="shuffles per variant")
    nl.add_argument("--n-benign", type=int, default=18000, help="benign variants sampled (all pathogenic kept); 0 = all")
    nl.add_argument("--n-boot", type=int, default=1000)
    nl.add_argument("--module-fields", default="module,classification", help="engine result fields used as the module label; 'module' = derived from primary_tf + direction (splice_5SS/3SS, polyA, TF motif, gain/loss)")
    nl.add_argument("--top", type=int, default=12)
    nl.add_argument("--out", default="eval_null")
    a = ap.parse_args()
    if a.cmd in ("sweep", "null") and not (a.pool and a.fasta):
        raise SystemExit("--pool and --fasta are required")
    {"discover": cmd_discover, "sweep": cmd_sweep, "null": cmd_null}[a.cmd](a)


if __name__ == "__main__":
    main()