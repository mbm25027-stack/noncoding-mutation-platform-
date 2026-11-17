# -*- coding: utf-8 -*-
"""
my_engine_template.py  --  the ONLY file you have to edit.

Copy it to  my_engine.py  (same folder as the scripts) and make score() call your
NC-MutScan v3.3 build.  The experiment scripts call it once per variant:

    score(wt, mut, meta) -> dict

    wt, mut : str   wild-type and mutant sequence, 25 bp either side of the variant
                    (51 bp for a SNV, variant at index meta["center"])
    meta    : dict  chrom, pos, ref, alt, gene, center

    returns : dict  s_raw     float   pre-threshold score (max |dPWM| over modules)
                    s_unified float   thresholded score (0 if sub-threshold or masked)
                    call      bool    the frozen binary rule (S_unified >= 4.2); optional
                    module    str     'splice' | 'polyA' | 'TF' | 'none'   (optional but useful)
                    features  dict    any extra numbers, e.g. {'splice_dpwm': 3.1, 'wt_splice': 6.2,
                                      'polya_dpwm': 0, 'tf_dpwm': 1.4}  (optional; used in exp2
                                      as extra model features)

Find the function to call in your code:
    python nc_exp_common.py probe server
(and look at how run_locked_blind.py / evaluate_v32_dev.py already call the engine).

Use the SAME build that produced the manuscript numbers (v3.3, SHA-256 b276d0d44045), not the
live web page, and run it in-process: ~2,000 variants/s locally versus ~15 ms per HTTP call.
"""

THRESHOLD = 4.2        # S_unified >= 4.2 bits -> rule-in call (manuscript definition)


def score(wt: str, mut: str, meta: dict) -> dict:
    # ---------------- OPTION A: call your engine in-process (recommended) ----------------
    # from server import analyze_pair              # <- the name you found with `probe`
    # r = analyze_pair(wt, mut)                    # adapt the arguments to its signature
    # return dict(
    #     s_raw=r["s_raw"],                        # <- adapt the keys to your return schema
    #     s_unified=r["s_unified"],
    #     call=r["s_unified"] >= THRESHOLD,
    #     module=r.get("module", ""),
    #     features={"splice_dpwm": r.get("splice_dpwm", 0.0), "tf_dpwm": r.get("tf_dpwm", 0.0)},
    # )

    # ---------------- OPTION B: local HTTP server (slow: ~15 ms/variant) -----------------
    # import requests
    # r = requests.post("http://127.0.0.1:5000/<your-endpoint>", json={"sequence": mut, ...}, timeout=30).json()
    # return dict(s_raw=..., s_unified=..., call=...)

    raise NotImplementedError("Edit my_engine.py: wire score() to your NC-MutScan v3.3 engine.")