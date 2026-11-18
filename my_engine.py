# -*- coding: utf-8 -*-
# my_engine.py -- adapter between the exp scripts and the frozen NC-MutScan v3.3 engine
import server  # loads nc_mutscan_research_suite.py once, at import time

EXPECTED_BUILD = "b276d0d44045"
if server.CORE_SHORT_SHA256 != EXPECTED_BUILD:
    raise SystemExit(f"Engine build {server.CORE_SHORT_SHA256} != frozen build {EXPECTED_BUILD}")


def score(wt, mut, meta):
    res = server.run_v33_prediction({
        "wt_seq": wt,
        "mut_seq": mut,
        "flank_bp": int(meta["center"]),
        "locus": f"{meta['chrom']}:{meta['pos']}",
        "title": f"{meta['ref']}>{meta['alt']}",
    })
    if res is None:
        return {"s_raw": 0.0, "s_unified": 0.0, "call": False, "module": "", "features": {}}

    feats = {k: v for k, v in res.items()
             if isinstance(v, (int, float)) and not isinstance(v, bool)
             and k not in ("S_raw", "S_unified")}
    return {
        "s_raw": res["S_raw"],
        "s_unified": res["S_unified"],
        "call": bool(res["is_driver"]),
        "module": res.get("module", ""),
        "features": feats,
    }
