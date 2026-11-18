# sv_engine.py -- rule-based screen for structural variants affecting CTCF anchors / regulatory elements
try:
    import epigenomics_engine
except Exception:
    epigenomics_engine = None

SUPPORTED = {"TRANSLOCATION", "INVERSION", "DELETION", "DUPLICATION"}
SCORE_NOTE = "Rule-based heuristic score, not a calibrated probability."


def _norm_ctx(ctx):
    ctx = ctx if isinstance(ctx, dict) else {}
    et = ctx.get("element_type", "")
    if isinstance(et, (list, tuple, set)):
        et = ",".join(str(x) for x in et)
    et = str(et or "")
    return {**ctx, "in_regulatory_region": bool(ctx.get("in_regulatory_region", False)), "element_type": et}


def _has_ctcf(ctx):
    return "CTCF" in ctx["element_type"].upper()


def analyze_structural_variant(sv_type, chrom_a, pos_a, chrom_b, pos_b):
    sv_type = str(sv_type or "").strip().upper()
    chrom_a, chrom_b = str(chrom_a or "").strip(), str(chrom_b or "").strip()
    try:
        pos_a, pos_b = int(pos_a), int(pos_b)
    except (TypeError, ValueError):
        return {"error": "pos_a and pos_b must be integers."}

    if sv_type not in SUPPORTED:
        return {"error": f"Unsupported sv_type '{sv_type}'. Use one of: {', '.join(sorted(SUPPORTED))}."}
    if not chrom_a or not chrom_b or pos_a <= 0 or pos_b <= 0:
        return {"error": "Provide chrom_a, chrom_b and positive pos_a, pos_b."}
    if sv_type != "TRANSLOCATION" and chrom_a != chrom_b:
        return {"error": f"{sv_type} breakpoints must be on the same chromosome."}
    if sv_type != "TRANSLOCATION" and pos_a >= pos_b:
        return {"error": "For intra-chromosomal events pos_a must be smaller than pos_b."}
    if epigenomics_engine is None:
        return {"error": "epigenomics_engine is not available; cannot annotate breakpoints."}

    ctx_a = _norm_ctx(epigenomics_engine.query_locus_context(chrom_a, pos_a))
    ctx_b = _norm_ctx(epigenomics_engine.query_locus_context(chrom_b, pos_b))
    bps = [(f"{chrom_a}:{pos_a}", ctx_a, f"{chrom_b}:{pos_b}"),
           (f"{chrom_b}:{pos_b}", ctx_b, f"{chrom_a}:{pos_a}")]

    is_driver = False
    verdict = "NO REGULATORY DISRUPTION DETECTED"
    mechanism = "Breakpoints do not fall on annotated CTCF anchors or regulatory elements."
    score = 0.15

    ctcf_hits = [(loc, ctx) for loc, ctx, _ in bps if _has_ctcf(ctx)]
    reg_hits = [(loc, ctx, other) for loc, ctx, other in bps if ctx["in_regulatory_region"]]

    if sv_type == "DELETION":
        if ctcf_hits:
            loc = ctcf_hits[0][0]
            is_driver, score = True, 0.95
            verdict = "CANDIDATE INSULATOR DELETION"
            mechanism = f"Deletion breakpoint at {loc} hits a CTCF boundary element, which may allow ectopic enhancer-promoter contact."
        elif reg_hits:
            loc = reg_hits[0][0]
            is_driver, score = True, 0.60
            verdict = "CANDIDATE REGULATORY ELEMENT DELETION"
            mechanism = f"Deletion breakpoint at {loc} lies in a regulatory element ({reg_hits[0][1]['element_type'] or 'unspecified'})."
    elif sv_type in ("TRANSLOCATION", "INVERSION", "DUPLICATION"):
        if ctcf_hits:
            loc = ctcf_hits[0][0]
            is_driver, score = True, 0.92
            verdict = "CANDIDATE TAD BOUNDARY DISRUPTION"
            mechanism = f"Breakpoint at {loc} falls on a CTCF anchor, which may alter 3D chromatin domains."
        elif reg_hits:
            loc, ctx, other = reg_hits[0]
            is_driver, score = True, 0.88
            verdict = "CANDIDATE ENHANCER HIJACKING"
            mechanism = (f"Rearrangement places the regulatory element at {loc} "
                         f"({ctx['element_type'] or 'unspecified'}) next to the locus at {other}.")

    return {
        "sv_type": sv_type,
        "is_driver": is_driver,
        "verdict": verdict,
        "mechanism": mechanism,
        "pathogenicity_score": round(score, 4),
        "score_note": SCORE_NOTE,
        "breakpoint_a_context": ctx_a,
        "breakpoint_b_context": ctx_b,
    }