#!/usr/bin/env python3
"""
vep_spliceai_helper.py - use the Ensembl VEP website to get SpliceAI scores for the pool,
then convert them into the key<TAB>score file that nc_spliceai_complementarity.py reads
(--spliceai-tsv). Standard library only.

STEP 1 (before VEP):   python vep_spliceai_helper.py prep  pool_for_spliceai.vcf  vep_input.vcf
   - makes a VCF whose ID column = chrom:pos:ref:alt (no "chr"), sorted by position.
STEP 2: upload vep_input.vcf at https://www.ensembl.org/Homo_sapiens/Tools/VEP (GRCh38),
        tick SpliceAI under "Predictions", run, and download the full results
        (VCF or the tab-delimited "VEP" text format both work).
STEP 3 (after VEP):    python vep_spliceai_helper.py convert  vep_results.txt  spliceai_scores.tsv
   - score = max(DS_AG, DS_AL, DS_DG, DS_DL); max taken over all rows/genes of a variant.
Then:  python nc_spliceai_complementarity.py ... --spliceai-tsv spliceai_scores.tsv --out ...
"""
import gzip, sys

ORDER = {str(i): i for i in range(1, 23)}
ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
DS = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, encoding="utf-8", errors="replace")


def prep(src, dst):
    rows = []
    with opener(src) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            c = p[0][3:] if p[0].lower().startswith("chr") else p[0]
            rows.append((ORDER.get(c, 99), c, int(p[1]), p[3], p[4]))
    rows = sorted(set(rows))
    with open(dst, "w") as o:
        o.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, c, pos, r, a in rows:
            o.write(f"{c}\t{pos}\t{c}:{pos}:{r}:{a}\t{r}\t{a}\t.\t.\t.\n")
    print(f"Wrote {dst}: {len(rows):,} variants (ID = chrom:pos:ref:alt)")


def num(x):
    try:
        return float(x)
    except ValueError:
        return None


def best(vals):
    v = [x for x in vals if x is not None]
    return max(v) if v else None


def convert(src, dst):
    out, n_rows = {}, 0
    csq_fields = None
    tab_cols = None
    with opener(src) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = line.split("Format: ")[1].rstrip('">').split("|")
                continue
            if line.startswith("#Uploaded_variation"):
                tab_cols = line.lstrip("#").split("\t")
                continue
            if line.startswith("#") or not line:
                continue
            p = line.split("\t")
            n_rows += 1
            if tab_cols:                                    # tab-delimited VEP format
                key = p[0]
                extra = dict(kv.split("=", 1) for kv in p[-1].split(";") if "=" in kv)
                v = best([num(extra.get("SpliceAI_pred_" + d, "")) for d in DS])
            elif csq_fields:                                # VCF output with CSQ
                key = p[2]
                info = dict(kv.split("=", 1) for kv in p[7].split(";") if "=" in kv)
                vals = []
                for ent in info.get("CSQ", "").split(","):
                    f_ = ent.split("|")
                    d = dict(zip(csq_fields, f_))
                    vals += [num(d.get("SpliceAI_pred_" + x, "")) for x in DS]
                v = best(vals)
            else:
                sys.exit("Could not recognise the VEP output format (need the VCF or tab 'VEP' format).")
            if v is not None:
                out[key] = max(v, out.get(key, -1.0))
    with open(dst, "w") as o:
        o.write("key\tscore\n")
        for k, v in out.items():
            o.write(f"{k}\t{v}\n")
    print(f"Read {n_rows:,} rows; wrote {len(out):,} variants with a SpliceAI score to {dst}")
    if not out:
        print("WARNING: no SpliceAI columns found. Was SpliceAI ticked in VEP, and is this the full results file?")


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("prep", "convert"):
        sys.exit(__doc__)
    (prep if sys.argv[1] == "prep" else convert)(sys.argv[2], sys.argv[3])