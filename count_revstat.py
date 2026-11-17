import collections
c = collections.Counter()
for line in open(r"check_blind\pool_all.vcf", encoding="utf-8", errors="replace"):
    if line.startswith("#"):
        continue
    info = dict(kv.split("=", 1) for kv in line.split("\t")[7].strip().split(";") if "=" in kv)
    cls = "P" if "pathogenic" in info.get("CLNSIG", "").lower() else "B"
    c[(cls, info.get("CLNREVSTAT", ""))] += 1
for k, v in sorted(c.items()):
    print(k, v)