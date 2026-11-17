import collections, numpy as np

PATH = r"check_blind\pool_all.vcf"
TARGET_P, TARGET_B = 3472, 193490
STARS = {"reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_single_submitter": 1}

agg = collections.Counter()
term_cls = collections.Counter()
for line in open(PATH, encoding="utf-8", errors="replace"):
    if line.startswith("#"):
        continue
    t = line.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
    cls = "P" if "pathogenic" in info.get("CLNSIG", "").lower() else "B"
    stars = STARS.get(info.get("CLNREVSTAT", ""), 0)
    snv = len(t[3]) == 1 and len(t[4]) == 1
    terms = frozenset(x.split("|")[-1] for x in info.get("MC", "").split(","))
    agg[(cls, stars, snv, terms)] += 1
    for x in terms:
        term_cls[(x, cls)] += 1

print("MC term counts (P / B):")
names = sorted({x for x, _ in term_cls}, key=lambda x: -(term_cls[(x, "P")] + term_cls[(x, "B")]))
for x in names:
    print(f"  {x:40s} P={term_cls[(x,'P')]:6d}  B={term_cls[(x,'B')]:6d}")

use = names[:15]
keys = list(agg)
masks = np.array([sum(1 << i for i, tm in enumerate(use) if tm in k[3]) for k in keys])
isP = np.array([k[0] == "P" for k in keys])
st = np.array([k[1] for k in keys]); snv = np.array([k[2] for k in keys])
cnt = np.array([agg[k] for k in keys])

res = []
for ms in range(4):
    for so in (False, True):
        sel = (st >= ms) & (snv if so else True)
        m, p, c = masks[sel], isP[sel], cnt[sel]
        for S in range(1 << len(use)):
            kept = (m & S) == 0
            P = int(c[kept & p].sum()); B = int(c[kept & ~p].sum())
            res.append((abs(P - TARGET_P) + abs(B - TARGET_B), ms, so, S, P, B))
res.sort()
print("\nClosest matches (distance, min_stars, snv_only, P, B, excluded MC terms):")
for d, ms, so, S, P, B in res[:8]:
    ex = [use[i] for i in range(len(use)) if S >> i & 1]
    print(d, ms, so, P, B, ex)