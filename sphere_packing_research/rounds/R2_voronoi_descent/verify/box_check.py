"""Second, fully independent (rigorous box brute force, no Fincke-Pohst) check of selected saved forms incl. the top V per n."""
import numpy as np, json, glob, os, math, time
from vlib import box_bruteforce, kappa
here = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(99)
recs = {}
for f in sorted(glob.glob(os.path.join(here, '..', 'data', 'runs_n*_*.json'))):
    d = json.load(open(f))
    recs.setdefault(d['n'], []).extend(d['recs'])
out = []
t0 = time.time()
for n in sorted(recs):
    R = recs[n]
    # box size of each record
    def bsize(r):
        A = np.array(r['A_ext']); Ai = np.linalg.inv(A)
        return float(np.prod(2 * np.floor(np.sqrt((1 + 1e-6) * np.diag(Ai)) + 1e-9) + 1))
    sizes = np.array([bsize(r) for r in R])
    Vs = np.array([r['V_ext'] for r in R])
    order = np.argsort(-Vs)
    chosen = []
    # top-V record with the smallest box, then distinct classes (by V) with their smallest-box representative, largest V first
    seen = set()
    for i in order:
        key = (round(Vs[i], 4), R[i]['pairs_ext'])
        if key in seen: continue
        seen.add(key)
        cls = [j for j in range(len(R)) if (round(Vs[j], 4), R[j]['pairs_ext']) == key]
        j = min(cls, key=lambda j: sizes[j])
        chosen.append(j)
    limit = 25 if n <= 8 else 22
    chosen = chosen[:limit]
    extra = [int(j) for j in rng.permutation(len(R))[:max(0, 20 - len(chosen))]] if n <= 7 else []
    nchk = 0
    for j in chosen + extra:
        if time.time() - t0 > 200: break
        A = np.array(R[j]['A_ext'])
        res = box_bruteforce(A, 1 + 1e-6, maxsize=6e7)
        if res is None:
            out.append(dict(n=n, j=j, V=R[j]['V_ext'], skipped=True, size=sizes[j])); continue
        mn, cnt, size = res
        ok = (abs(mn - 1) < 1e-8) and (cnt == 2 * R[j]['pairs_ext'])
        out.append(dict(n=n, j=j, V=R[j]['V_ext'], min=mn, count=cnt, pairs_saved=R[j]['pairs_ext'], size=size, ok=bool(ok)))
        nchk += 1
    done = [o for o in out if o['n'] == n and not o.get('skipped')]
    print(f"n={n}: box-checked {len(done)} forms (distinct classes first, top V {max(o['V'] for o in done) if done else None}), "
          f"all ok: {all(o['ok'] for o in done)}, min range [{min(o['min'] for o in done):.10f},{max(o['min'] for o in done):.10f}], "
          f"skipped (box too big): {sum(1 for o in out if o['n']==n and o.get('skipped'))}, t={time.time()-t0:.0f}s", flush=True)
json.dump(out, open(os.path.join(here, 'box_check_results.json'), 'w'))
