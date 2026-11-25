"""Recompute the summary statistics from the raw json (independent of summary.py)."""
import json, glob, os, math, statistics as st
import numpy as np
from collections import Counter, defaultdict
here = os.path.dirname(os.path.abspath(__file__))
def kappa(n): return math.pi ** (n / 2) / math.gamma(n / 2 + 1)
cd = {2: 1 / (2 * 3 ** .5), 3: 1 / (4 * 2 ** .5), 4: 1 / 8, 5: 1 / (8 * 2 ** .5), 6: 1 / (8 * 3 ** .5), 7: 1 / 16, 8: 1 / 16,
      9: 1 / (16 * 2 ** .5), 10: 1 / (16 * 3 ** .5)}
Vb = {n: 2 ** n * kappa(n) * cd[n] for n in cd}
D = defaultdict(list); files = defaultdict(list)
for f in sorted(glob.glob(os.path.join(here, '..', 'data', 'runs_n*_*.json'))):
    d = json.load(open(f))
    D[(d['mode'], d['n'])].extend(d['recs']); files[(d['mode'], d['n'])].append((os.path.basename(f), len(d['recs'])))
rows = {}
print('mode n N runs | Vp mean/med/max | Ve mean/med/max | pairs_ext mean [min,max] | distinct perf/ext | freq best perf->ext | LP steps | nonedge | min eut | V/n^2 med perf, ext')
for (mode, n) in sorted(D):
    R = D[(mode, n)]
    Vp = [r['V_perfect'] for r in R]; Ve = [r['V_ext'] for r in R]
    pe = [r['pairs_ext'] for r in R]
    dp = len(set((round(r['V_perfect'], 5), r['pairs_perfect']) for r in R))
    de = len(set((round(r['V_ext'], 5), r['pairs_ext']) for r in R))
    fbp = sum(abs(v - Vb[n]) < 1e-5 * Vb[n] for v in Vp) / len(R)
    fbe = sum(abs(v - Vb[n]) < 1e-5 * Vb[n] for v in Ve) / len(R)
    steps = sum(r['descent_steps'] for r in R) / len(R)
    rows[(mode, n)] = dict(runs=len(R), Vp_med=st.median(Vp), Ve_med=st.median(Ve), Vp_mean=sum(Vp) / len(R), Ve_mean=sum(Ve) / len(R))
    print(f"{mode:8s} {n:2d} {n*(n+1)//2:2d} {len(R):4d} | {sum(Vp)/len(R):.3f} / {st.median(Vp):.3f} / {max(Vp):.3f} | "
          f"{sum(Ve)/len(R):.3f} / {st.median(Ve):.3f} / {max(Ve):.3f} | {sum(pe)/len(R):.1f} [{min(pe)},{max(pe)}] | {dp} / {de} | "
          f"{fbp:.3f} -> {fbe:.3f} | {steps:.2f} | {sum(r['nonedge'] for r in R)} | {min(r['eut_lp'] for r in R):.2e} | "
          f"{st.median(Vp)/n**2:.3f}, {st.median(Ve)/n**2:.3f} | Vp>Ve violations {sum(r['V_perfect']>r['V_ext']*(1+1e-9) for r in R)} | steps0 but Vp!=Ve {sum((r['descent_steps']==0) and abs(r['V_perfect']-r['V_ext'])>1e-7*r['V_ext'] for r in R)}")
print('\nTotal runs', sum(len(v) for v in D.values()))
print('files', dict(files))
# fits
print('\nFits (rms in log2 units)')
def fits(ns, y):
    ns = np.array(ns, float); y = np.array(y, float)
    pe = np.polyfit(ns, np.log2(y), 1); re = np.sqrt(np.mean((np.polyval(pe, ns) - np.log2(y)) ** 2))
    pp = np.polyfit(np.log2(ns), np.log2(y), 1); rp = np.sqrt(np.mean((np.polyval(pp, np.log2(ns)) - np.log2(y)) ** 2))
    return f"exp 2^({pe[0]:.3f} n + {pe[1]:.2f}) rms {re:.3f} | power {2**pp[1]:.3f} n^{pp[0]:.3f} rms {rp:.3f}"
for mode in ('descent', 'random'):
    for key, lab in (('Ve_med', 'extreme median'), ('Vp_med', 'perfect median')):
        for lo in (3, 5):
            ns = [n for n in range(lo, 11) if (mode, n) in rows]
            print(mode, lab, f'{lo}..10:', fits(ns, [rows[(mode, n)][key] for n in ns]))
for lo in (3, 5):
    ns = list(range(lo, 11)); print('best known', f'{lo}..10:', fits(ns, [Vb[n] for n in ns]))
# local slopes
print('\nlocal slopes of log2(median V_perfect) and log2(median V_ext):')
for mode in ('descent', 'random'):
    for a, b in ((7, 8), (8, 9), (9, 10)):
        if (mode, a) in rows and (mode, b) in rows:
            print(mode, f'{a}->{b}', round(math.log2(rows[(mode, b)]['Vp_med'] / rows[(mode, a)]['Vp_med']), 3),
                  round(math.log2(rows[(mode, b)]['Ve_med'] / rows[(mode, a)]['Ve_med']), 3), ' n^2 slope', round(2 * math.log2(b / a), 3))
# bootstrap error of medians and of slopes n=9,10
rng = np.random.default_rng(7)
print('\nbootstrap (2000 resamples) 16-84% interval of medians:')
bs = {}
for (mode, n) in sorted(D):
    if n < 8: continue
    for key in ('V_perfect', 'V_ext'):
        v = np.array([r[key] for r in D[(mode, n)]])
        m = np.median(v[rng.integers(0, len(v), (2000, len(v)))], axis=1)
        bs[(mode, n, key)] = m
        print(mode, n, key, 'median', round(float(np.median(v)), 3), 'CI68', np.round(np.percentile(m, [16, 84]), 2), 'CI95', np.round(np.percentile(m, [2.5, 97.5]), 2))
for mode in ('descent', 'random'):
    for a, b in ((8, 9), (9, 10)):
        s = np.log2(bs[(mode, b, 'V_perfect')] / bs[(mode, a, 'V_perfect')])
        print('slope', mode, f'{a}->{b}', 'V_perfect: sd', round(float(s.std()), 3), 'CI95', np.round(np.percentile(s, [2.5, 97.5]), 3))
# trajectory table
print('\nTrajectory: median V when contact rank first reaches fraction f of N')
for (mode, n) in sorted(D):
    if n < 6: continue
    N = n * (n + 1) // 2
    for lab, R in (('all', D[(mode, n)]),) + ((('first16', D[(mode, n)][:16]),) if (mode, n) == ('descent', 10) else ()):
        cols = []
        for f in (0.25, 0.5, 0.75, 0.9, 1.0):
            vals = []
            for r in R:
                v = [t[2] for t in r['traj'] if t[0] >= f * N]
                if v: vals.append(v[0])
            cols.append(st.median(vals) if vals else float('nan'))
        print(mode, n, lab, len(R), ' '.join(f'{c:.3f}' for c in cols))
# top forms per n
print('\nTop extreme forms (V, pairs, q): count')
for n in (8, 9, 10):
    C = Counter()
    for mode in ('descent', 'random'):
        for r in D.get((mode, n), []):
            C[(round(r['V_ext'], 3), r['pairs_ext'], r['denom_ext'], round(r['det_ext'] * r['denom_ext'] ** n))] += 1
    print(n, sorted(C.items(), reverse=True)[:8])
# pairs stats
for n in (9, 10):
    for mode in ('descent', 'random'):
        pe = [r['pairs_ext'] for r in D[(mode, n)]]
        print('pairs_ext', mode, n, 'mean', round(sum(pe) / len(pe), 1), 'min', min(pe), 'max', max(pe), ' /N:', round(sum(pe) / len(pe) / (n * (n + 1) / 2), 2))
