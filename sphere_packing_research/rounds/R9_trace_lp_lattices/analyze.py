#!/usr/bin/env python3
"""Tables for R9 (trace-LP on genuine random lattices)."""
import json, glob, os, math, sys
import numpy as np
from collections import Counter
here = os.path.dirname(os.path.abspath(__file__))
def ms(x):
    x = np.array([v for v in x if v is not None], float)
    if len(x) == 0: return '   -   '
    if len(x) == 1: return '%.3f (1)' % x[0]
    return '%.3f+-%.3f' % (x.mean(), x.std(ddof=1) / math.sqrt(len(x)))
def frac(recs, key, Vb):
    x = [r[key] for r in recs if key in r]
    if not x: return '-'
    return '%d/%d' % (sum(v >= Vb * (1 - 1e-6) for v in x), len(x))
rows = []; types = {}; summ = {}
files = sorted(glob.glob(os.path.join(here, 'data', 'n*_*.json')), key=lambda f: (int(os.path.basename(f).split('_')[0][1:]), f))
for f in files:
    d = json.load(open(f)); n = d['n']; N = d['N']; Vb = d['Vbest']; kind = d['kind']
    recs = [r for r in d['recs'] if 'V_LP' in r]
    nfail = len(d['recs']) - len(recs)
    g = lambda k, sc=1.0: [r[k] / sc for r in recs if k in r]
    s = dict(n=n, kind=kind, N=N, nlat=len(recs), lp_fail=nfail, Vbest_over_N=Vb / N,
             V0=ms(g('V0', N)), VLP=ms(g('V_LP', N)), vLP=ms(g('v_LP', N)), B0=ms(g('B0sq')), er=ms(g('eig_ratio')),
             pairs=ms(g('pairs_LP')), pairs_over_N=ms(g('pairs_LP', N)), simple='%d/%d' % (sum(r['pairs_LP'] == N for r in recs), len(recs)),
             VK=ms(g('V_K', N)), vK=ms(g('v_K', N)), VKd=ms(g('V_Kd', N)), VextLP=ms(g('V_extLP', N)), VextK=ms(g('V_extK', N)),
             pairsK=ms(g('pairs_K')), nK=len(g('V_K')), nKd=len(g('V_Kd')), nextLP=len(g('V_extLP')), nextK=len(g('V_extK')),
             fLP=frac(recs, 'V_LP', Vb), fK=frac(recs, 'V_K', Vb), fKd=frac(recs, 'V_Kd', Vb), fextLP=frac(recs, 'V_extLP', Vb),
             fextK=frac(recs, 'V_extK', Vb), rounds=ms([r['lp_info']['rounds'] for r in recs]),
             cact=ms([r['lp_info']['act_qlen_max'] for r in recs]), ncons=ms([r['lp_info']['ncons'] for r in recs]),
             stepsLP=ms(g('steps_extLP')), stepsK=ms(g('steps_extK')),
             maxV=max(max(g(k), default=0) for k in ('V_LP', 'V_K', 'V_Kd', 'V_extLP', 'V_extK')),
             minerr=max(max((abs(r[k] - 1) for r in recs if k in r), default=0) for k in ('min_LP', 'min_K', 'min_Kd', 'min_extLP', 'min_extK')),
             errs=sum(1 for r in d['recs'] for k in r if 'error' in k),
             LPgeK='%d/%d' % (sum(r['V_LP'] > r['V_K'] * (1 + 1e-9) for r in recs if 'V_K' in r), len(g('V_K'))),
             LPeqK='%d/%d' % (sum(abs(r['V_LP'] - r['V_K']) <= 1e-9 * r['V_K'] for r in recs if 'V_K' in r), len(g('V_K'))))
    summ[(n, kind)] = s
    for tag, kV, kp in (('LP', 'V_LP', 'pairs_LP'), ('K', 'V_K', 'pairs_K'), ('extLP', 'V_extLP', 'pairs_extLP'), ('extK', 'V_extK', 'pairs_extK')):
        types[(n, kind, tag)] = Counter((r[kp], round(r[kV] / N, 4)) for r in recs if kV in r)
out = []
P = out.append
P('# R9 tables: trace-LP vs sticky walk vs greedy descent on genuine random lattices (all numbers from tracelp.py runs)\n')
P('Normalisation: det Q = 1, V = kappa_n det(A)^(-1/2), N = n(n+1)/2. Entries are mean+-standard error over lattices. Poisson-model predictions: V_LP/N -> 2, v_LP/N -> 2/e = 0.736, V_K/N -> 1.\n')
P('## Table 1: trace-LP vertex\n')
P('| n | family | #lat | Vbest/N | V_0/N | V_LP/N | v_LP/N | B0sq | eig ratio | active pairs | pairs/N | simple vertex (pairs=N) | LP = best known | LP rounds | #constraints | max Q-length^2 of active vec / lambda_1^2 |')
P('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for (n, kind), s in summ.items():
    P(f"| {n} | {kind} | {s['nlat']} | {s['Vbest_over_N']:.3f} | {s['V0']} | {s['VLP']} | {s['vLP']} | {s['B0']} | {s['er']} | {s['pairs']} | {s['pairs_over_N']} | {s['simple']} | {s['fLP']} | {s['rounds']} | {s['ncons']} | {s['cact']} |")
P('\n## Table 2: comparison on the same lattices (V/N)\n')
P('| n | family | V_LP/N | V_K/N (sticky walk, random sign) [#] | V_Kd/N (sticky walk, descent sign) [#] | V_ext/N from LP vertex [#] | V_ext/N from walk end [#] | pairs at walk end | greedy steps from LP | greedy steps from walk |')
P('|---|---|---|---|---|---|---|---|---|---|')
for (n, kind), s in summ.items():
    P(f"| {n} | {kind} | {s['VLP']} | {s['VK']} [{s['nK']}] | {s['VKd']} [{s['nKd']}] | {s['VextLP']} [{s['nextLP']}] | {s['VextK']} [{s['nextK']}] | {s['pairsK']} | {s['stepsLP']} | {s['stepsK']} |")
P('\n## Table 3: fraction of cases reaching the best known lattice (V = Vbest up to 1e-6), and LP vs walk\n')
P('| n | family | trace-LP | sticky walk (random) | sticky walk (descent) | greedy from LP | greedy from walk | V_LP > V_K | V_LP = V_K | max V/Vbest over all forms | max abs(min A - 1) | #errors |')
P('|---|---|---|---|---|---|---|---|---|---|---|---|')
for (n, kind), s in summ.items():
    P(f"| {n} | {kind} | {s['fLP']} | {s['fK']} | {s['fKd']} | {s['fextLP']} | {s['fextK']} | {s['LPgeK']} | {s['LPeqK']} | {s['maxV'] / (s['Vbest_over_N'] * s['N']):.6f} | {s['minerr']:.1e} | {s['errs']} |")
P('\n## Table 4: which perfect forms are hit: (number of pairs of minimal vectors, V/N) : count\n')
for (n, kind, tag), c in types.items():
    if c:
        P(f"- n={n} {kind} {tag}: " + ', '.join(f"({p}, {v:.4f}): {k}" for (p, v), k in sorted(c.items(), key=lambda t: -t[1])[:8]) + (' ...' if len(c) > 8 else ''))

# ---- Table 5: LP-only scan (independent lattices), pooled with the stage-1 LPs of the main runs
pool = {}
for f in sorted(glob.glob(os.path.join(here, 'data', '*.json'))):
    b = os.path.basename(f)
    if not (b.startswith('lpscan_n') or (b.startswith('n') and b[1].isdigit())):
        continue
    d = json.load(open(f))
    for r in d['recs']:
        if 'V_LP' in r:
            pool.setdefault((d['n'], d['kind']), dict(N=d['N'], Vb=d['Vbest'], recs=[], fail=0))['recs'].append(r)
        else:
            pool.setdefault((d['n'], d['kind']), dict(N=d['N'], Vb=d['Vbest'], recs=[], fail=0))['fail'] += 1
P('\n## Table 5: trace-LP vertex, all lattices pooled (main runs + LP-only scan lpscan_*.json), per n and family\n')
P('| n | family | #lat (LP failures) | Vbest/N | V_LP/N mean+-se | sample std (Poisson/Gamma law: 2/sqrt N) | median | min | max | v_LP/N | pairs/N mean+-se | min pairs/N | simple (pairs=N) | B0sq | eig ratio | max abs(minA-1) | LP = best known | mean s per LP |')
P('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
scan = {}
for (n, kind) in sorted(pool):
    p = pool[(n, kind)]; N = p['N']; R = p['recs']; Vb = p['Vb']
    V = np.array([r['V_LP'] for r in R]) / N; pr = np.array([r['pairs_LP'] for r in R]) / N
    assert V.max() * N <= Vb * (1 + 1e-7)
    scan[f'{n}_{kind}'] = dict(n=n, kind=kind, nlat=len(R), fail=p['fail'], Vbest_over_N=Vb / N, VLP_mean=float(V.mean()),
                               VLP_se=float(V.std(ddof=1) / math.sqrt(len(V))) if len(V) > 1 else None, VLP_median=float(np.median(V)), VLP_std=float(V.std(ddof=1)) if len(V) > 1 else None,
                               VLP_min=float(V.min()), VLP_max=float(V.max()), pairs_over_N_mean=float(pr.mean()), pairs_over_N_min=float(pr.min()),
                               vLP_mean=float(np.mean([r['v_LP'] for r in R]) / N), nsimple=int(sum(r['pairs_LP'] == N for r in R)))
    P(f"| {n} | {kind} | {len(R)} ({p['fail']}) | {Vb/N:.3f} | {ms(V)} | {V.std(ddof=1):.3f} ({2/math.sqrt(N):.3f}) | {np.median(V):.3f} | {V.min():.3f} | {V.max():.3f} | {ms([r['v_LP']/N for r in R])} | {ms(pr)} | {pr.min():.3f} | "
      f"{sum(r['pairs_LP'] == N for r in R)}/{len(R)} | {ms([r['B0sq'] for r in R])} | {ms([r['eig_ratio'] for r in R])} | {max(abs(r['min_LP']-1) for r in R):.1e} | "
      f"{sum(r['V_LP'] >= Vb*(1-1e-6) for r in R)}/{len(R)} | {np.mean([r['t_lp'] for r in R]):.2f} |")
json.dump(scan, open(os.path.join(here, 'summary_lp_pooled.json'), 'w'), indent=1)
txt = '\n'.join(out) + '\n'
open(os.path.join(here, 'results_table.md'), 'w').write(txt)
json.dump({f'{n}_{k}': {a: b for a, b in s.items()} for (n, k), s in summ.items()}, open(os.path.join(here, 'summary.json'), 'w'), indent=1)
print(txt)
