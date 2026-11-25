import json, glob, os, csv, math, sys
import numpy as np
import eq
here = os.path.dirname(os.path.abspath(__file__))
best = {}
kiss = {}
with open(os.path.join(here, '..', 'R1_lit_perfect', 'best_known_lattices.csv')) as f:
    for row in csv.DictReader(f):
        n = int(row['n'])
        best[n] = float(row['V_eq_2^n_density'])
        kiss[n] = int(row['kissing_number']) if row['kissing_number'] else None
# more precise values of V_best from center density: V = delta * kappa_n * 2^n
with open(os.path.join(here, '..', 'R1_lit_perfect', 'best_known_lattices.csv')) as f:
    for row in csv.DictReader(f):
        n = int(row['n'])
        best[n] = float(row['center_density']) * eq.kappa(n) * 2 ** n
order = {'eis': 0, 'gauss': 1, 'hur': 2, 'zeta5': 3, 'zeta7': 4, 'zeta8': 5, 'zeta12': 6, 'zeta13': 7, 'ico': 8}
rows = []
bestforms = {}
maxratio = 0
for f in sorted(glob.glob(os.path.join(here, 'data', '*.json'))):
    r = json.load(open(f))
    recs = r['recs']
    if not recs:
        continue
    n, w, D = r['n'], r['w'], r['D']
    Vp = np.array([x['V_perfect'] for x in recs])
    Ve = np.array([x.get('V_ext', np.nan) for x in recs])
    op = np.array([x['orbits_perfect'] for x in recs])
    oe = np.array([x.get('orbits_ext', np.nan) for x in recs])
    Vb = best[n]
    okv = all(x['verify_perfect'][0] == 0 and x['verify_perfect'][1] == x['pairs_perfect'] and abs(x['verify_perfect'][2] - 1) < 1e-8
              for x in recs) and all(('verify_ext' not in x) or (x['verify_ext'][0] == 0 and x['verify_ext'][1] == x['pairs_ext']
                                                                  and abs(x['verify_ext'][2] - 1) < 1e-8) for x in recs)
    hitp = np.mean(np.abs(Vp / Vb - 1) < 1e-6)
    hite = np.nanmean(np.abs(Ve / Vb - 1) < 1e-6)
    Vall = np.nanmax(np.concatenate([Vp, Ve]))
    maxratio = max(maxratio, Vall / Vb)
    nd = sum(1 for x in recs if 'descent_error' in x)
    rows.append(dict(fam=r['name'], k=r['k'], n=n, w=w, D=D, mode=r['mode'], runs=len(recs), fails=r['fails'],
                     Vp_min=Vp.min(), Vp_med=np.median(Vp), Vp_mean=Vp.mean(), Vp_max=Vp.max(),
                     Ve_med=np.nanmedian(Ve), Ve_max=np.nanmax(Ve), orbP_mean=op.mean(), orbP_min=op.min(), orbP_max=op.max(),
                     simple_frac=float(np.mean(np.abs(op - D) < 1e-9)), orbE_mean=np.nanmean(oe),
                     ratio_med=np.median(Vp) / (w * D / 2), ratio_mean=Vp.mean() / (w * D / 2), ratioE_med=np.nanmedian(Ve) / (w * D / 2),
                     Vbest=Vb, hitp=hitp, hite=hite, best_ratio=Vall / Vb, verified=okv, desc_err=nd,
                     dsteps=np.nanmean([x.get('descent_steps', np.nan) for x in recs]),
                     tper=r['elapsed'] / len(recs)))
    # best form
    ib = int(np.nanargmax(np.where(np.isnan(Ve), Vp, Ve)))
    x = recs[ib]
    key = (r['name'], r['k'])
    Vx = x.get('V_ext', x['V_perfect'])
    if key not in bestforms or Vx > bestforms[key]['V'] * (1 + 1e-9):
        bestforms[key] = dict(V=Vx, A=x['A_final'], pairs=x.get('pairs_ext', x['pairs_perfect']), n=n, mode=r['mode'])
rows.sort(key=lambda d: (order.get(d['fam'], 99), d['k'], d['mode'] != 'random'))
out = []
out.append('| family | k | n | w | D | mode | runs | fails | V_perf min / med / mean / max | V_ext med / max | orbits at perfect mean [min,max] | frac exactly D | med V_perf/(wD/2) | med V_ext/(wD/2) | V_best known | hit freq perf / ext | best/V_best | LP steps | verified | s/run |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for d in rows:
    out.append(f"| {d['fam']} | {d['k']} | {d['n']} | {d['w']} | {d['D']} | {d['mode']} | {d['runs']} | {d['fails']} | "
               f"{d['Vp_min']:.2f} / {d['Vp_med']:.2f} / {d['Vp_mean']:.2f} / {d['Vp_max']:.2f} | {d['Ve_med']:.2f} / {d['Ve_max']:.2f} | "
               f"{d['orbP_mean']:.1f} [{d['orbP_min']:.0f},{d['orbP_max']:.0f}] | {d['simple_frac']:.2f} | {d['ratio_med']:.3f} | {d['ratioE_med']:.3f} | "
               f"{d['Vbest']:.2f} | {d['hitp']:.2f} / {d['hite']:.2f} | {d['best_ratio']:.4f} | {d['dsteps']:.2f} | {d['verified']} | {d['tper']:.2f} |")
out.append('')
out.append(f'max over all runs of V/V_best_known = {maxratio:.6f} (must be <= 1)')
out.append('')
# explicit forms
out.append('## Best form per (family,k): integer Gram data')
out.append('| family | k | n | V | V/V_best | kissing (2*pairs) | q (q*A integral, q = min of integral form) | det(q*A) |')
out.append('|---|---|---|---|---|---|---|---|')
expl = {}
for key in sorted(bestforms, key=lambda t: (order.get(t[0], 99), t[1])):
    b = bestforms[key]
    A = np.array(b['A'])
    q = eq.denominator(A)
    n = b['n']
    detq = np.linalg.det(q * A) if q else float('nan')
    out.append(f"| {key[0]} | {key[1]} | {n} | {b['V']:.4f} | {b['V']/best[n]:.4f} | {2*b['pairs']} | {q} | {detq:.4f} |")
    expl[f'{key[0]}_k{key[1]}'] = dict(V=b['V'], q=q, kissing=2 * b['pairs'], det_qA=detq,
                                         qA=np.rint(q * A).astype(int).tolist() if q else None, mode=b['mode'])
    # Hermitian matrix
    if key[0] in ('eis', 'gauss', 'hur') and q:
        gens0, w, info = eq.family(key[0])
        beta = np.array(info['beta'], dtype=float)
        d = len(beta)
        k = key[1]
        M = beta.copy()
        M[:, 1:] *= -1
        H = np.zeros((k, k, d))
        for i in range(k):
            for j in range(k):
                H[i, j] = np.linalg.solve(M, A[i * d, j * d:(j + 1) * d])
        # reconstruct check
        def mul(p, r):
            if d == 2:
                z = complex(*p) * complex(*r)
                return np.array([z.real, z.imag])
            return eq.qmul(p, r)
        def conj(p):
            c = np.array(p, dtype=float).copy(); c[1:] *= -1; return c
        A2 = np.zeros_like(A)
        for i in range(k):
            for j in range(k):
                for a in range(d):
                    for bb in range(d):
                        A2[i * d + a, j * d + bb] = mul(mul(conj(beta[a]), beta[bb]), H[i, j])[0]
        err = np.abs(A2 - A).max()
        expl[f'{key[0]}_k{key[1]}']['H_times_q'] = np.round(q * H, 6).tolist()
        expl[f'{key[0]}_k{key[1]}']['H_reconstruct_err'] = float(err)
json.dump(expl, open(os.path.join(here, 'best_forms.json'), 'w'), indent=0)
out.append('')
out.append('## Hermitian matrices q*H of the best forms (components w.r.t. 1,i resp. 1,i,j,k; form = Re sum conj(x_i) x_j h_ij type)')
for key, v in expl.items():
    if 'H_times_q' in v:
        out.append(f"* {key}: V={v['V']:.4f}, q={v['q']}, kissing={v['kissing']}, reconstruct err={v['H_reconstruct_err']:.1e}")
        H = np.array(v['H_times_q'])
        for i in range(H.shape[0]):
            out.append('    ' + '  '.join('(' + ','.join(f'{c:g}' for c in H[i, j]) + ')' for j in range(H.shape[1])))
# per-dimension comparison
out.append('')
out.append('## Per dimension: which symmetry is most efficient (random rule)')
out.append('| n | family (k) | w | D | wD/2 | med V_perf | med V_ext | max V | V_best | hit freq ext |')
out.append('|---|---|---|---|---|---|---|---|---|---|')
for d in sorted([d for d in rows if d['mode'] == 'random'], key=lambda d: (d['n'], -d['Vp_med'])):
    out.append(f"| {d['n']} | {d['fam']} ({d['k']}) | {d['w']} | {d['D']} | {d['w']*d['D']/2:.0f} | {d['Vp_med']:.2f} | {d['Ve_med']:.2f} | {max(d['Vp_max'], d['Ve_max']):.2f} | {d['Vbest']:.2f} | {d['hite']:.2f} |")
open(os.path.join(here, 'results_tables.md'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))

# ---- scaling fit: log(median V_perf) vs log(w D / 2), random rule, n >= 10
sel = [d for d in rows if d['mode'] == 'random' and d['n'] >= 10]
if len(sel) >= 3:
    x = np.log([d['w'] * d['D'] / 2 for d in sel])
    y = np.log([d['Vp_med'] for d in sel])
    p = np.polyfit(x, y, 1)
    res = y - np.polyval(p, x)
    ye = np.log([d['Ve_med'] for d in sel])
    pe = np.polyfit(x, ye, 1)
    y2 = np.log([d['Vp_med'] for d in sel]) - np.log([d['n'] * (d['n'] + 1) / 2 for d in sel])
    extra = ['', '## Scaling fit (random rule, n >= 10, %d (family,k) points)' % len(sel),
             f'log med V_perf = {p[0]:.3f} * log(wD/2) + {p[1]:.3f}  (rms residual {res.std():.3f}); '
             f'ratio med V_perf/(wD/2): min {np.exp((y-x).min()):.3f}, geometric mean {np.exp((y-x).mean()):.3f}, max {np.exp((y-x).max()):.3f}',
             f'log med V_ext  = {pe[0]:.3f} * log(wD/2) + {pe[1]:.3f}; ratio med V_ext/(wD/2): min {np.exp((ye-x).min()):.3f}, geometric mean {np.exp((ye-x).mean()):.3f}, max {np.exp((ye-x).max()):.3f}',
             f'for comparison med V_perf / N (N = n(n+1)/2, the unsymmetric count): min {np.exp(y2.min()):.3f}, geometric mean {np.exp(y2.mean()):.3f}, max {np.exp(y2.max()):.3f}']
    open(os.path.join(here, 'results_tables.md'), 'a').write('\n'.join(extra) + '\n')
    print('\n'.join(extra))
