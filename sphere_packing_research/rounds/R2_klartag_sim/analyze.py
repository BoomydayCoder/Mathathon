#!/usr/bin/env python3
"""Analyze the JSONL output of sim.py: statistics table, fits, contact-pairs-vs-volume table and plot."""
import json, glob, os, math, collections
import numpy as np
from math import gamma, pi
D = os.path.dirname(os.path.abspath(__file__))
def kappa(n): return pi ** (n / 2) / gamma(n / 2 + 1)
def zeta(n): return float(sum(1.0 / k ** n for k in range(1, 200000))) if n > 1 else float('inf')
# center densities of the best known (for n<=8 proven optimal) lattices
delta = {2: 1 / (2 * 3 ** .5), 3: 1 / (4 * 2 ** .5), 4: 1 / 8, 5: 1 / (8 * 2 ** .5), 6: 1 / (8 * 3 ** .5), 7: 1 / 16, 8: 1 / 16,
         9: 1 / (16 * 2 ** .5), 10: 1 / (16 * 3 ** .5), 11: 1 / (18 * 3 ** .5), 12: 1 / 27}
def Vbest(n): return 2 ** n * delta[n] * kappa(n)
rows = []; summ = {}; lines = []
for kind in ['gm', 'zn']:
    for f in sorted(glob.glob(os.path.join(D, 'data', kind + '_n*.jsonl')), key=lambda s: int(s.split('_n')[-1].split('.')[0])):
        R = [json.loads(l) for l in open(f) if l.strip()]
        if not R: continue
        nfail = sum(1 for l in open(f + '.failures')) if os.path.exists(f + '.failures') else 0
        n = R[0]['n']; d = n * (n + 1) // 2
        Fz = [r for r in R if r['frozen']]
        Vf = np.array([r['Vfinal'] for r in Fz]); vb = Vbest(n)
        best = np.mean(np.abs(Vf - vb) < 1e-4 * vb)
        cnt = collections.Counter([round(v, 3) for v in Vf])
        npairs = np.array([r['final_npairs'] for r in Fz])
        tf = np.array([r['tfreeze'] for r in Fz]); TK = R[0]['TK']
        VT = np.array([r['VT_scaled'] for r in R]); pT = np.array([r['pairsT'] for r in R])
        lam = np.array([r['lam1'] for r in Fz]); V0 = np.array([r['V0'] for r in Fz])
        big = lam >= 1.0
        s = dict(kind=kind, n=n, nfail=nfail, runs=len(R), frozen=len(Fz), warn=int(sum(r['warn'] for r in R)),
                 min_final_minq=float(min(r['final_minq'] for r in R)),
                 tracked_eq_enum=float(np.mean([r['final_npairs'] == r['npairs_tracked'] for r in Fz])),
                 V_mean=float(Vf.mean()), V_median=float(np.median(Vf)), V_max=float(Vf.max()), V_min=float(Vf.min()),
                 ElogV=float(np.log(Vf).mean()), se_logV=float(np.log(Vf).std() / math.sqrt(len(Vf))),
                 two_zeta=2 * zeta(n), n2=n * n, V_over_n2=float(Vf.mean() / n ** 2), Vbest=vb, frac_best=float(best),
                 mean_over_best=float(Vf.mean() / vb), n_distinct=len(cnt),
                 pairs_mean=float(npairs.mean()), pairs_min=int(npairs.min()), pairs_max=int(npairs.max()), d=d,
                 frac_pairs_eq_d=float(np.mean(npairs == d)),
                 tfreeze_mean=float(tf.mean()), tfreeze_median=float(np.median(tf)), TK=TK, frac_frozen_by_TK=float(np.mean(tf <= TK)),
                 VT_mean=float(VT.mean()), VT_median=float(np.median(VT)), pairsT_mean=float(pT.mean()),
                 V0_mean=float(V0.mean()), lam1_mean=float(lam.mean()),
                 V_mean_lam_ge1=float(Vf[big].mean()) if big.any() else None, n_lam_ge1=int(big.sum()),
                 top_values=cnt.most_common(8))
        # contact pairs vs V(t): mean V at the moment the rank reaches k
        byk = collections.defaultdict(list)
        for r in Fz:
            seen = set()
            for (t, V, npair, rk) in r['events']:
                if rk not in seen:
                    seen.add(rk); byk[rk].append((V, npair, t))
        s['V_at_rank'] = {int(k): [float(np.mean([x[0] for x in v])), float(np.mean([x[1] for x in v])), float(np.mean([x[2] for x in v])), len(v)] for k, v in sorted(byk.items())}
        summ[kind + str(n)] = s; rows.append(s)
json.dump(summ, open(os.path.join(D, 'summary.json'), 'w'), indent=1)
L = []
L.append('| lattice | n | runs (frozen) | mean V | median V | max V | min V | E[log V] | 2zeta(n) | mean V/n^2 | best known V_n | mean/best | frac at densest | #distinct final V | pairs mean [min,max] (d) | frac pairs=d |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for s in rows:
    L.append(f"| {s['kind']} | {s['n']} | {s['runs']} ({s['frozen']}) | {s['V_mean']:.4f} | {s['V_median']:.4f} | {s['V_max']:.4f} | {s['V_min']:.4f} | {s['ElogV']:.4f}+-{s['se_logV']:.4f} | {s['two_zeta']:.4f} | {s['V_over_n2']:.4f} | {s['Vbest']:.4f} | {s['mean_over_best']:.4f} | {s['frac_best']:.3f} | {s['n_distinct']} | {s['pairs_mean']:.2f} [{s['pairs_min']},{s['pairs_max']}] ({s['d']}) | {s['frac_pairs_eq_d']:.3f} |")
L.append('')
L.append('| lattice | n | T_K=16 log n/n^2 | mean freeze time | median freeze time | frac frozen by T_K | mean V(T_K) | mean pairs at T_K | mean V(0)=lam1^n/1.001^(n/2) | mean lam1 | mean V_final given lam1>=1 (count) | V_final/d | failed runs | warnings | min final min A x.x | tracked==enumerated contacts |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for s in rows:
    cond = f"{s['V_mean_lam_ge1']:.4f} ({s['n_lam_ge1']})" if s['V_mean_lam_ge1'] is not None else 'n/a'
    L.append(f"| {s['kind']} | {s['n']} | {s['TK']:.4f} | {s['tfreeze_mean']:.4f} | {s['tfreeze_median']:.4f} | {s['frac_frozen_by_TK']:.3f} | {s['VT_mean']:.4f} | {s['pairsT_mean']:.2f} | {s['V0_mean']:.4f} | {s['lam1_mean']:.4f} | {cond} | {s['V_mean']/s['d']:.4f} | {s['nfail']} | {s['warn']} | {s['min_final_minq']:.12f} | {s['tracked_eq_enum']:.3f} |")
L.append('')
# fits
g = [s for s in rows if s['kind'] == 'gm' and s['frozen'] >= 5]
ns = np.array([s['n'] for s in g], float); y = np.array([s['ElogV'] for s in g]); yb = np.array([math.log(s['Vbest']) for s in g])
def fit(x, y):
    Am = np.vstack([np.ones_like(x), x]).T; c, res, *_ = np.linalg.lstsq(Am, y, rcond=None)
    return c, float(np.sqrt(np.mean((Am @ c - y) ** 2)))
if len(g) >= 3:
    for name, yy in [('E[log V_final] (Klartag process, random lattice)', y), ('log(best known V_n)', yb)]:
        c1, r1 = fit(np.log(ns), yy); c2, r2 = fit(ns, yy)
        L.append(f"Fit of {name} over n={int(ns.min())}..{int(ns.max())}: a+b*log n: a={c1[0]:.4f}, b={c1[1]:.4f}, rms={r1:.4f};  a+b*n: a={c2[0]:.4f}, b={c2[1]:.4f}, rms={r2:.4f}")
    summ['fits_n'] = ns.tolist()
L.append('')
L.append('Most frequent final V values (value: count) per n:')
for s in rows:
    L.append(f"- {s['kind']} n={s['n']}: " + ', '.join(f"{v}: {c}" for v, c in s['top_values']))
L.append('')
L.append('Contact pairs versus volume: mean over runs of V(t) at the moment the rank of the contact set reaches k (k, mean V, mean #pairs, mean t, V/k):')
for s in rows:
    if s['kind'] != 'gm': continue
    ks = sorted(s['V_at_rank'].keys()); d = s['d']
    sel = sorted(set([k for k in ks if k in (1, 2, 3, 5, 10, 15, 20, 28, 36, 45, 55) or k == d or k == d - 1 or k == d // 2]))
    L.append(f"- n={s['n']} (d={d}): " + '; '.join(f"k={k}: V={s['V_at_rank'][k][0]:.3f}, pairs={s['V_at_rank'][k][1]:.2f}, t={s['V_at_rank'][k][2]:.4f}, V/k={s['V_at_rank'][k][0]/k:.3f}" for k in sel))
open(os.path.join(D, 'results_table.md'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for s in rows:
        if s['kind'] != 'gm': continue
        ks = sorted(s['V_at_rank'].keys())
        ax[0].plot([s['V_at_rank'][k][0] for k in ks], ks, marker='.', label=f"n={s['n']}")
    ax[0].set_xlabel('mean V(t) when contact rank reaches k'); ax[0].set_ylabel('k (independent contact pairs)'); ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].legend(fontsize=7); ax[0].set_title('contact pairs vs swept volume')
    ax[1].plot(ns, np.exp(y), 'o-', label='exp E[log V_final], Klartag process'); ax[1].plot(ns, np.exp(yb), 's-', label='best known V_n'); ax[1].plot(ns, ns ** 2, ':', label='n^2'); ax[1].plot(ns, [2 * zeta(int(k)) for k in ns], '--', label='2 zeta(n)')
    ax[1].set_yscale('log'); ax[1].set_xlabel('n'); ax[1].legend(fontsize=7); ax[1].set_title('final V versus n')
    fig.tight_layout(); fig.savefig(os.path.join(D, 'klartag_sim.png'), dpi=110)
except Exception as e:
    print('plot failed', e)
