import json, glob, os, math, sys
import numpy as np
from collections import Counter
from mpmath import zeta

here = os.path.dirname(os.path.abspath(__file__))


def kappa(n):
    return math.pi ** (n / 2) / math.gamma(n / 2 + 1)


# best known lattice center densities (Conway-Sloane): V = 2^n kappa_n delta
delta = {2: 1 / (2 * math.sqrt(3)), 3: 1 / (4 * math.sqrt(2)), 4: 1 / 8, 5: 1 / (8 * math.sqrt(2)),
         6: 1 / (8 * math.sqrt(3)), 7: 1 / 16, 8: 1 / 16, 9: 1 / (16 * math.sqrt(2)), 10: 1 / (16 * math.sqrt(3))}
Vbest = {n: 2 ** n * kappa(n) * d for n, d in delta.items()}

data = {}
for f in sorted(glob.glob(os.path.join(here, 'data', 'runs_n*_*.json'))):
    d = json.load(open(f))
    data.setdefault((d['mode'], d['n']), []).extend(d['recs'])

lines = []
rows = []
for mode in sorted(set(k[0] for k in data)):
    lines.append(f'\n### mode = {mode} (sign rule of the sticky walk)\n')
    lines.append('| n | N=n(n+1)/2 | #samples | V start (med) | V perfect mean / med / max | V extreme mean / med / max | pairs perf (mean) | pairs ext mean [min,max] | #distinct perfect | #distinct extreme | freq best known | descent steps (mean) | 2 zeta(n) | n^2 | best known V_n |')
    lines.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for (md, n) in sorted(data):
        if md != mode:
            continue
        R = data[(md, n)]
        Vs = np.array([r['V_start'] for r in R])
        Vp = np.array([r['V_perfect'] for r in R])
        Ve = np.array([r['V_ext'] for r in R])
        pp = np.array([r['pairs_perfect'] for r in R])
        pe = np.array([r['pairs_ext'] for r in R])
        dp = Counter((round(r['det_perfect'] / kappa(n) ** 0, 7 if n < 9 else 8), r['pairs_perfect']) for r in R)
        de = Counter((round(r['V_ext'], 5), r['pairs_ext']) for r in R)
        dpp = Counter((round(r['V_perfect'], 5), r['pairs_perfect']) for r in R)
        fb = np.mean(np.abs(Ve - Vbest[n]) < 1e-5 * Vbest[n])
        fbp = np.mean(np.abs(Vp - Vbest[n]) < 1e-5 * Vbest[n])
        ds = np.mean([r['descent_steps'] for r in R])
        lines.append(f'| {n} | {n*(n+1)//2} | {len(R)} | {np.median(Vs):.3f} | {Vp.mean():.3f} / {np.median(Vp):.3f} / {Vp.max():.3f} | '
                     f'{Ve.mean():.3f} / {np.median(Ve):.3f} / {Ve.max():.3f} | {pp.mean():.1f} | {pe.mean():.1f} [{pe.min()},{pe.max()}] | '
                     f'{len(dpp)} | {len(de)} | {fb:.3f} (perfect stage: {fbp:.3f}) | {ds:.2f} | {2*float(zeta(n)):.3f} | {n*n} | {Vbest[n]:.3f} |')
        rows.append(dict(mode=mode, n=n, N=n * (n + 1) // 2, samples=len(R), Vstart_med=float(np.median(Vs)),
                         Vp_mean=float(Vp.mean()), Vp_med=float(np.median(Vp)), Vp_max=float(Vp.max()),
                         Ve_mean=float(Ve.mean()), Ve_med=float(np.median(Ve)), Ve_max=float(Ve.max()),
                         Ve_min=float(Ve.min()), pairs_perf=float(pp.mean()), pairs_ext=float(pe.mean()),
                         distinct_perfect=len(dpp), distinct_extreme=len(de), freq_best=float(fb),
                         freq_best_perfect=float(fbp), descent_steps=float(ds), Vbest=Vbest[n],
                         eut_min=float(np.nanmin([r['eut_lp'] for r in R])),
                         nonedge=int(sum(r['nonedge'] for r in R)),
                         denoms=Counter(r['denom_ext'] for r in R).most_common(8),
                         top_forms=[(k, v) for k, v in de.most_common(8)],
                         time_per_sample=float(np.mean([r['time'] for r in R]))))

# fits
lines.append('\n### Growth fits (median V of extreme forms and of first perfect forms)\n')
for mode in sorted(set(k[0] for k in data)):
    rr = [r for r in rows if r['mode'] == mode and r['n'] >= 3]
    if len(rr) < 4:
        continue
    ns = np.array([r['n'] for r in rr], float)
    for key, lab in [('Ve_med', 'extreme median'), ('Vp_med', 'perfect median'), ('Ve_max', 'extreme max'), ('Vbest', 'best known')]:
        y = np.array([r[key] for r in rr])
        for lo in (3, 5):
            m = ns >= lo
            if m.sum() < 3:
                continue
            pe_ = np.polyfit(ns[m], np.log2(y[m]), 1)
            re_ = np.sqrt(np.mean((np.polyval(pe_, ns[m]) - np.log2(y[m])) ** 2))
            pp_ = np.polyfit(np.log(ns[m]), np.log(y[m]), 1)
            rp_ = np.sqrt(np.mean((np.polyval(pp_, np.log(ns[m])) - np.log(y[m])) ** 2)) / math.log(2)
            lines.append(f'- mode={mode}, {lab}, n in [{lo},{int(ns.max())}]: exponential fit V ~ 2^({pe_[0]:.3f} n + {pe_[1]:.2f}), rms(log2)={re_:.3f}; '
                         f'power fit V ~ {math.exp(pp_[1]):.3f} n^{pp_[0]:.3f}, rms(log2)={rp_:.3f}')

lines.append('\n### Details per n\n')
for r in rows:
    lines.append(f"- mode={r['mode']} n={r['n']}: min eutaxy coeff over all final forms = {r['eut_min']:.3e}; non-edge LP directions = {r['nonedge']}; "
                 f"denominators q (q*A integral) top counts = {r['denoms']}; top extreme forms (V, pairs):count = {r['top_forms']}; "
                 f"V_ext min = {r['Ve_min']:.4f}; time/sample = {r['time_per_sample']:.2f}s")

# trajectory: mean V as function of rank during the sticky walk
lines.append('\n### Sticky walk trajectory: median V when rank of contact set first reaches fraction f of N\n')
lines.append('| mode | n | f=0 (start) | 0.25 | 0.5 | 0.75 | 0.9 | 1.0 (perfect) |')
lines.append('|---|---|---|---|---|---|---|---|')
for (md, n) in sorted(data):
    R = data[(md, n)]
    N = n * (n + 1) // 2
    cols = []
    for f in (0.25, 0.5, 0.75, 0.9, 1.0):
        vals = []
        for r in R:
            v = [t[2] for t in r['traj'] if t[0] >= f * N]
            if v:
                vals.append(v[0])
        cols.append(np.median(vals) if vals else float('nan'))
    lines.append(f"| {md} | {n} | {np.median([r['V_start'] for r in R]):.3f} | " + ' | '.join(f'{c:.3f}' for c in cols) + ' |')

open(os.path.join(here, 'results_tables.md'), 'w').write('\n'.join(lines) + '\n')
json.dump(rows, open(os.path.join(here, 'data', 'summary_rows.json'), 'w'), indent=1, default=str)
print('\n'.join(lines))
