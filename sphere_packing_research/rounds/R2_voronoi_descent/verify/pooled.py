import json, glob, os, math
import numpy as np
def kappa(n): return math.pi ** (n / 2) / math.gamma(n / 2 + 1)
Vb = {9: 2 ** 9 * kappa(9) / (16 * 2 ** .5), 10: 2 ** 10 * kappa(10) / (16 * 3 ** .5), 8: 2 ** 8 * kappa(8) / 16}
orig = {}
for f in sorted(glob.glob('../data/runs_n*_*.json')):
    d = json.load(open(f)); orig.setdefault((d['n'], d['mode']), []).extend(d['recs'])
fresh = json.load(open('fresh_control.json'))
rng = np.random.default_rng(5)
def ci(v):
    v = np.array(v); m = np.median(v[rng.integers(0, len(v), (4000, len(v)))], axis=1)
    return np.round(np.percentile(m, [2.5, 97.5]), 2).tolist()
print('n mode | reported: runs, med Vp, med Ve, freq best (perf->ext) | fresh seeds: same | pooled: runs, med Vp [CI95], med Ve [CI95], freq best perf -> ext')
for n in (9, 10):
    for mode in ('descent', 'random'):
        O = orig[(n, mode)]; F = [o for o in fresh if o['n'] == n and o['mode'] == mode]
        vpo = [r['V_perfect'] for r in O]; veo = [r['V_ext'] for r in O]
        vpf = [o['Vp'] for o in F]; vef = [o['Ve'] for o in F]
        fb = lambda v: float(np.mean(np.abs(np.array(v) - Vb[n]) < 1e-5 * Vb[n])) if len(v) else float('nan')
        vp = vpo + vpf; ve = veo + vef
        print(f"{n} {mode:7s} | {len(O)}, {np.median(vpo):.2f}, {np.median(veo):.2f}, {fb(vpo):.3f}->{fb(veo):.3f} | "
              f"{len(F)}, {np.median(vpf) if F else float('nan'):.2f}, {np.median(vef) if F else float('nan'):.2f}, {fb(vpf):.3f}->{fb(vef):.3f} | "
              f"{len(vp)}, {np.median(vp):.2f} {ci(vp)}, {np.median(ve):.2f} {ci(ve)}, {fb(vp):.3f}->{fb(ve):.3f} (+-{math.sqrt(max(fb(ve)*(1-fb(ve)),1e-9)/len(ve)):.3f})"
              f" | Vp/n^2 {np.median(vp)/n**2:.3f} Ve/n^2 {np.median(ve)/n**2:.3f} | fresh mean LP steps {np.mean([o['steps'] for o in F]) if F else float('nan'):.2f} | fresh checks ok {all(o['okp'] and o['oke'] for o in F)}")
# pooled local slopes for perfect-stage medians (random rule), using reported n=8 median plus sensitivity control
s = json.load(open('sampler_sensitivity.json'))
print('sampler sensitivity table:')
for d in s:
    print(f"  n={d['n']} {d['mode']:7s} {d['variant']:5s} runs {d['runs']:3d} freq_best {d['freq_best']:.3f} +- {d['se']:.3f} median {d['median']:.3f} mean {d['mean']:.3f} distinct {d['distinct']} bad {d['bad_min']}")
