"""Summary statistics recomputed from the verifier's own re-derived quantities (V from det of stored forms, pairs from own enumeration)."""
import json, numpy as np, math, sys
sys.path.insert(0, '.')
from vlib import Vbest, kap
R = json.load(open('recheck_main.json')) + json.load(open('recheck_scan.json'))
F2 = {(o['file'], o['idx'], o['tag']): o for o in json.load(open('recheck2_forms.json'))}
def ms(x):
    x = np.array(x, float)
    if len(x) == 0: return '-'
    if len(x) == 1: return '%.3f [1]' % x[0]
    return '%.3f+-%.3f' % (x.mean(), x.std(ddof=1) / math.sqrt(len(x)))
def val(r, tag, key):
    if tag == 'LP': return r['LP'][key]
    return F2[(r['file'], r['idx'], tag)][key]      # re-polished values for walk / greedy forms
pool = {}
for r in R: pool.setdefault((r['n'], r['kind']), []).append(r)
print('### pooled LP table (verifier)')
print('| n | family | #lat | Vbest/N | V_LP/N mean+-se | std | v_LP/N | pairs/N | simple | LP=best | B0sq mean | eig ratio mean | V0 mean |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
tot = 0
for (n, kind) in sorted(pool):
    P = pool[(n, kind)]; N = n * (n + 1) // 2; Vb = Vbest(n); tot += len(P)
    V = np.array([r['LP']['V'] for r in P]) / N
    print(f"| {n} | {kind} | {len(P)} | {Vb/N:.3f} | {ms(V)} | {V.std(ddof=1):.3f} | {np.mean([r['LP']['v'] for r in P])/N:.3f} | {np.mean([r['LP']['pairs'] for r in P])/N:.3f} | "
          f"{sum(r['LP']['pairs'] == N for r in P)}/{len(P)} | {sum(r['LP']['V'] >= Vb*(1-1e-6) for r in P)}/{len(P)} | {np.mean([r['LP']['B0sq'] for r in P]):.2f} | {np.mean([r['LP']['eig_ratio'] for r in P]):.2f} | {np.mean([r['V0'] for r in P]):.2f} |")
print('total LP vertices', tot)
print('\n### comparison table (main runs only; verifier)')
print('| n | family | #lat | V_LP/N | V_K/N [#] | V_Kd/N [#] | V_ext/N from LP [#] | V_ext/N from walk [#] | best: LP / K / Kd / extLP / extK | V_LP>V_K | V_LP=V_K | mean steps extLP | mean steps extK |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
main = {}
for r in R:
    if not r['file'].startswith('lpscan'): main.setdefault((r['n'], r['kind']), []).append(r)
for (n, kind) in sorted(main):
    P = main[(n, kind)]; N = n * (n + 1) // 2; Vb = Vbest(n)
    g = lambda tag: [val(r, tag, 'V') / N for r in P if tag in r]
    fr = lambda tag: '%d/%d' % (sum(val(r, tag, 'V') >= Vb * (1 - 1e-6) for r in P if tag in r), sum(tag in r for r in P))
    gt = sum(r['LP']['V'] > val(r, 'K', 'V') * (1 + 1e-9) for r in P if 'K' in r); eq = sum(abs(r['LP']['V'] - val(r, 'K', 'V')) <= 1e-9 * r['LP']['V'] for r in P if 'K' in r)
    nk = sum('K' in r for r in P)
    st = lambda tag: ('%.1f' % np.mean([val(r, tag, 'steps') for r in P if tag in r])) if any(tag in r for r in P) else '-'
    print(f"| {n} | {kind} | {len(P)} | {ms(g('LP'))} | {ms(g('K'))} [{len(g('K'))}] | {ms(g('Kd'))} [{len(g('Kd'))}] | {ms(g('extLP'))} [{len(g('extLP'))}] | {ms(g('extK'))} [{len(g('extK'))}] | "
          f"{fr('LP')} / {fr('K')} / {fr('Kd')} / {fr('extLP')} / {fr('extK')} | {gt}/{nk} | {eq}/{nk} | {st('extLP')} | {st('extK')} |")
# specific claims
print('\n### specific claims')
def cnt(n, pairs, kinds=None):
    P = [r for r in R if r['n'] == n and (kinds is None or r['kind'] in kinds)]
    return sum(r['LP']['pairs'] == pairs and r['LP']['V'] >= Vbest(n) * (1 - 1e-6) for r in P), len(P)
for n, p, nm in ((6, 36, 'E6'), (7, 63, 'E7'), (8, 120, 'E8'), (9, 136, 'Lambda9'), (10, 168, 'Lambda10'), (12, 378, 'K12')):
    print(nm, 'LP hits (pairs=%d & V=Vbest): %d/%d' % ((p,) + cnt(n, p)), ' gm: %d/%d' % cnt(n, p, ['gm']), ' gauss: %d/%d' % cnt(n, p, ['gauss']))
P = [r for r in R if 8 <= r['n'] <= 12]
print('simple vertices n=8..12: %d/%d' % (sum(r['LP']['pairs'] == r['n'] * (r['n'] + 1) // 2 for r in P), len(P)))
P = [r for r in R if r['n'] == 13]; print('simple n=13: %d/%d' % (sum(r['LP']['pairs'] == 91 for r in P), len(P)))
P = [r for r in R if r['n'] == 12 and r['kind'] == 'gm']
print('n=12 gm forms with >= 200 pairs:', sorted([(r['LP']['pairs'], round(r['LP']['V'] / 78, 3)) for r in P if r['LP']['pairs'] >= 200]))
P = [r for r in R if r['n'] == 12]
print('K12 reached by greedy from LP: gauss %d/%d, gm %d/%d' % (sum(val(r, 'extLP', 'V') >= Vbest(12) * (1 - 1e-6) for r in P if r['kind'] == 'gauss' and 'extLP' in r), sum(r['kind'] == 'gauss' and 'extLP' in r for r in P),
      sum(val(r, 'extLP', 'V') >= Vbest(12) * (1 - 1e-6) for r in P if r['kind'] == 'gm' and 'extLP' in r), sum(r['kind'] == 'gm' and 'extLP' in r for r in P)))
P = [r for r in R if r['n'] == 10]
print('Lambda10 reached: extLP %d/%d, extK %d/%d' % (sum(val(r, 'extLP', 'V') >= Vbest(10) * (1 - 1e-6) for r in P if 'extLP' in r), sum('extLP' in r for r in P), sum(val(r, 'extK', 'V') >= Vbest(10) * (1 - 1e-6) for r in P if 'extK' in r), sum('extK' in r for r in P)))
P = [r for r in R if r['n'] == 8]
print('E8: LP %d/%d ; sticky walk random sign %d/%d ; descent sign %d/%d' % (sum(r['LP']['V'] >= Vbest(8) * (1 - 1e-6) for r in P), len(P), sum(val(r, 'K', 'V') >= Vbest(8) * (1 - 1e-6) for r in P if 'K' in r), sum('K' in r for r in P), sum(val(r, 'Kd', 'V') >= Vbest(8) * (1 - 1e-6) for r in P if 'Kd' in r), sum('Kd' in r for r in P)))
for n in (8, 10, 12):
    P = [r for r in R if r['n'] == n and 'K' in r]
    print('n=%d: V_LP > V_K strictly %d/%d, equal %d/%d, V_LP < V_K %d/%d ; mean ratio V_LP/V_K = %.3f' % (n, sum(r['LP']['V'] > val(r, 'K', 'V') * (1 + 1e-9) for r in P), len(P), sum(abs(r['LP']['V'] - val(r, 'K', 'V')) <= 1e-9 * r['LP']['V'] for r in P), len(P),
          sum(r['LP']['V'] < val(r, 'K', 'V') * (1 - 1e-9) for r in P), len(P), np.mean([r['LP']['V'] for r in P]) / np.mean([val(r, 'K', 'V') for r in P])))
P = [r for r in R if r['n'] == 14 and 'extLP' in r]
print('n=14 greedy from LP:', [(round(val(r, 'extLP', 'V') / 105, 3), val(r, 'extLP', 'steps')) for r in P])
allV = [(r['n'], r[t]['V_over_Vbest']) for r in R for t in ('LP', 'K', 'Kd', 'extLP', 'extK') if t in r]
print('max V/Vbest over all %d stored forms: %.9f ; for n>=13: max V/N = %.3f' % (len(allV), max(v for n, v in allV), max(r[t]['V'] / (r['n'] * (r['n'] + 1) / 2) for r in R for t in ('LP', 'K', 'Kd', 'extLP', 'extK') if t in r and r['n'] >= 13)))
rd = [r['LP']['rounds'] for r in R]; nc = [r['LP']['ncons'] for r in R]
print('LP rounds: min %d max %d, distribution %s ; ncons min %d max %d ; #ncons>3500: %d' % (min(rd), max(rd), dict(zip(*np.unique(rd, return_counts=True))), min(nc), max(nc), sum(c > 3500 for c in nc)))
print('by n: max rounds / max ncons / #fallback', {n: (max(r['LP']['rounds'] for r in R if r['n'] == n), max(r['LP']['ncons'] for r in R if r['n'] == n), sum(r['LP']['fallback'] > 0 for r in R if r['n'] == n)) for n in sorted(set(r['n'] for r in R))})
