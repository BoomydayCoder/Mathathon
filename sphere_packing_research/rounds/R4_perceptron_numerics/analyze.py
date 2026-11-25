#!/usr/bin/env python3
"""Tables (mean +- standard error, number of samples) from data/*.jsonl -> results_table.md"""
import json, glob, math, os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
out = []


def ms(x):
    x = np.asarray([v for v in x if v is not None], float)
    if len(x) == 0: return 'n/a'
    if len(x) == 1: return '%.4g (1 sample)' % x[0]
    return '%.4g +- %.2g' % (x.mean(), x.std(ddof=1) / math.sqrt(len(x)))


def load(pat):
    rows = []
    for f in sorted(glob.glob(os.path.join(here, 'data', pat))):
        rows += [json.loads(l) for l in open(f)]
    return rows


A = load('A_N*.jsonl')
out.append('## Part A: LP (maximise t), abstract Gaussian/Poisson model.  Predictions: v*/N = 2/e = 0.7358, rho^2 = 8, active = N+1, V/N = 2, active shifted offsets: mean 2, sd sqrt(12) = 3.46 (see notes D2)\n')
out.append('| N | samples (failed) | v*/N | rho^2=|theta*|^2 | active | v* e^{rho^2/8}/N | mean act. sigma | sd act. sigma | max act. sigma / min over samples of (cutoff - max act. sigma) | #near-cutoff constraints with slack<1 |')
out.append('|---|---|---|---|---|---|---|---|---|---|')
for N in sorted(set(r['N'] for r in A)):
    L = [r for r in A if r['N'] == N and r['part'] == 'A_LP']
    F = [r for r in A if r['N'] == N and r['part'] == 'A_LP_FAILED']
    if not L: continue
    out.append('| %d | %d (%d) | %s | %s | %s | %s | %s | %s | %.2f / %.2f | %d |' % (
        N, len(L), len(F), ms([r['vstar_over_N'] for r in L]), ms([r['rho2'] for r in L]), ms([r['nactive'] for r in L]),
        ms([r['V_over_N'] for r in L]), ms([r['act_sig_mean'] for r in L]), ms([r['act_sig_std'] for r in L]),
        max(r['act_sig_max'] for r in L), min(r['cutoff_shifted'] - r['act_sig_max'] for r in L), sum(r['n_nearcut_slack_lt1'] for r in L)))
L = [r for r in A if r['part'] == 'A_LP']
if L:
    out.append('| all | %d | %s | %s | (N+1 in %d of %d samples) | %s | %s | %s | | |' % (len(L), ms([r['vstar_over_N'] for r in L]), ms([r['rho2'] for r in L]),
        sum(1 for r in L if r['nactive'] == r['N'] + 1), len(L), ms([r['V_over_N'] for r in L]), ms([r['act_sig_mean'] for r in L]), ms([r['act_sig_std'] for r in L])))
tp = os.path.join(here, 'trunc_pred.txt')
if os.path.exists(tp):
    out.append('\nGordon prediction with the actual cutoff (quadrature): ' + ' ; '.join(
        'N=%d: v*/N=%.5f' % (d['N'], d['vstar_over_N_trunc']) for d in (json.loads(l) for l in open(tp) if l.startswith('{'))))
out.append('\n## Part A: QP (minimum norm at fixed v).  Predictions: rho^2 = root of x e^{-x/8}/4 = v/N, c/N = rho^2/8, V/N = rho^2/4\n')
out.append('| N | v/N | samples | rho^2 measured | rho^2 pred | c/N measured | c/N pred | V/N = (v/N) e^{rho^2/8} measured | V/N pred | V/(2c) measured |')
out.append('|---|---|---|---|---|---|---|---|---|---|')
for N in sorted(set(r['N'] for r in A)):
    for vf in sorted(set(r['v_over_N'] for r in A if r['part'] == 'A_QP')):
        L = [r for r in A if r['N'] == N and r['part'] == 'A_QP' and r['v_over_N'] == vf and not r.get('failed')]
        if not L: continue
        out.append('| %d | %.2f | %d | %s | %.4f | %s | %.4f | %s | %.4f | %s |' % (
            N, vf, len(L), ms([r['rho2'] for r in L]), L[0]['rho2_pred'], ms([r['c_over_N'] for r in L]), L[0]['c_over_N_pred'],
            ms([r['V_over_N'] for r in L]), L[0]['V_over_N_pred'], ms([r['V_over_N'] / (2 * r['c_over_N']) for r in L if r['c_over_N'] > 0])))

B = load('B_LP_n*.jsonl')
out.append('\n## Part B(i): trace-LP in the geometric Poisson model.  Predictions (n -> infinity): v_LP/N = 0.7358, V_LP/N = 2, |B0|^2 = 4, active = N\n')
out.append('| n | N | samples | all PD | V_LP/N | v_LP/N | |B0|_HS^2 | log(V/v) | active (min..max) | max active point volume/N (cutoff) | max E[#unsampled pairs inside] | lam_max/lam_min | mean, sd of log(vol_active/v_LP) (pred 1, 1.73) |')
out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for n in sorted(set(r['n'] for r in B)):
    L = [r for r in B if r['n'] == n and r['part'] == 'B_LP']
    P = [r for r in L if r['pd']]
    out.append('| %d | %d | %d | %s | %s | %s | %s | %s | %d..%d | %.0f (%.0f) | %.3g | %s | %s ; %s |' % (
        n, L[0]['N'], len(L), all(r['pd'] for r in L), ms([r['V_over_N'] for r in P]), ms([r['v_over_N'] for r in P]),
        ms([r['B0sq'] for r in P]), ms([math.log(r['V_over_N'] / r['v_over_N']) for r in P]),
        min(r['nactive'] for r in L), max(r['nactive'] for r in L),
        max(r['act_vol_max_over_N'] for r in L), L[0]['Vmax_over_N'], max(r['missed_mean'] for r in P),
        ms([r['lam_max'] / r['lam_min'] for r in P]), ms([r['act_logvolratio_mean'] for r in P]), ms([r['act_logvolratio_std'] for r in P])))

E = load('B_EXT_n*.jsonl')
if E:
    out.append('\n## Part B(iii): pivoting ascent of the volume from the LP vertex (same point sets as the B_LP rows of these files)\n')
    out.append('| n | samples | status counts | pivots | V_LP/N (start) | V/N at stop | |B0|^2 at stop | lam_max/lam_min at stop | E[#unsampled pairs inside] at stop |')
    out.append('|---|---|---|---|---|---|---|---|---|')
    for n in sorted(set(r['n'] for r in E)):
        L = [r for r in E if r['n'] == n and r['part'] == 'B_EXT' and 'V_over_N' in r]
        S = [r for r in E if r['n'] == n and r['part'] == 'B_LP']
        if not L: continue
        st = {}
        for r in L: st[r['status']] = st.get(r['status'], 0) + 1
        out.append('| %d | %d | %s | %s | %s | %s | %s | %s | %s |' % (
            n, len(L), st, ms([r['npiv'] for r in L]), ms([r['V_over_N'] for r in S]), ms([r['V_over_N'] for r in L]),
            ms([r['B0sq'] for r in L]), ms([r['lam_max'] / r['lam_min'] for r in L]), ms([r['missed_mean'] for r in L])))
    out.append('\nTrajectories (pivot number, V/N, |B0|^2, circumscribed ball volume/N, lam_max/lam_min), first sample of each n, subsampled:')
    for n in sorted(set(r['n'] for r in E)):
        L = [r for r in E if r['n'] == n and r['part'] == 'B_EXT' and r.get('traj')]
        if L:
            tr = L[0]['traj']; step = max(1, len(tr) // 8)
            out.append('- n=%d: ' % n + '; '.join('(%d, %.3f, %.2f, %.3g, %.1f)' % tuple(x) for x in tr[::step]))

K = load('K_n*.jsonl')
if K:
    out.append('\n## Part B(ii): Klartag sticky process in the infinite Poisson model (exact lazy revelation, partK.py), start = largest empty ball, dt = 1.25e-3/n (units r_1 = 1)\n')
    out.append('| n | N | runs | frozen | capped/stuck | V_K/N (frozen) | median V_K/N | contact pairs at freeze | freezing time t | lam_max/lam_min (median) | |B0|^2 | revealed points | min q over revealed points |')
    out.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for n in sorted(set(r['n'] for r in K)):
        L = [r for r in K if r['n'] == n]
        Fz = [r for r in L if r['status'] == 'frozen']
        if not Fz:
            out.append('| %d | %d | %d | 0 | %d | n/a | | | | | | | |' % (n, L[0]['N'], len(L), len(L)))
            continue
        out.append('| %d | %d | %d | %d | %d | %s | %.3f | %s | %s | %.1f | %s | %s | %.12f |' % (
            n, L[0]['N'], len(L), len(Fz), len(L) - len(Fz), ms([r['V_over_N'] for r in Fz]), float(np.median([r['V_over_N'] for r in Fz])),
            ms([r['npairs'] for r in Fz]), ms([r['t'] for r in Fz]), float(np.median([r['lam_ratio'] for r in Fz])),
            ms([r['B0sq'] for r in Fz]), ms([r['n_points'] for r in Fz]), min(r['min_q'] for r in Fz)))
    out.append('\nV(t)/k along the path (k = rank of the contact set), mean over frozen runs of V/k at k = N/4, N/2, 3N/4, N:')
    for n in sorted(set(r['n'] for r in K)):
        Fz = [r for r in K if r['n'] == n and r['status'] == 'frozen']
        if not Fz: continue
        N = Fz[0]['N']; res = []
        for frac in (0.25, 0.5, 0.75, 1.0):
            vals = []
            for r in Fz:
                tr = r['traj']; kk = min(tr, key=lambda x: abs(x[0] - frac * N))
                vals.append(kk[2] * N / kk[0])
            res.append('%.2f: %s' % (frac, ms(vals)))
        out.append('- n=%d: ' % n + ' ; '.join(res))
Kd = load('Kdt_n*.jsonl')
if Kd:
    Fz = [r for r in Kd if r['status'] == 'frozen']
    out.append('\nTime-step test (n=%d, dt/4): %d runs, %d frozen, V_K/N = %s' % (Kd[0]['n'], len(Kd), len(Fz), ms([r['V_over_N'] for r in Fz])))

txt = '\n'.join(out)
open(os.path.join(here, 'results_table.md'), 'w').write(txt + '\n')
print(txt)
