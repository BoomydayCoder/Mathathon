#!/usr/bin/env python3
"""Verifier's independent recomputation of the summary statistics from the raw JSONL records."""
import json, math, glob, collections, statistics as st
import numpy as np
D = '../data/'
def kappa(n): return math.pi ** (n / 2) / math.gamma(n / 2 + 1)
# best known lattices: (min norm m, Gram det) -> V = kappa_n m^(n/2)/sqrt(det)
best = {2: (2, 3), 3: (2, 4), 4: (2, 4), 5: (2, 4), 6: (2, 3), 7: (2, 2), 8: (2, 1), 9: (4, 512), 10: (4, 768)}  # A2,D3,D4,D5,E6,E7,E8,Lambda9,Lambda10
def Vb(n):
    m, dd = best[n]; return kappa(n) * m ** (n / 2) / math.sqrt(dd)
def zeta(n): return sum(k ** -float(n) for k in range(1, 100000))
rep = {  # reported: runs, meanV, medV, maxV, ElogV, fracbest, distinct, pairs_mean
 2: (207, 3.6276, 3.6276, 3.6276, 1.2886, 1.000, 1, 3.00), 3: (225, 5.9238, 5.9238, 5.9238, 1.7790, 1.000, 1, 6.00),
 4: (220, 9.3155, 8.8276, 9.8696, 2.2301, 0.468, 2, 10.94), 5: (290, 14.4829, 14.8882, 14.8882, 2.6706, 0.845, 3, 19.22),
 6: (223, 22.6742, 23.8686, 23.8686, 3.1178, 0.673, 5, 33.57), 7: (112, 35.1996, 37.7981, 37.7981, 3.5491, 0.786, 8, 56.79),
 8: (163, 44.7917, 36.3733, 64.9394, 3.7546, 0.331, 53, 71.36), 9: (100, 45.2659, 44.7534, 60.9406, 3.8081, 0.000, 94, 60.62),
 10: (50, 56.8051, 54.5642, 75.3837, 4.0333, 0.000, 49, 76.14)}
rep2 = {  # TK, mean tf, med tf, frac frozen by TK, mean V(TK), pairs at TK, mean V0, mean lam1
 2: (2.7726, 1.9333, 0.6801, 0.836, 3.2014, 2.69, 1.6185, 1.1846), 3: (1.9531, 0.8338, 0.6256, 0.933, 5.7351, 5.81, 1.6093, 1.0870),
 4: (1.3863, 0.8254, 0.6278, 0.873, 9.1638, 10.60, 1.7032, 1.0566), 5: (1.0300, 0.8448, 0.6751, 0.762, 14.1404, 18.16, 1.7337, 1.0436),
 6: (0.7963, 1.1253, 0.8645, 0.466, 20.6985, 27.83, 1.7043, 1.0299), 7: (0.6354, 1.2699, 0.8907, 0.295, 28.8431, 39.00, 1.7189, 1.0244),
 8: (0.5199, 1.0752, 0.6901, 0.362, 36.6234, 47.47, 1.9357, 1.0343), 9: (0.4340, 0.5545, 0.4885, 0.460, 43.0074, 54.17, 1.7759, 1.0136),
 10: (0.3684, 0.5000, 0.4217, 0.340, 53.2572, 63.80, 2.1974, 1.0396)}
repz = {3: (83, 5.9238, 5.9238, 5.9238, 1.000, 6.00, 0.4432), 4: (63, 9.3404, 8.8276, 9.8696, 0.492, 10.98, 0.5298),
        5: (55, 14.4909, 14.8882, 14.8882, 0.855, 19.27, 0.5844), 6: (46, 22.7240, 23.8686, 23.8686, 0.717, 33.72, 0.8030),
        8: (41, 45.2347, 37.7892, 64.9394, 0.341, 72.44, 0.6634)}
def cluster(vals, rel=1e-6):
    v = sorted(vals); groups = [[v[0]]]
    for x in v[1:]:
        if x - groups[-1][-1] <= rel * x: groups[-1].append(x)
        else: groups.append([x])
    return groups
res = {}; maxdev = 0.0
print('Vbest check:', {n: round(Vb(n), 4) for n in best})
ElogV = {}; anomalies = []
for kind in ['gm', 'zn']:
    for n in range(2, 11):
        try: R = [json.loads(l) for l in open(D + f'{kind}_n{n}.jsonl') if l.strip()]
        except FileNotFoundError: continue
        d = n * (n + 1) // 2
        nfroz = sum(1 for r in R if r['frozen'] and r['rank'] == d)
        V = np.array([r['Vfinal'] for r in R if r['frozen']])
        vb = Vb(n)
        over = [(i, r['Vfinal']) for i, r in enumerate(R) if r['Vfinal'] > vb * (1 + 1e-7)]
        fb = float(np.mean(np.abs(V - vb) < 1e-5 * vb)); nb = int(np.sum(np.abs(V - vb) < 1e-5 * vb))
        groups = cluster(V.tolist())
        npairs = np.array([r['final_npairs'] for r in R if r['frozen']])
        tf = np.array([r['tfreeze'] for r in R if r['frozen']]); TK = 16 * math.log(n) / n ** 2
        VT = np.array([r['VT_scaled'] for r in R]); pT = np.array([r['pairsT'] for r in R])
        lam = np.array([r['lam1'] for r in R]); V0 = np.array([r['V0'] for r in R])
        # consistency of V0 with lam1: V0 = lam1^n / 1.001^(n/2)
        v0dev = float(np.max(np.abs(V0 - lam ** n / 1.001 ** (n / 2)) / V0))
        # event sanity: rank increments, pairs >= rank, final event rank = d
        jumps = 0; nonmono = 0; lastsurplus = 0; evtot = 0; pairs_lt_rank = 0; multi_before_last = 0
        for r in R:
            ev = r['events']; pr = 0; pp = 0
            for j, (t, v, npr, rk) in enumerate(ev):
                evtot += 1
                if rk - pr != 1: jumps += 1
                if npr < rk: pairs_lt_rank += 1
                if (npr - pp) > 1 and j < len(ev) - 1: multi_before_last += 1
                pr = rk; pp = npr
            if ev and ev[-1][3] != d: nonmono += 1
        s = dict(kind=kind, n=n, runs=len(R), frozen=nfroz, meanV=float(V.mean()), medV=float(np.median(V)), maxV=float(V.max()), minV=float(V.min()),
                 ElogV=float(np.log(V).mean()), seLogV=float(np.log(V).std(ddof=1) / math.sqrt(len(V))), seV=float(V.std(ddof=1) / math.sqrt(len(V))),
                 twozeta=2 * zeta(n), V_n2=float(V.mean() / n ** 2), V_d=float(V.mean() / d), Vbest=vb, mean_best=float(V.mean() / vb),
                 fracbest=fb, nbest=nb, distinct=len(groups), pairs_mean=float(npairs.mean()), pairs_min=int(npairs.min()), pairs_max=int(npairs.max()),
                 frac_pairs_eq_d=float(np.mean(npairs == d)), TK=TK, tf_mean=float(tf.mean()), tf_med=float(np.median(tf)),
                 frac_frozen_TK=float(np.mean(tf <= TK)), VT_mean=float(VT.mean()), pT_mean=float(pT.mean()), V0_mean=float(V0.mean()), lam_mean=float(lam.mean()),
                 min_minq=float(min(r['final_minq'] for r in R)), max_minq=float(max(r['final_minq'] for r in R)),
                 tracked_ne_enum=int(sum(r['final_npairs'] != r['npairs_tracked'] for r in R)), warn_runs=int(sum(1 for r in R if r['warn'] > 0)),
                 over_best=over, v0dev=v0dev, rank_jumps_ne1=jumps, events=evtot, last_rank_ne_d=nonmono, pairs_lt_rank=pairs_lt_rank,
                 multi_pair_events_before_last=multi_before_last, n_lam_ge1=int((lam >= 1).sum()),
                 meanV_lam_ge1=float(V[lam >= 1].mean()) if (lam >= 1).any() else None, meanV_lam_lt1=float(V[lam < 1].mean()) if (lam < 1).any() else None,
                 top=[(round(g[0], 4), len(g)) for g in sorted(groups, key=lambda g: -len(g))[:8]])
        res[f'{kind}{n}'] = s
        if kind == 'gm':
            ElogV[n] = s['ElogV']
            mine = (s['runs'], s['meanV'], s['medV'], s['maxV'], s['ElogV'], s['fracbest'], s['distinct'], s['pairs_mean'])
            dev = max(abs(a - b) for a, b in zip(mine, rep[n])); maxdev = max(maxdev, dev)
            mine2 = (s['TK'], s['tf_mean'], s['tf_med'], s['frac_frozen_TK'], s['VT_mean'], s['pT_mean'], s['V0_mean'], s['lam_mean'])
            dev2 = max(abs(a - b) for a, b in zip(mine2, rep2[n]))
            print(f"gm n={n}: runs {s['runs']} frozen {s['frozen']} meanV {s['meanV']:.4f} (se {s['seV']:.3f}) med {s['medV']:.4f} max {s['maxV']:.4f} ElogV {s['ElogV']:.4f} (se {s['seLogV']:.4f}) "
                  f"V/n2 {s['V_n2']:.4f} V/d {s['V_d']:.4f} mean/best {s['mean_best']:.4f} fracbest {s['fracbest']:.3f} ({s['nbest']}) distinct {s['distinct']} pairs {s['pairs_mean']:.2f} [{s['pairs_min']},{s['pairs_max']}] frac=d {s['frac_pairs_eq_d']:.3f} | maxdev vs report tbl1 {dev:.2e}, tbl2 {dev2:.2e}")
            print(f"      TK {s['TK']:.4f} tf mean {s['tf_mean']:.4f} med {s['tf_med']:.4f} (med/TK {s['tf_med']/s['TK']:.3f}) frozen by TK {s['frac_frozen_TK']:.3f} V(TK) {s['VT_mean']:.4f} pairsT {s['pT_mean']:.2f} V0 {s['V0_mean']:.4f} lam1 {s['lam_mean']:.4f} "
                  f"| min minq {s['min_minq']:.10f} max minq {s['max_minq']:.10f} tracked!=enum {s['tracked_ne_enum']} warnruns {s['warn_runs']} overbest {len(over)} v0dev {v0dev:.1e} | events {evtot} rank-jumps!=1 {jumps} lastrank!=d {nonmono} multi-pair events before last {multi_before_last} | lam>=1: {s['n_lam_ge1']} meanV {s['meanV_lam_ge1']} ; lam<1 meanV {s['meanV_lam_lt1']}")
        else:
            mine = (s['runs'], s['meanV'], s['medV'], s['maxV'], s['fracbest'], s['pairs_mean'], s['tf_med'])
            dev = max(abs(a - b) for a, b in zip(mine, repz[n]))
            print(f"zn n={n}: runs {s['runs']} frozen {s['frozen']} meanV {s['meanV']:.4f} med {s['medV']:.4f} max {s['maxV']:.4f} fracbest {s['fracbest']:.3f} ({s['nbest']}) pairs {s['pairs_mean']:.2f} tf med {s['tf_med']:.4f} | maxdev vs report {dev:.2e} | min minq {s['min_minq']:.10f} overbest {len(over)} rank-jumps!=1 {jumps}/{evtot} tracked!=enum {s['tracked_ne_enum']}")
        print('      top values:', s['top'])
ns = np.array(sorted(ElogV)); y = np.array([ElogV[n] for n in ns]); yb = np.array([math.log(Vb(n)) for n in ns])
def fit(x, y):
    Am = np.vstack([np.ones_like(x), x]).T; c = np.linalg.lstsq(Am, y, rcond=None)[0]; return c, float(np.sqrt(np.mean((Am @ c - y) ** 2)))
for nm, yy in [('ElogV', y), ('log Vbest', yb)]:
    c1, r1 = fit(np.log(ns.astype(float)), yy); c2, r2 = fit(ns.astype(float), yy)
    print(f'fit {nm} n=2..10: a+b log n: a={c1[0]:.4f} b={c1[1]:.4f} rms={r1:.4f}; a+b n: a={c2[0]:.4f} b={c2[1]:.4f} rms={r2:.4f}')
    res['fit_' + nm] = dict(logfit=[float(c1[0]), float(c1[1]), r1], linfit=[float(c2[0]), float(c2[1]), r2])
c1, r1 = fit(np.log(ns[:5].astype(float)), y[:5]); c2, r2 = fit(ns[:5].astype(float), y[:5])
print(f'fit ElogV n=2..6: log rms {r1:.4f}, lin rms {r2:.4f}')
c1, r1 = fit(np.log(ns[2:].astype(float)), y[2:]); c2, r2 = fit(ns[2:].astype(float), y[2:])
print(f'fit ElogV n=4..10: a+b log n: b={c1[1]:.3f} rms {r1:.4f}; a + b n: b={c2[1]:.3f} rms {r2:.4f}')
# V at rank k / k table
print('V at first event with rank k, divided by k (mean over runs that have an event with exactly that rank; count):')
for n in range(4, 11):
    R = [json.loads(l) for l in open(D + f'gm_n{n}.jsonl') if l.strip()]; d = n * (n + 1) // 2
    byk = collections.defaultdict(list)
    for r in R:
        seen = set()
        for (t, v, npr, rk) in r['events']:
            if rk not in seen: seen.add(rk); byk[rk].append(v)
    ks = [k for k in (1, 2, 3, 5, 10, 15, 20, 28, 36, 45, d - 1, d) if k in byk]
    print(f' n={n}: ' + '; '.join(f'k={k}: {np.mean(byk[k])/k:.3f} ({len(byk[k])})' for k in ks))
    res[f'Vk{n}'] = {int(k): [float(np.mean(byk[k]) / k), len(byk[k])] for k in sorted(byk)}
json.dump(res, open('recomputed_stats.json', 'w'), indent=1, default=str)
print('max abs deviation from reported table 1 (gm):', maxdev)
