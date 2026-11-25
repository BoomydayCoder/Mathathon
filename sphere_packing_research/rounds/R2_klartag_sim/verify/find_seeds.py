#!/usr/bin/env python3
"""Identify the RNG seeds of the saved batches by matching lambda_1 of the first lattice of each batch."""
import sys, json, math
import numpy as np
sys.dont_write_bytecode = True
import sim_r2_copy as S
D = '../data/'
starts = {2: [0, 107], 3: [0], 4: [0, 200], 5: [0, 250], 6: [0, 198], 7: [0], 8: [0, 19, 89, 93], 9: [0], 10: [0]}
out = {}
for n in range(2, 11):
    R = [json.loads(l) for l in open(D + f'gm_n{n}.jsonl') if l.strip()]
    targets = {i: R[i]['lam1'] for i in starts[n] if i < len(R)}
    cands = sorted(set(list(range(0, 130)) + [k * 100 + n for k in range(0, 100)] + [k * 1000 + n for k in range(0, 10)] + [k*10+n for k in range(0,100)]))
    for s in cands:
        rng = np.random.default_rng(s)
        B = S.make_lattice(n, 'gm', rng)
        Bc = B.copy(); S.lll(Bc); G = Bc.T @ Bc
        Z = S.enum_points(G, G.diagonal().min())
        lam1 = math.sqrt(S.quad(G, Z).min())
        for i, t in targets.items():
            if abs(lam1 - t) < 1e-9 * t:
                print('n', n, 'batch start', i, 'seed', s, flush=True); out[f'{n}:{i}'] = s
json.dump(out, open('seeds_found.json', 'w'))
