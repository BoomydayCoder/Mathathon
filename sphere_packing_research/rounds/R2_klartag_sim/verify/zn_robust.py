#!/usr/bin/env python3
"""Robustness check with R2's own simulator: start from scaled Z^n at n=9,10 (R2 only ran Z^n up to n=8), to test whether the
n=9,10 statistics depend on the coarse Goldstein-Mayer ensemble (p=100003).  Final forms re-verified with the verifier's enumeration."""
import sys, json, math, time
import numpy as np
sys.dont_write_bytecode = True
import sim_r2_copy as S
import myenum as E
n = int(sys.argv[1]); nruns = int(sys.argv[2]); seed = int(sys.argv[3]); cap = float(sys.argv[4])
rng = np.random.default_rng(seed); out = open(f'zn_robust_n{n}.jsonl', 'w'); t0 = time.time(); k = 0
while k < nruns and time.time() - t0 < cap:
    B = S.make_lattice(n, 'zn', rng); sm = S.Sim(n, B, rng, 1.25e-3 / n); steps = 0
    while sm.Qc.shape[1] < sm.d and steps < 400000:
        sm.step(); steps += 1
    sm.fix_drift(); chk = E.check_form(sm.A, B / sm.lam1)
    out.write(json.dumps(dict(V=sm.V(), my_V=chk['V'], my_minq=chk['minq'], my_nmin=chk['nmin'], my_rank=chk['rank'], t=sm.t, warn=sm.warn)) + '\n'); out.flush(); k += 1
print('zn robust n', n, 'runs', k, 'time', round(time.time() - t0, 1))
