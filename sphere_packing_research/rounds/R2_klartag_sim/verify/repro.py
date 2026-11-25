#!/usr/bin/env python3
"""Reproduce saved runs of R2 (same code, same seed), capture the final form A and the lattice basis, and re-verify
them with the verifier's own enumeration (myenum.py).  Usage: repro.py n kind seed offset nruns datafile tag"""
import sys, json, math, time
import numpy as np
sys.dont_write_bytecode = True
import sim_r2_copy as S
import myenum as E

n = int(sys.argv[1]); kind = sys.argv[2]; seed = int(sys.argv[3]); offset = int(sys.argv[4]); nruns = int(sys.argv[5])
datafile = sys.argv[6]; tag = sys.argv[7]
saved = [json.loads(l) for l in open(datafile) if l.strip()]
rng = np.random.default_rng(seed)
out = open(f'repro_{tag}.jsonl', 'w')
forms = {}
t0 = time.time()
for i in range(nruns):
    if offset + i >= len(saved): break
    try:
        B = S.make_lattice(n, kind, rng)
        dt = 1.25e-3 / n
        sm = S.Sim(n, B, rng, dt)
        steps = 0
        while sm.Qc.shape[1] < sm.d and steps < 400000:
            sm.step(); steps += 1
        sm.fix_drift()
    except Exception as e:
        out.write(json.dumps(dict(i=offset + i, error=repr(e))) + '\n'); out.flush(); break
    A = sm.A.copy(); Bc = B / sm.lam1
    rec = saved[offset + i]
    chk = E.check_form(A, Bc)
    box = chk['box']
    r = dict(i=offset + i, n=n, seed=seed, lam1_match=bool(abs(sm.lam1 - rec['lam1']) < 1e-12 * rec['lam1']),
             Vsim=sm.V(), Vsaved=rec['Vfinal'], V_reldiff=abs(sm.V() - rec['Vfinal']) / rec['Vfinal'],
             nevents_match=bool(len(sm.events) == len(rec['events'])),
             saved_npairs=rec['final_npairs'], saved_minq=rec['final_minq'],
             my_minq=chk['minq'], my_nmin=chk['nmin'], my_rank=chk['rank'], my_V=chk['V'], my_Vhonest=chk['Vhonest'],
             box_minq=(box['minq'] if box else None), box_nmin=(box['nmin'] if box else None), boxsize=chk['boxsize'],
             symm_err=float(np.abs(A - A.T).max()), eigmin=float(np.linalg.eigvalsh((A + A.T) / 2).min()))
    out.write(json.dumps(r) + '\n'); out.flush()
    forms[f'A_{offset+i}'] = A; forms[f'B_{offset+i}'] = Bc
np.savez_compressed(f'repro_{tag}_forms.npz', **forms)
print(tag, 'done', len(forms) // 2, 'runs, time', round(time.time() - t0, 1), flush=True)
