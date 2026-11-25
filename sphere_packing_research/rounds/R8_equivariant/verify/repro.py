"""(b) reproduce runs with the original code (same seeds), capture the FIRST Q-perfect form (not stored in the raw data)
and verify it with my own enumerator; compare V_perfect / pairs_perfect with the stored raw data."""
import numpy as np, json, sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
HERE = os.path.dirname(os.path.abspath(__file__)); R8 = os.path.dirname(HERE)
sys.path.insert(0, R8)
import eq as their
cap = []
orig = their.Eq.greedy
def g2(self, c, S, rng, maxit=200):
    cap.append((self.mat(c).copy(), len(S)))
    return orig(self, c, S, rng, maxit)
their.Eq.greedy = g2
os.makedirs(os.path.join(HERE, 'repro'), exist_ok=True)
cases = [l.split(':') for l in sys.argv[2:]]
res = []
for name, k, mode, seed, nr in cases:
    k = int(k); seed = int(seed); nr = int(nr)
    cap.clear()
    t0 = time.time()
    r = their.run_family(name, k, nr, 1e9, mode, seed, outdir=os.path.join(HERE, 'repro'))
    st = json.load(open(os.path.join(R8, 'data', f'{name}_k{k}_{mode}_s{seed}.json')))
    n, D, w = r['n'], r['D'], r['w']
    els, wm, n0, extra = vlib.base_group(name)
    QB = vlib.Qbasis_reynolds(els, k)
    gs = [np.kron(np.eye(k), g.astype(float)) for g in els]
    nrep = 0; nbad = 0; maxrel = 0
    for i, (x, (A, npairs)) in enumerate(zip(r['recs'], cap)):
        s = st['recs'][i]
        same = abs(x['V_perfect'] / s['V_perfect'] - 1) < 1e-9 and x['pairs_perfect'] == s['pairs_perfect'] and \
               abs(x.get('V_ext', 0) / s.get('V_ext', 1) - 1) < 1e-9
        nrep += same
        # my verification of the first perfect form
        inv = max(np.abs(g.T @ A @ g - A).max() for g in gs)
        U = vlib.my_lll(A); Uf = U.astype(float)
        Ar = Uf.T @ A @ Uf
        Y, _ = vlib.enum((Ar + Ar.T) / 2, 1 + 1e-8)
        X = Y @ U.T; Xf = X.astype(float)
        v = np.einsum('ij,jk,ik->i', Xf, A, Xf)
        V = vlib.kappa(n) * math.exp(-0.5 * np.linalg.slogdet(A)[1])
        Phi = np.einsum('mi,dij,mj->md', Xf, QB, Xf, optimize=True)
        sv = np.linalg.svd(Phi, compute_uv=False)
        rank = int((sv > 1e-8 * sv[0]).sum())
        ok = (inv < 1e-9 and v.min() > 1 - 1e-9 and len(v) == s['pairs_perfect'] and abs(V / s['V_perfect'] - 1) < 1e-9 and rank == D
              and V <= vlib.Vbest(n) * (1 + 1e-9))
        maxrel = max(maxrel, abs(V / s['V_perfect'] - 1))
        nbad += (not ok)
        if not ok:
            print('   BAD', name, k, i, inv, v.min(), len(v), s['pairs_perfect'], V, s['V_perfect'], rank)
    d = dict(name=name, k=k, mode=mode, seed=seed, runs=len(r['recs']), reproduced=int(nrep), first_perfect_verified=int(len(r['recs']) - nbad),
             max_V_relerr=float(maxrel), time=time.time() - t0)
    res.append(d)
    print('REPRO', json.dumps(d), flush=True)
    json.dump(res, open(os.path.join(HERE, f'repro_{sys.argv[1]}.json'), 'w'), indent=1)
