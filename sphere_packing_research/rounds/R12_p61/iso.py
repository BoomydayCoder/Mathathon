#!/usr/bin/env python3
# R12 tester: tensor-stability test on isotypic components of subgroups H < SL_2(p) acting on the Weil lattice.
# Exact integer core (reuses builder B's weil.py: Weil, hnf_mod, lll_int, det_exact).
# KEY FORMULA (derived in notes.md): Lambda^# = (eta/(r p)) Lambda  (verified exactly below)  =>
#   det_B(Lambda cap U) = D^{k/2} * g_U,   g_U = [Lambda : (Lambda cap U) + (Lambda cap U^perp)]   (glue index)
#   R_U = (n/k)^k / g_U ,  g_U = |rowspace(Et) mod N| where e_U = Et/N (integer matrix Et).
import sys, os, math, json, time, itertools
from math import gcd
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from weil import *

def mmul(A, B, p):
    return ((A[0]*B[0]+A[1]*B[2]) % p, (A[0]*B[1]+A[1]*B[3]) % p, (A[2]*B[0]+A[3]*B[2]) % p, (A[2]*B[1]+A[3]*B[3]) % p)
def mpow(A, k, p):
    X = (1, 0, 0, 1)
    for _ in range(k): X = mmul(X, A, p)
    return X
def morder(A, p):
    I = (1, 0, 0, 1); X = A; k = 1
    while X != I: X = mmul(X, A, p); k += 1
    return k
def closure(gens, p, cap=5000):
    I = (1, 0, 0, 1); seen = {I}; fr = [I]
    while fr:
        new = []
        for x in fr:
            for g in gens:
                y = mmul(x, g, p)
                if y not in seen: seen.add(y); new.append(y)
        fr = new
        assert len(seen) <= cap
    return sorted(seen)

class Rep:
    def __init__(s, Wl):
        s.Wl = Wl; s.p = Wl.p; n = Wl.n
        s.PS = Wl.word('S'); s.PT = Wl.word('T')
        I = np.eye(n, dtype=np.int64); s.I = I
        assert np.array_equal(s.PS @ s.PS, -I)
        s.Tpow = [I]
        for k in range(1, s.p): s.Tpow.append(s.Tpow[-1] @ s.PT)
        assert np.array_equal(s.Tpow[-1] @ s.PT, I)
        ST = s.PS @ s.PT; assert np.array_equal(ST @ ST @ ST, -I)
        # isometry checks
        assert np.array_equal(s.PS @ Wl.B @ s.PS.T, Wl.B) and np.array_equal(s.PT @ Wl.B @ s.PT.T, Wl.B)
        assert np.array_equal(s.PS @ Wl.E @ s.PS.T, Wl.E) and np.array_equal(s.PT @ Wl.E @ s.PT.T, Wl.E)
    def op(s, M):
        p = s.p; a, b, c, d = M
        if c == 0:
            return s.op(mmul(M, (0, 1, p - 1, 0), p)) @ s.PS
        ci = pow(c, p - 2, p); x = (a - 1) * ci % p; z = (d - 1) * ci % p
        return s.Tpow[x] @ s.PS @ s.Tpow[(-c) % p] @ (-s.PS) @ s.Tpow[z]
    def test_hom(s, ntest=6, seed=1):
        rng = np.random.default_rng(seed); p = s.p
        def rnd():
            X = (1, 0, 0, 1)
            for _ in range(6):
                X = mmul(X, mpow((1, 1, 0, 1), int(rng.integers(0, p)), p), p); X = mmul(X, (0, p - 1, 1, 0), p)
            return X
        for _ in range(ntest):
            A, B = rnd(), rnd()
            assert np.array_equal(s.op(A) @ s.op(B), s.op(mmul(A, B, p))), 'not a homomorphism'
        # also c == 0 elements
        A = (3, 5, 0, pow(3, p - 2, p)); B = rnd()
        assert np.array_equal(s.op(A) @ s.op(B), s.op(mmul(A, B, p)))
        assert np.array_equal(s.op((p - 1, 0, 0, p - 1)), -s.I)
        return True

def setup(p, use_saved=True):
    Wl = Weil(p)
    Bdir = '/home/user/Mathathon/sphere_packing_research/rounds/R10_weil_lattices/B/'
    f = os.path.join(Bdir, 'Ared_p%d.npy' % p)
    if use_saved and os.path.exists(f):
        Wl.Ar = np.load(f); Wl.B = np.load(os.path.join(Bdir, 'B_p%d.npy' % p)); Wl.E = np.load(os.path.join(Bdir, 'E_p%d.npy' % p))
        Dn = Wl.Dn
        assert np.array_equal(exdiv(Wl.Ar @ Wl.GB @ Wl.Ar.T, Dn * Dn), Wl.B) and np.array_equal(exdiv(Wl.Ar @ Wl.EB @ Wl.Ar.T, Dn * Dn), Wl.E)
        Wl.Dexact = Wl.m ** 2 * p; Wl.D = float(Wl.Dexact)
        sgn, ld = np.linalg.slogdet(Wl.B.astype(float))
        log('loaded saved basis; float logdet(B)/(p-1) = %.9f vs log D = %.9f' % (ld / (p - 1), math.log(Wl.D)))
        assert abs(ld / (p - 1) - math.log(Wl.D)) < 1e-8
    else:
        Wl.build()
    return Wl

def check_modular(Wl):
    """exact check  B X B = r E with X integral  (=> Lambda^# = eta/(r p) Lambda, given det identities)"""
    n = Wl.n; r = Wl.m; Bf = Wl.B.astype(np.longdouble) if False else Wl.B.astype(float); Ef = Wl.E.astype(float)
    Y = np.linalg.solve(Bf, Ef)            # B^-1 E
    Xf = r * np.linalg.solve(Bf, Y.T).T    # r B^-1 E B^-1
    X = np.rint(Xf); res = np.abs(Xf - X).max()
    Xo = X.astype(np.int64).astype(object); Bo = Wl.B.astype(object); Eo = Wl.E.astype(object)
    ok = np.array_equal(Bo @ Xo @ Bo, r * Eo)
    sgn, ld = np.linalg.slogdet(X)
    log('modularity check: rounding residual %.2e, exact B X B == r E: %s, float log|det X| = %.3e' % (res, ok, ld))
    return ok and abs(ld) < 1e-6

def ramanujan(d, k):
    # c_d(k) = sum over primitive d-th roots of unity z of z^k
    s = 0
    for j in range(d):
        if gcd(j, d) == 1: s += np.exp(2j * np.pi * j * k / d)
    v = int(round(s.real)); assert abs(s - v) < 1e-9
    return v

def factor_smooth(g, primes=(2, 3, 5, 7, 11, 13, 31, 37, 61)):
    out = {}
    for q in primes:
        while g % q == 0: g //= q; out[q] = out.get(q, 0) + 1
    return out, g

def comp_stats(Wl, Et, N, label, results, extra=None):
    """Et integer matrix with e_U = Et/N (row-coordinate action).  Returns dict with exact glue index and R_U."""
    n = Wl.n; I = np.eye(n, dtype=np.int64)
    assert np.array_equal(Et @ Et, N * Et), 'not idempotent'
    tr = int(np.trace(Et)); assert tr % N == 0; k = tr // N
    assert np.array_equal(Et @ Wl.B, (Et @ Wl.B).T), 'projector not B-self-adjoint'
    assert np.array_equal(Et @ Wl.E, (Et @ Wl.E).T), 'projector not E-self-adjoint (not eta-stable)'
    if k == 0 or k == n:
        return None
    Hn = hnf_mod(Et, N); prod = 1
    for v in np.diag(Hn): prod *= int(v)
    assert N ** n % prod == 0
    g = N ** n // prod
    logR = k * math.log(n / k) - math.log(g)
    logRc = (n - k) * math.log(n / (n - k)) - math.log(g)
    fac, rest = factor_smooth(g)
    r = dict(p=Wl.p, label=label, dimU=k, N=N, glue_index=str(g), glue_factor=fac, glue_rest=str(rest), logR=logR, R=math.exp(logR),
             logR_complement=logRc, R_complement=math.exp(logRc) if logRc < 700 else float('inf'), delta_rank=k / 8)
    if extra: r.update(extra)
    log('COMP p=%d %-34s dim %3d  glue g=%s  log R_U=%+.6f R_U=%.6g   [complement dim %d: log R=%+.4f]' %
        (Wl.p, label, k, '*'.join('%d^%d' % (q, e) for q, e in fac.items()) + ('' if rest == 1 else '*%d' % rest), logR, math.exp(logR), n - k, logRc))
    results.append(r)
    return r

def sat_basis(Wl, Et, N, k, g, lam=10**6):
    """LLL-reduced basis (rows, coords in reduced basis) of Lambda cap U, certified by det = g * D^{k/2}."""
    n = Wl.n; K = N * np.eye(n, dtype=np.int64) - Et
    gg = int(np.gcd.reduce(K.ravel())); K = K // gg
    KK = K @ K.T
    lam = int(min(lam, 2 ** 58 // int(np.abs(KK).max())))
    Gk = np.eye(n, dtype=np.int64) + lam * KK
    U, _ = lll_int(Gk, 0.99)
    ker = [i for i in range(n) if not np.any(U[i] @ K)]
    assert len(ker) == k, ('kernel rank mismatch', len(ker), k)
    C = U[ker]
    U2, Bs = lll_int(C @ Wl.B @ C.T, 0.99); C = U2 @ C
    sgn, ld = np.linalg.slogdet(Bs.astype(float))
    target = math.log(g) + (k / 2) * math.log(Wl.D)
    ok = abs(ld - target) < 1e-6
    exact = None
    if k <= 40:
        d = det_exact(Bs); exact = (d == g * Wl.Dexact ** (k // 2)) if k % 2 == 0 else None
    log('   saturated basis rank %d: float logdet %.6f vs log(g D^{k/2}) = %.6f ok=%s exact=%s; diag(B) min %d max %d' % (k, ld, target, ok, exact, Bs.diagonal().min(), Bs.diagonal().max()))
    assert ok
    return C, exact

def cyclic_components(rep, A, name, results, sums=True, want=None):
    """A 2x2 matrix generating cyclic group; all Q-isotypic components (by exact order of eigen-character) and (optionally) all sums."""
    Wl = rep.Wl; p = Wl.p; n = Wl.n; N = morder(A, p)
    P = rep.op(A); pw = [np.eye(n, dtype=np.int64)]
    for k in range(1, N): pw.append(pw[-1] @ P)
    assert np.array_equal(pw[-1] @ P, np.eye(n, dtype=np.int64))
    atoms = {}
    for d in range(1, N + 1):
        if N % d: continue
        Et = sum(ramanujan(d, k) * pw[k] for k in range(N))
        tr = int(np.trace(Et)); assert tr % N == 0
        if tr == 0: continue
        atoms[d] = Et
    log('%s (cyclic order %d): eigen-orders present with dims: %s' % (name, N, {d: int(np.trace(E)) // N for d, E in atoms.items()}))
    out = {}
    ds = sorted(atoms)
    subsets = [c for r_ in range(1, len(ds)) for c in itertools.combinations(ds, r_)] if sums else [(d,) for d in ds]
    for c in subsets:
        Et = sum(atoms[d] for d in c)
        lab = '%s: eigen-orders {%s}' % (name, ','.join(map(str, c)))
        r = comp_stats(Wl, Et, N, lab, results, extra=dict(subgroup=name, orders=list(c)))
        out[c] = (Et, N, r)
    return out

def group_components(rep, gens, name, comps, results, expect_order=None, relcheck=None):
    Wl = rep.Wl; p = Wl.p; n = Wl.n
    els = closure(gens, p); Hn = len(els)
    if expect_order: assert Hn == expect_order, (name, Hn)
    ords = [morder(M, p) for M in els]
    ops = [rep.op(M) for M in els]
    hist = {}
    for o in ords: hist[o] = hist.get(o, 0) + 1
    log('%s: order %d, element orders %s' % (name, Hn, dict(sorted(hist.items()))))
    if relcheck: relcheck(rep, gens)
    out = {}
    for lab, (deg, val) in comps.items():
        Et = sum(val[o] * P for o, P in zip(ords, ops))
        gg = gcd(deg, Hn); Et = (deg // gg) * Et; N = Hn // gg
        r = comp_stats(Wl, Et, N, '%s: %s' % (name, lab), results, extra=dict(subgroup=name, component=lab))
        out[lab] = (Et, N, r)
    return out

# character data (trace over the Galois orbit, as function of element order; all elements of H here are +-1 or split/non-split semisimple)
CH_2I = {'icosian 2a+2b (Q-irr dim 8)': (2, {1: 4, 2: -4, 4: 0, 3: -2, 6: 2, 5: -1, 10: 1}),
         'degree-4 faithful': (4, {1: 4, 2: -4, 4: 0, 3: 1, 6: -1, 5: -1, 10: 1}),
         'degree-6 faithful': (6, {1: 6, 2: -6, 4: 0, 3: 0, 6: 0, 5: 1, 10: -1})}
CH_2T = {'Hurwitz 2 (rational char)': (2, {1: 2, 2: -2, 4: 0, 3: -1, 6: 1}),
         "2'+2'' (Q(omega))": (2, {1: 4, 2: -4, 4: 0, 3: 1, 6: -1})}

def find_gens(p, kind):
    """generators s,t of binary polyhedral group <s,t | (st)^2 = s^3 = t^k = -1> inside SL_2(F_p); k=3 (2T), 5 (2I)"""
    k = {'2T': 3, '2I': 5}[kind]
    if k == 3: taus = [1]
    else:
        r5 = [x for x in range(p) if (x * x - 5) % p == 0]; assert r5
        taus = [((1 + r5[0]) * pow(2, p - 2, p)) % p]
    out = []
    for tau in taus:
        t = (0, p - 1, 1, tau)   # trace tau, det 1
        for a in range(p):
            for b in range(p):
                # s = [[a,b],[c,1-a]], det = a(1-a) - bc = 1
                rhs = (a * (1 - a) - 1) % p
                if b == 0:
                    if rhs != 0: continue
                    cs = range(p)
                else: cs = [rhs * pow(b, p - 2, p) % p]
                for c in cs:
                    s_ = (a, b, c, (1 - a) % p)
                    st = mmul(s_, t, p)
                    if (st[0] + st[3]) % p == 0:
                        out.append((s_, t))
                        if len(out) >= 3: return out
    return out

def relcheck_poly(k):
    def f(rep, gens):
        s_, t = gens; I = rep.I
        Ps, Pt = rep.op(s_), rep.op(t)
        st = Ps @ Pt
        assert np.array_equal(st @ st, -I) and np.array_equal(Ps @ Ps @ Ps, -I)
        X = I
        for _ in range(k): X = X @ Pt
        assert np.array_equal(X, -I)
        log('   operator relations (st)^2 = s^3 = t^%d = -1 verified exactly on Lambda' % k)
    return f

def run(p, do_search=True, tl=120.0):
    results = []; searches = []
    Wl = setup(p); n = Wl.n
    rep = Rep(Wl); rep.test_hom()
    log('p=%d: homomorphism SL_2(F_p) -> GL(Lambda) verified on random products' % p)
    okm = check_modular(Wl); assert okm
    g0 = Wl.g
    # split torus
    Asplit = (g0, 0, 0, pow(g0, p - 2, p))
    split = cyclic_components(rep, Asplit, 'split torus C_%d' % (p - 1), results, sums=True)
    # non-split torus: find element of order p+1
    rng = np.random.default_rng(5); Ans = None
    for _ in range(100000):
        a, b, c = [int(v) for v in rng.integers(0, p, 3)]
        if b == 0: continue
        # d from det: a d - b c = 1
        if a == 0: continue
        d = (1 + b * c) * pow(a, p - 2, p) % p
        M = (a, b, c, d)
        tr = (a + d) % p
        if legendre(tr * tr - 4, p) == -1 and morder(M, p) == p + 1: Ans = M; break
    ns = cyclic_components(rep, Ans, 'non-split torus C_%d' % (p + 1), results, sums=False)
    # unipotent
    un = cyclic_components(rep, (1, 1, 0, 1), 'unipotent C_%d' % p, results, sums=False)
    poly = {}
    gl = find_gens(p, '2T')
    for i, gens in enumerate(gl[:1]):
        poly['2T'] = group_components(rep, gens, 'SL_2(3)', CH_2T, results, expect_order=24, relcheck=relcheck_poly(3))
    if p % 10 in (1, 9):
        gl = find_gens(p, '2I')
        for i, gens in enumerate(gl[:1]):
            poly['2I'] = group_components(rep, gens, 'SL_2(5)', CH_2I, results, expect_order=120, relcheck=relcheck_poly(5))
    json.dump(results, open(os.path.join(HERE, 'iso_p%d%s.json' % (p, '' if do_search else '_ns')), 'w'), indent=1)
    if do_search:
        todo = []
        for key in poly:
            for lab, (Et, N, r) in poly[key].items():
                if r and r['dimU'] <= 40: todo.append((r['label'], Et, N, r))
        for c, (Et, N, r) in split.items():
            if r and r['dimU'] <= 40: todo.append((r['label'], Et, N, r))
        todo.sort(key=lambda x: (0 if 'SL_2' in x[0] else 1, x[3]['dimU']))
        for lab, Et, N, r in todo:
            try:
                C, exact = sat_basis(Wl, Et, N, r['dimU'], int(r['glue_index']))
            except AssertionError as ex:
                log('   sat_basis FAILED for %s: %s' % (lab, ex)); continue
            sr = Wl.tsearch(C, label='p%d %s' % (p, lab), tlim=tl, maxnodes=2 * 10**9)
            sr['det_certified_exact'] = exact
            searches.append(sr)
            json.dump(searches, open(os.path.join(HERE, 'search_p%d.json' % p), 'w'), indent=1, default=str)
    log('DONE p=%d' % p)
    return Wl, rep, results, searches

if __name__ == '__main__':
    p = int(sys.argv[1]); ds = (len(sys.argv) < 3 or sys.argv[2] != 'nosearch'); tl = float(sys.argv[3]) if len(sys.argv) > 3 else 120.0
    run(p, ds, tl)
