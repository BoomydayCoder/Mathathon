#!/usr/bin/env python3
# ADVERSARIAL VERIFIER (R12_p61/verify): independent certificate of R_U for isotypic components.
# Reuses from the tester ONLY: setup() (loads builder B's saved B,E + S,T integer matrices), hnf_mod and sat_basis (to PRODUCE a candidate basis).
# Everything certified here by own code: group words, relations, projector, membership, primitivity (rank mod q), exact det, exact R_U.
import sys, math, time, json
from fractions import Fraction
import numpy as np
R12 = '/home/user/Mathathon/sphere_packing_research/rounds/R12_p61'
sys.path.insert(0, R12)
import weil, iso
p = int(sys.argv[1]); t0 = time.time()
def P(*a): print('[%6.1fs]' % (time.time() - t0), *a, flush=True)
Wl = iso.setup(p); n = Wl.n; r = Wl.m; D = r * r * p
B = Wl.B; E = Wl.E; I = np.eye(n, dtype=np.int64)
PS = Wl.word('S'); PT = Wl.word('T')
assert (PS @ PS == -I).all(); ST = PS @ PT; assert (ST @ ST @ ST == -I).all()
for X in (PS, PT):
    assert (X @ B @ X.T == B).all() and (X @ E @ X.T == E).all()
P('p=%d n=%d r=%d D=%d  tr(S)=%d tr(T)=%d (expect 0, -2)' % (p, n, r, D, np.trace(PS), np.trace(PT)))
def bareiss(M):
    M = [[int(x) for x in row] for row in M]; m = len(M); prev = 1; sign = 1
    for k in range(m - 1):
        if M[k][k] == 0:
            for i in range(k + 1, m):
                if M[i][k] != 0: M[k], M[i] = M[i], M[k]; sign = -sign; break
            else: return 0
        Mk = M[k]; pk = Mk[k]
        for i in range(k + 1, m):
            Mi = M[i]; a = Mi[k]
            if a == 0 and prev == pk: continue
            for j in range(k + 1, m): Mi[j] = (Mi[j] * pk - a * Mk[j]) // prev
        prev = pk
    return sign * M[m - 1][m - 1]
def rank_mod(C, q):
    A = [[int(x) % q for x in row] for row in C]; rk = 0; rows = len(A); cols = len(A[0])
    for c in range(cols):
        piv = None
        for i in range(rk, rows):
            if A[i][c]: piv = i; break
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]; inv = pow(A[rk][c], q - 2, q); A[rk] = [x * inv % q for x in A[rk]]
        for i in range(rows):
            if i != rk and A[i][c]:
                f = A[i][c]; A[i] = [(x - f * y) % q for x, y in zip(A[i], A[rk])]
        rk += 1
        if rk == rows: break
    return rk
def detmod(M, q):
    A = [[int(x) % q for x in row] for row in M]; m = len(A); det = 1
    for c in range(m):
        piv = None
        for i in range(c, m):
            if A[i][c]: piv = i; break
        if piv is None: return 0
        if piv != c: A[c], A[piv] = A[piv], A[c]; det = -det
        det = det * A[c][c] % q; inv = pow(A[c][c], q - 2, q)
        for i in range(c + 1, m):
            if A[i][c]:
                f = A[i][c] * inv % q; A[i] = [(x - f * y) % q for x, y in zip(A[i], A[c])]
    return det % q
# ---- exact det(B) = D^(n/2)
if '--nodet' not in sys.argv:
    for q in (1000003, 998244353, 2147483647):
        assert detmod(B, q) == pow(D, n // 2, q), 'det mod q mismatch'
    P('det(B) == D^(n/2) mod three large primes: OK')
    dB = bareiss(B); P('EXACT det(B) == D^(n/2):', dB == D ** (n // 2))
# ---- eta
Mf = np.linalg.solve(B.astype(float), E.astype(float)).T   # M = E B^-1
M = np.rint(Mf).astype(np.int64); P('eta matrix rounding residual %.2e' % np.abs(Mf - M).max())
Mo = M.astype(object); Bo = B.astype(object); Eo = E.astype(object)
P('EXACT: M B == E:', bool((Mo @ Bo == Eo).all()), ' M^2 == p:', bool((Mo @ Mo == p * I.astype(object)).all()),
  ' M commutes with S,T:', bool((PS @ M == M @ PS).all() and (PT @ M == M @ PT).all()))
Xf = r * np.linalg.solve(B.astype(float), Mf.T).T if False else r * np.linalg.solve(B.astype(float), M.astype(float))
X = np.rint(Xf).astype(np.int64)
P('EXACT: B X == r M with X integral (=> Lambda^# = eta/(rp) Lambda given det(B)=D^(n/2)):', bool((Bo @ X.astype(object) == r * Mo).all()), 'resid %.1e' % np.abs(Xf - X).max())
# ---- group words (own decomposition  M = T^x S T^c S T^z, x=(a+1)/c, z=(d+1)/c)
def mm(A, C): return ((A[0]*C[0]+A[1]*C[2]) % p, (A[0]*C[1]+A[1]*C[3]) % p, (A[2]*C[0]+A[3]*C[2]) % p, (A[2]*C[1]+A[3]*C[3]) % p)
S2 = (0, p - 1, 1, 0); ID = (1, 0, 0, 1)
Tp = [I]
for k in range(1, p): Tp.append(Tp[-1] @ PT)
assert (Tp[-1] @ PT == I).all()
def op(Mx):
    a, b, c, d = Mx
    if c == 0: return op(mm(Mx, S2)) @ (-PS)
    ci = pow(c, p - 2, p); x = (a + 1) * ci % p; z = (d + 1) * ci % p
    w = mm(mm(mm(mm((1, x, 0, 1), S2), (1, c, 0, 1)), S2), (1, z, 0, 1)); assert w == Mx
    return Tp[x] @ PS @ Tp[c] @ PS @ Tp[z]
def order2(Mx):
    k = 1; Y = Mx
    while Y != ID: Y = mm(Y, Mx); k += 1
    return k
g0 = weil.prim_root(p)
def gen_group(s, t, expect):
    """parallel closure: 2x2 element -> operator as the same word in (rho(s), rho(t))"""
    Ps, Pt = op(s), op(t); els = {ID: I}; fr = [ID]
    while fr:
        new = []
        for x in fr:
            for gg, Pg in ((s, Ps), (t, Pt)):
                y = mm(x, gg)
                if y not in els: els[y] = els[x] @ Pg; new.append(y)
        fr = new; assert len(els) <= expect
    assert len(els) == expect
    return els
def poly_group(k, bsel):
    """binary polyhedral <s,t | s^3 = t^k = (st)^2 = -1>: t = -d(a), a of order k; s: trace 1, trace(st) = 0"""
    a = pow(g0, (p - 1) // k, p); ai = pow(a, p - 2, p); t = ((-a) % p, 0, 0, (-ai) % p)
    # s=(u,b,c,1-u): tr(st) = -(a u + ai (1-u)) = 0  => u (a - ai) = -ai
    u = (-ai) * pow((a - ai) % p, p - 2, p) % p; d = (1 - u) % p; bc = (u * d - 1) % p; assert bc != 0
    b = bsel % p; c = bc * pow(b, p - 2, p) % p; s = (u, b, c, d)
    Ps, Pt = op(s), op(t); mI = -I
    A = Ps @ Ps @ Ps; assert (A == mI).all(), 's^3'
    A = np.linalg.matrix_power(Pt, k); assert (A == mI).all(), 't^k'
    A = Ps @ Pt; assert (A @ A == mI).all(), '(st)^2'
    return s, t
def certify(Et, N, label):
    assert (Et @ Et == N * Et).all(), 'idempotent'
    tr = int(np.trace(Et)); k = tr // N; assert tr == N * k
    assert (Et @ B == (Et @ B).T).all(), 'B-selfadjoint'; assert (Et @ E == (Et @ E).T).all(), 'E-selfadjoint'
    assert (Et @ M == M @ Et).all(), 'commutes with eta'
    Hn = weil.hnf_mod(Et, N); prod = 1
    for v in np.diag(Hn): prod *= int(v)
    g = N ** n // prod
    C, _ = iso.sat_basis(Wl, Et, N, k, g)
    C = np.array(C, dtype=np.int64); assert C.shape == (k, n)
    assert not (C @ (N * I - Et)).any(), 'basis not in U'
    Gs = (C.astype(object) @ Bo) @ C.T.astype(object)
    d = bareiss(Gs); assert d > 0
    fac = {}; rest = d
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61):
        while rest % q == 0: rest //= q; fac[q] = fac.get(q, 0) + 1
    assert rest == 1, ('unfactored part', rest)
    prim = all(rank_mod(C, q) == k for q in fac)
    R = Fraction(n, k) ** k * Fraction(D) ** (k // 2) / d if k % 2 == 0 else None
    P('CERT %-40s dim %3d  primitive(rank mod q=k for q in %s): %s  det=%s  R_U = %s = %.6g   tester-g matches: %s' % (
        label, k, sorted(fac), prim, '*'.join('%d^%d' % kv for kv in sorted(fac.items())), R, float(R), d == g * D ** (k // 2)))
    return dict(label=label, dim=k, primitive=prim, R=str(R), Rfloat=float(R), glue_match=bool(d == g * D ** (k // 2))), C
out = []
# ---- split torus fixed lattices
for Nc in [q for q in (3, 5, 9, 15) if (p - 1) % q == 0]:
    a = pow(g0, (p - 1) // Nc, p); Pd = op((a, 0, 0, pow(a, p - 2, p)))
    assert (np.linalg.matrix_power(Pd, Nc) == I).all()
    Et = sum(np.linalg.matrix_power(Pd, j) for j in range(Nc))
    res, C = certify(Et, Nc, 'Fix(split C_%d)' % Nc); out.append(res)
    if p == 37 and Nc == 3:
        for ent in json.load(open(R12 + '/search_p37.json')):
            if '{4,12}' in ent['label'] or 'Hurwitz' in ent['label']:
                w = np.array(ent['w'], dtype=object); b = int(w @ Bo @ w); e = int(w @ Eo @ w)
                Gam = Fraction(b * b * p - e * e, 4 * p * D); th = Fraction((p - 1) ** 2, 4 * p)
                inU = not ((np.array(ent['w'], dtype=np.int64)) @ (Nc * I - Et)).any()
                P('p37 saved vector [%s]: b=%d e=%d Gamma/theta = %s ; in Fix(C_3): %s ; (theta/Gamma)^4 = %.4f' % (ent['label'], b, e, Gam / th, inU, float((th / Gam) ** 4)))
                out.append(dict(label='p37 vec ' + ent['label'], b=b, e=e, ratio=str(Gam / th), inFixC3=bool(inU)))
# ---- SL_2(3) Hurwitz component, SL_2(5) icosian etc.
chi3 = {1: 2, 2: -2, 4: 0, 3: -1, 6: 1}
for bsel in (1, g0):
    s, t = poly_group(3, bsel); els = gen_group(s, t, 24)
    Et = sum(chi3[order2(x)] * Px for x, Px in els.items())
    res, C = certify(Et, 12, 'SL_2(3) [b=%d] Hurwitz comp' % bsel); out.append(res)
    if p == 37:
        for ent in json.load(open(R12 + '/search_p37.json')):
            if 'Hurwitz' in ent['label']:
                inU = not ((np.array(ent['w'], dtype=np.int64)) @ (12 * I - Et)).any(); P('   saved Hurwitz violator lies in THIS Hurwitz component:', inU)
if p % 5 in (1, 4):
    chi_ico = {1: 4, 2: -4, 4: 0, 3: -2, 6: 2, 5: -1, 10: 1}
    chi4 = {1: 4, 2: -4, 4: 0, 3: 1, 6: -1, 5: -1, 10: 1}
    for bsel in (1, g0):
        s, t = poly_group(5, bsel); els = gen_group(s, t, 120)
        Et = sum(chi_ico[order2(x)] * Px for x, Px in els.items())
        res, C = certify(Et, 60, 'SL_2(5) [b=%d] icosian comp' % bsel); out.append(res)
        Et4 = sum(chi4[order2(x)] * Px for x, Px in els.items())
        res, C = certify(Et4, 30, 'SL_2(5) [b=%d] degree-4 comp' % bsel); out.append(res)
json.dump(out, open('verify_p%d.json' % p, 'w'), indent=1)
P('done')
