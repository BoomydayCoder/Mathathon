#!/usr/bin/env python3
# extras: delta_u matrices (action of the order O), exact R_U of the line through a violator, distinct violator counts, unit-twisted eigenlattices
import sys, json, math, os, time
import numpy as np
from weil import *
p = int(sys.argv[1]); mode = sys.argv[2]
Wl, out = run(p, []); n, m, Dn = Wl.n, Wl.m, Wl.Dn
res = {}
g = Wl.g; ginv2 = pow(pow(g, 2, p), p - 2, p)
Mdg = np.zeros((n, n), dtype=np.int64)
for j in range(1, p):
    Mdg[j - 1, (p - 1) + (j * ginv2) % p - 1] = 1; Mdg[(p - 1) + j - 1, (j * ginv2) % p - 1] = -1
def delta_R(us):
    rows = []
    for c in (0, 1):
        x = us.copy() if c == 0 else us @ Mdg
        for j in range(1, p):
            x = x @ Wl.MT; rows.append(x)
    return np.array(rows, dtype=np.int64)
US = Wl.CH @ Wl.Ar                      # scaled coords of basis of Lambda^H
RS = [delta_R(u) for u in US]
v0s = np.zeros(n, dtype=np.int64); v0s[:p - 1] = -Dn
assert np.array_equal(delta_R(v0s), Dn * np.eye(n, dtype=np.int64)), 'delta_{v0} != id'
for R in RS:
    assert np.array_equal(R @ Wl.NS, Wl.NS @ R) and np.array_equal(R @ Wl.MT, Wl.MT @ R), 'delta_u not G-equivariant'
PD = [Wl.isometry(exdiv(Wl.Ar @ R, Dn)) for R in RS]     # delta_{u_i} as integer matrices on the saved basis (=> O preserves Lambda)
log('delta_u matrices: delta_{v0}=1, commute with S,T (exact), preserve Lambda_0 (exact) for a Z-basis of Lambda^H  => O = End_G(Lambda_0) acts')
def smith_prod(M):
    M = [[int(x) for x in r] for r in M]; R, Cn = len(M), len(M[0]); prod = 1
    for t in range(R):
        while True:
            best = None
            for i in range(t, R):
                for j in range(t, Cn):
                    if M[i][j] != 0 and (best is None or abs(M[i][j]) < abs(M[best[0]][best[1]])): best = (i, j)
            if best is None: return 0
            i, j = best; M[t], M[i] = M[i], M[t]
            for r in M: r[t], r[j] = r[j], r[t]
            pv = M[t][t]; clean = True
            for i in range(t + 1, R):
                q = M[i][t] // pv
                if q: M[i] = [a - q * b for a, b in zip(M[i], M[t])]
                if M[i][t]: clean = False
            for j in range(t + 1, Cn):
                q = M[t][j] // pv
                if q:
                    for r in M: r[j] -= q * r[t]
                if M[t][j]: clean = False
            if clean: break
        prod *= abs(M[t][t])
    return prod
def line_RU(w, label):
    w = np.array(w, dtype=np.int64)
    Ow = np.array([w @ P for P in PD])                 # O w  (8 x n, saved-basis coords)
    Bs = Ow @ Wl.B @ Ow.T; d = det_exact(Bs)
    b = int(w @ Wl.B @ w); e = int(w @ Wl.E @ w); N = p * b * b - e * e
    idx = smith_prod(Ow)
    pred = N ** 4; ok = (d * 16 ** 4 == pred)
    num = m ** 8 * Wl.Dexact ** 4 * idx * idx
    r = dict(label=label, b=b, e=e, N=N, det_Ow=str(d), det_formula_ok=ok, saturation_index=idx, R_U_exact='%d/%d' % (num, d), R_U=num / d, lower_bound=(Wl.theta * 4 * Wl.D * p / N) ** 4)
    log('LINE', r); return r
if mode == 'line':
    res['v0'] = line_RU((-np.ones(1, dtype=np.int64) * 0).tolist() * 0 or Wl.isometry(v0s[None, :])[0], 'v0')
    for cf in sys.argv[3:]:
        c = json.load(open(os.path.join(HERE, cf))); res[cf] = line_RU(c['w_coords_in_saved_basis'], cf)
if mode == 'count':
    q = int(sys.argv[3]); C = Wl.fixed_perm(q); Bs, Es = Wl.grams(C); k = len(C); sp = math.sqrt(p)
    le = math.log(Wl.eps); step = 0.1; nt = int(math.ceil(2 * le / step)); bound = 2 * math.sqrt(Wl.theta * Wl.D) * math.cosh(step / 2) * (1 + 1e-9)
    thr = Wl.Dexact * (p - 1) ** 2; found = {}; U = np.eye(k, dtype=np.int64)
    for i in range(nt):
        t = -le + step * (i + 0.5)
        def Qof(U):
            Bc = (U @ Bs @ U.T).astype(float); Ec = (U @ Es @ U.T).astype(float)
            return math.exp(t) * (Bc + Ec / sp) / 2 + math.exp(-t) * (Bc - Ec / sp) / 2
        Q = Qof(U); U2, _ = lll_int(np.rint(Q * (2.0 ** 36 / np.abs(Q).max())).astype(np.int64)); U = U2 @ U; Q = Qof(U)
        X, vals, flag, nodes = fp(Q, bound); assert flag == 0
        Wv = X @ U; b, e, N = Wl.gamma_data(Wv, Bs, Es)
        for j in np.nonzero(N <= thr)[0]:
            tstar = 0.5 * math.log((b[j] - e[j] / sp) / (b[j] + e[j] / sp))
            if -le <= tstar < le:
                key = tuple(Wv[j]) if tuple(Wv[j]) > tuple(-Wv[j]) else tuple(-Wv[j]); found[key] = (int(b[j]), int(e[j]), int(N[j]))
    from collections import Counter
    cnt = Counter(v[2] for v in found.values())
    tab = sorted((N, c, N / thr) for N, c in cnt.items())
    res['count_C%d' % q] = dict(rank=k, distinct_vectors_mod_pm_and_eps=len(found), table_N_count_ratio=tab, threshold=thr)
    log('COUNT split C_%d fixed lattice (rank %d): distinct vectors mod +-eps^Z with Gamma<=theta: %d ; (N, count, Gamma/theta): %s' % (q, k, len(found), [(a, c, round(r, 5)) for a, c, r in tab]))

if mode == 'twist':
    # unit-twisted eigenlattices ker(g - delta), g = generator of split C_3, delta in O^x of order 3 (not Delta-subspaces, but natural places for short vectors)
    maxd = int(sys.argv[3]); tl = float(sys.argv[4])
    BH, EH = Wl.grams(Wl.CH)
    X, vals, flag, nodes = fp(BH.astype(float), 4 * m * m + 0.5)
    bb, ee, NN = Wl.gamma_data(X, BH, EH)
    units = X[(bb == 4 * m * m) & (ee == 0)]
    I = np.eye(n, dtype=np.int64); ord3 = []
    for x in units:
        for sgn in (1, -1):
            Pd = sum(int(sgn * x[i]) * PD[i] for i in range(8))
            if not np.array_equal(Pd, I) and np.array_equal(Pd @ Pd @ Pd, I): ord3.append(Pd)
    res['n_units_nrd1_mod_pm'] = int(len(units)); res['n_order3'] = len(ord3)
    log('units of O with nrd=1 (b=4m^2,e=0) mod +-: %d ; elements of order 3 in O^x: %d' % (len(units), len(ord3)))
    Pg = Wl.torus_perm(3); assert np.array_equal(Pg @ Pg @ Pg, I)
    k3 = len(Wl.fixed_perm(3)); expected = (n - k3) // 2
    res['twist'] = []
    for di, Pd in enumerate(ord3[:maxd]):
        assert np.array_equal(Pd @ Pg, Pg @ Pd)
        K = Pg - Pd
        Gk = Wl.B + 10 ** 5 * (K @ Wl.B @ K.T)
        U, _ = lll_int(Gk, 0.99)
        ker = np.array([not np.any(U[i] @ K) for i in range(n)])
        log('twist %d: kernel rank found %d, expected %d' % (di, int(ker.sum()), expected))
        if int(ker.sum()) != expected: continue
        C = U[ker]; U2, _ = lll_int(C @ Wl.B @ C.T, 0.99); C = U2 @ C
        Bs, Es = Wl.grams(C); sgn_, ld = np.linalg.slogdet(Bs.astype(float))
        kk = len(C); Rlike = math.exp(kk * math.log(n / kk) + (kk / 2) * math.log(Wl.D) - ld)
        r = Wl.tsearch(C, label='p%d twisted eigenlattice ker(g3 - delta3) #%d' % (p, di), maxnodes=int(2e9), tlim=tl)
        r['R_like'] = Rlike; res['twist'].append(r)
        log('   (density indicator (n/k)^k D^(k/2)/det = %.6g)' % Rlike)
        json.dump(res, open(os.path.join(HERE, 'extras_p%d_twist.json' % p), 'w'), indent=1, default=str)

if mode == 'lines':
    # exact R_U of the Delta-line through EVERY vector w with Gamma(w) <= M*theta in the split C_q-fixed sublattice (handles non-principal lines up to index^2 <= M^4)
    from fractions import Fraction
    q = int(sys.argv[3]); Mf = float(sys.argv[4]); cap = int(sys.argv[5]); C = Wl.fixed_perm(q); Bs, Es = Wl.grams(C); k = len(C); sp = math.sqrt(p)
    le = math.log(Wl.eps); step = 0.1; nt = int(math.ceil(2 * le / step)); bound = 2 * math.sqrt(Mf * Wl.theta * Wl.D) * math.cosh(step / 2) * (1 + 1e-9)
    thr = Wl.Dexact * (p - 1) ** 2; found = {}; U = np.eye(k, dtype=np.int64); complete = True
    for i in range(nt):
        t = -le + step * (i + 0.5)
        def Qof(U):
            Bc = (U @ Bs @ U.T).astype(float); Ec = (U @ Es @ U.T).astype(float)
            return math.exp(t) * (Bc + Ec / sp) / 2 + math.exp(-t) * (Bc - Ec / sp) / 2
        Q = Qof(U); U2, _ = lll_int(np.rint(Q * (2.0 ** 36 / np.abs(Q).max())).astype(np.int64)); U = U2 @ U; Q = Qof(U)
        X, vals, flag, nodes = fp(Q, bound, maxc=2000000)
        if flag: complete = False
        Wv = X @ U; b, e, N = Wl.gamma_data(Wv, Bs, Es)
        for j in np.nonzero(N <= Mf * thr)[0]:
            tstar = 0.5 * math.log((b[j] - e[j] / sp) / (b[j] + e[j] / sp))
            if -le <= tstar < le:
                key = tuple(Wv[j]) if tuple(Wv[j]) > tuple(-Wv[j]) else tuple(-Wv[j]); found[key] = int(N[j])
    items = sorted(found.items(), key=lambda kv: kv[1]); log('lines: %d distinct vectors mod +-eps^Z with Gamma <= %.2f theta (enumeration complete=%s); processing %d' % (len(items), Mf, complete, min(cap, len(items))))
    from collections import Counter
    stat = Counter(); maxR = Fraction(0); argmax = None
    for key, N in items[:cap]:
        w = np.array(key, dtype=np.int64) @ C
        Ow = np.array([w @ P for P in PD]); idx = smith_prod(Ow); R = Fraction(thr ** 4 * idx * idx, N ** 4)
        stat[(round(N / thr, 5), idx, 'R>1' if R > 1 else ('R=1' if R == 1 else 'R<1'))] += 1
        if R > maxR: maxR = R; argmax = (w.tolist(), N, idx)
    res['lines_C%d' % q] = dict(rank=k, M=Mf, complete=complete, n_vectors=len(items), n_processed=min(cap, len(items)), max_R_U=str(maxR), max_R_U_float=float(maxR),
                                 stats=[(list(kk), v) for kk, v in sorted(stat.items())], argmax=argmax)
    log('LINES split C_%d (rank %d), all w with Gamma<=%.2f theta: max exact R_U = %s = %.6f ; (Gamma/theta, saturation index, class): count = %s' % (q, k, Mf, maxR, float(maxR), sorted(stat.items())))
json.dump(res, open(os.path.join(HERE, 'extras_p%d_%s.json' % (p, '_'.join(sys.argv[2:4]).replace('.json', ''))), 'w'), indent=1, default=str)
log('DONE')
