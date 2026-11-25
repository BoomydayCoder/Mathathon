#!/usr/bin/env python3
# Referee script (independent of builders' code): exact checks of saved Gram data and certificates.
import numpy as np, json, sys, time, traceback
from fractions import Fraction
from math import gcd
R = '/home/user/Mathathon/sphere_packing_research/rounds/R10_weil_lattices/'
p = int(sys.argv[1]); h = (p - 1) // 2; n = 2 * (p - 1); r = (p - 1) // 4
OUT = {}
t0 = time.time()
def log(*a):
    print('[%6.1fs]' % (time.time() - t0), *a, flush=True)
def leg(a):
    v = pow(a % p, (p - 1) // 2, p); return -1 if v == p - 1 else int(v)
def bareiss(M, pivot=False):
    nn = len(M); A = [list(map(int, row)) for row in M]; prev = 1; minors = []; sgn = 1
    for k in range(nn):
        if A[k][k] == 0:
            if not pivot: return 0, minors
            for i in range(k + 1, nn):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]; sgn = -sgn; break
            else: return 0, minors
        akk = A[k][k]; Ak = A[k]; minors.append(akk)
        for i in range(k + 1, nn):
            Ai = A[i]; aik = Ai[k]
            if aik == 0:
                for j in range(k + 1, nn): Ai[j] = (Ai[j] * akk) // prev
            else:
                for j in range(k + 1, nn): Ai[j] = (Ai[j] * akk - aik * Ak[j]) // prev
        prev = akk
    return sgn * A[nn - 1][nn - 1], minors
def iroot(N, k):
    lo, hi = 0, 1 << (N.bit_length() // k + 2)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= N: lo = mid
        else: hi = mid - 1
    return lo
def O(a): return np.array(a).astype(object)
def check_gram(tag, Bf, Ef):
    res = {}
    B = np.load(R + Bf); E = np.load(R + Ef)
    res['shape'] = list(B.shape); res['dtype'] = str(B.dtype)
    Bo = O(B); Eo = O(E)
    res['symmetric'] = bool((Bo == Bo.T).all() and (Eo == Eo.T).all())
    det, minors = bareiss(Bo.tolist())
    res['posdef_exact_leading_minors'] = bool(len(minors) == B.shape[0] and all(m > 0 for m in minors))
    D = iroot(det, p - 1)
    res['det'] = str(det); res['D'] = D; res['det_is_D_pow_p-1'] = bool(D ** (p - 1) == det)
    M = np.rint(np.linalg.solve(B.astype(float), E.astype(float))).astype(np.int64)
    Mo = O(M); I = O(np.eye(B.shape[0], dtype=np.int64))
    res['M=B^-1E_integral_exact(BM==E)'] = bool((Bo.dot(Mo) == Eo).all())
    res['M^2==pI_exact'] = bool((Mo.dot(Mo) == p * I).all())
    res['E B^-1 E == pB (follows)'] = res['M=B^-1E_integral_exact(BM==E)'] and res['M^2==pI_exact'] and bool((Eo.dot(Mo) == p * Bo).all())
    res['(I+M)/2 integral (O_F acts)'] = bool((((Mo + I) % 2) == 0).all())
    res['B even diag'] = bool(all(int(Bo[i, i]) % 2 == 0 for i in range(B.shape[0])))
    res['min diag B'] = int(min(int(Bo[i, i]) for i in range(B.shape[0])))
    log(tag, json.dumps(res)[:900])
    return Bo, Eo, Mo, D, det, res
def viol(tag, Bo, Eo, D, w):
    w = O(w)
    b = int(w.dot(Bo).dot(w)); e = int(w.dot(Eo).dot(w))
    N = p * b * b - e * e; THR = (p - 1) ** 2 * D
    mp = (b + e / np.sqrt(p)) / 2; mm = (b - e / np.sqrt(p)) / 2
    res = dict(b=b, e=e, N=N, THR=THR, violator=bool(N < THR), Gamma=mp * mm / D, theta=(p - 1) ** 2 / (4 * p),
               ratio=str(Fraction(N, THR)), m_plus=mp, m_minus=mm, RU_lower=(THR / N) ** 4)
    log(tag, json.dumps(res))
    return res
# ---------------- A ----------------
try:
    BoA, EoA, MoA, DA, detA, OUT['A_gram'] = check_gram('A gram', 'A/B_p%d_L0.npy' % p, 'A/E_p%d_L0.npy' % p)
except Exception as ex:
    traceback.print_exc(); BoA = None
try:
    BoB, EoB, MoB, DB, detB, OUT['B_gram'] = check_gram('B gram', 'B/B_p%d.npy' % p, 'B/E_p%d.npy' % p)
    if BoA is not None:
        OUT['A_vs_B_same_det'] = bool(detA == detB); log('A vs B same det:', detA == detB)
except Exception as ex:
    log('B gram not available/failed:', repr(ex)); BoB = None; EoB = None; DB = None
# violators
for (who, f, key, Bo_, Eo_, D_) in [('A', 'A/violator_p%d.json' % p, 'x', BoA, EoA, DA if BoA is not None else None),
                                     ('B', 'B/violator_p%d.json' % p, 'w_coords_in_saved_basis', BoB, EoB, DB if BoB is not None else None),
                                     ('B', 'B/violator_p%d_C3.json' % p, 'w_coords_in_saved_basis', BoB, EoB, DB if BoB is not None else None),
                                     ('B', 'B/p%d_ns5_minvec.json' % p, 'w_coords_in_saved_basis', BoB, EoB, DB if BoB is not None else None)]:
    try:
        d = json.load(open(R + f))
    except Exception:
        continue
    if Bo_ is None: continue
    OUT['viol_' + f] = viol('violator ' + f, Bo_, Eo_, D_, d[key])
    if 'Ow_rows' in d:
        Wr = O(d['Ow_rows']); G8 = Wr.dot(Bo_).dot(Wr.T); d8, _ = bareiss(G8.tolist(), pivot=True)
        log('  builder Ow_rows: det_B =', d8, ' r^8 D^4 =', r ** 8 * D_ ** 4, ' R_U>=', float(Fraction(r ** 8 * D_ ** 4, d8)))
# ---------------- independent G-structure check on A's lattice ----------------
try:
    BasC = np.load(R + 'A/BasC_p%d_L0.npy' % p)
    log('BasC shape', BasC.shape)
    x = np.arange(1, h + 1); xy = np.outer(x, x)
    def ops(kappa):
        S = np.zeros((h, h, h), complex); T = np.zeros((h, h, h), complex)
        for ai in range(h):
            a = ai + 1
            z = lambda k: np.exp(2j * np.pi * ((a * k) % p) / p)
            S[ai] = kappa * (z(2 * xy) - z(-2 * xy)) / (leg(a) * np.sqrt(p)); T[ai] = np.diag(z(x * x))
        return S, T
    for kappa in (1, -1):
        S, T = ops(kappa)
        ST = np.einsum('axy,ayz->axz', S, T); ST3 = np.einsum('axy,ayz,azw->axw', ST, ST, ST)
        S2 = np.einsum('axy,ayz->axz', S, S)
        Id = np.broadcast_to(np.eye(h), (h, h, h))
        log('kappa', kappa, '|S^2+I|', abs(S2 + Id).max(), '|(ST)^3-I|', abs(ST3 - Id).max(), '|(ST)^3+I|', abs(ST3 + Id).max())
    S, T = ops(-1)
    def act(Mop, V): return np.einsum('axy,iay->iax', Mop, V)
    def remb(V): f = V.reshape(V.shape[0], -1); return np.concatenate([f.real, f.imag], axis=1)
    Xb = remb(BasC); Pinv = np.linalg.pinv(Xb)
    sv = np.linalg.svd(Xb, compute_uv=False); log('rank of saved basis (real emb):', int((sv > 1e-8 * sv[0]).sum()), 'smallest sv', sv[-1])
    signs = np.array([leg(a) for a in range(1, h + 1)], float)
    Bnum = 2 * Xb @ Xb.T
    W = BasC * signs[None, :, None]; Enum = 2 * np.sqrt(p) * remb(W) @ Xb.T
    OUT['A_gram_matches_embedding'] = [float(abs(Bnum - BoA.astype(float)).max()), float(abs(Enum - EoA.astype(float)).max())]
    log('max|B_num - B_saved|, max|E_num - E_saved| =', OUT['A_gram_matches_embedding'])
    def coords(V):
        c = remb(V) @ Pinv; res = abs(c @ Xb - remb(V)).max(); ci = np.rint(c); return ci.astype(np.int64), float(abs(c - ci).max()), float(res)
    CS, i1, r1 = coords(act(S, BasC)); CT, i2, r2 = coords(act(T, BasC))
    log('S,T images of basis: integrality defect', i1, i2, ' residual (in span)', r1, r2)
    CSo, CTo = O(CS), O(CT); I = O(np.eye(n, dtype=np.int64))
    ex = {}
    ex['S isometry B,E'] = bool((CSo.dot(BoA).dot(CSo.T) == BoA).all() and (CSo.dot(EoA).dot(CSo.T) == EoA).all())
    ex['T isometry B,E'] = bool((CTo.dot(BoA).dot(CTo.T) == BoA).all() and (CTo.dot(EoA).dot(CTo.T) == EoA).all())
    ex['S^2=-I'] = bool((CSo.dot(CSo) == -I).all())
    STo = CSo.dot(CTo); ST3o = STo.dot(STo).dot(STo); ex['(ST)^3=-I'] = bool((ST3o == -I).all()); ex['(ST)^3=+I'] = bool((ST3o == I).all())
    Tp = I
    for _ in range(p): Tp = Tp.dot(CTo)
    ex['T^p=I'] = bool((Tp == I).all())
    MoA = MoA.T  # row-vector convention: eta acts on coordinate rows by M^T
    ex['M^T (eta) commutes with S,T'] = bool((MoA.dot(CSo) == CSo.dot(MoA)).all() and (MoA.dot(CTo) == CTo.dot(MoA)).all())
    OUT['A_G_structure_exact'] = ex; OUT['A_G_integrality_defect'] = [i1, i2, r1, r2]
    log('exact integer checks of S,T on saved basis:', json.dumps(ex))
    # ---- torus element d(a) as a word, fixed lattice F_H
    FH = np.load(R + 'A/F_p%d_H.npy' % p)
    g = next(g for g in range(2, p) if len({pow(g, k, p) for k in range(p - 1)}) == p - 1)
    a0 = pow(g, 4, p) if p % 8 == 5 else None
    if a0 is not None:
        ainv = pow(a0, p - 2, p)
        def Tpow(k): return np.stack([np.diag(np.diag(T[ai]) ** k) for ai in range(h)])
        def mul(*Ms):
            Z = Ms[0]
            for Mx in Ms[1:]: Z = np.einsum('axy,ayz->axz', Z, Mx)
            return Z
        Wd = mul(Tpow(a0), S, Tpow(ainv), S, Tpow(a0), S)
        Dm = np.zeros((h, h))
        for xx in range(1, h + 1):
            y = (a0 * xx) % p
            if y <= h: Dm[xx - 1, y - 1] = 1
            else: Dm[xx - 1, p - y - 1] = -1
        dd = min(abs(Wd - s * Mx[None]).max() for s in (1, -1) for Mx in (Dm, Dm.T))
        log('torus word T^a S T^(1/a) S T^a S equals +-monomial d(a): defect', dd, ' (a=%d of order %d)' % (a0, r))
        VF = np.tensordot(FH.astype(float), BasC, axes=1)
        fx = max(abs(np.einsum('xy,iay->iax', Dm, VF) - VF).max(), 0)
        log('F_H rows fixed by d(a): defect', fx, ' rank F_H', np.linalg.matrix_rank(FH.astype(float)))
        OUT['torus_word_defect'] = float(dd); OUT['F_H_fixed_defect'] = float(fx)
    FHo = O(FH); GH = FHo.dot(BoA).dot(FHo.T); dH, _ = bareiss(GH.tolist(), pivot=True)
    log('det_B(F_H) =', dH, ' r^8 D^4 =', r ** 8 * DA ** 4, ' equal:', dH == r ** 8 * DA ** 4)
    OUT['detB_FH'] = str(dH); OUT['detB_FH_equals_r8D4'] = bool(dH == r ** 8 * DA ** 4)
    # ---- border vector v in Lambda^H
    THR = (p - 1) ** 2 * DA
    best = None
    import itertools
    for cf in itertools.product((-1, 0, 1), repeat=FH.shape[0]):
        if not any(cf): continue
        v = np.array(cf, dtype=np.int64) @ FH
        vo = O(v); b = int(vo.dot(BoA).dot(vo)); e = int(vo.dot(EoA).dot(vo))
        if p * b * b - e * e == THR:
            best = v; log('border vector in F_H lattice: coeffs', cf, 'b', b, 'e', e); break
    if best is None:
        log('no border vector among {-1,0,1} combos of F_H rows; using row 0 (X_j may be non-integral)'); best = FH[0].copy()
    v = best
    # ---- commutant via cyclic vector: X_j (g v) = g u_j
    rng = np.random.default_rng(12345)
    Z = np.vstack([v[None, :], FH]).astype(np.int64)
    Vl = []; Ul = []
    for it in range(6 * n):
        if rng.random() < 0.5: Z = Z @ CS
        else:
            for _ in range(int(rng.integers(1, p))): Z = Z @ CT
        Vl.append(Z[0].copy()); Ul.append(Z[1:].copy())
    Vall = np.array(Vl); Uall = np.array(Ul)  # (N,n), (N,8,n)
    log('orbit sample rank', np.linalg.matrix_rank(Vall.astype(float)), 'max coord', abs(Vall).max())
    Xs = []
    okc = True; maxdef = 0
    for j in range(FH.shape[0]):
        Xf, *_ = np.linalg.lstsq(Vall.astype(float), Uall[:, j, :].astype(float), rcond=None)
        Xi = np.rint(Xf).astype(np.int64); maxdef = max(maxdef, abs(Xf - Xi).max())
        Xo = O(Xi)
        okc &= bool((O(Vall).dot(Xo) == O(Uall[:, j, :])).all())
        okc &= bool((Xo.dot(CSo) == CSo.dot(Xo)).all() and (Xo.dot(CTo) == CTo.dot(Xo)).all())
        Xs.append(Xo)
    log('commutant elements X_j: integrality defect', maxdef, ' exact (X_j commute with S,T; X_j(g v)=g u_j on sample):', okc)
    OUT['commutant_integral_defect'] = float(maxdef); OUT['commutant_exact_ok'] = bool(okc)
    flat = np.array([np.array(Xo, dtype=float).ravel() for Xo in Xs]); log('rank of span{X_j}:', np.linalg.matrix_rank(flat))
    # is M (sqrt p) in span, closed under multiplication?
    def inspan(Y):
        c, *_ = np.linalg.lstsq(flat.T, np.array(Y, dtype=float).ravel(), rcond=None)
        return float(abs(flat.T @ c - np.array(Y, dtype=float).ravel()).max()), float(abs(c - np.rint(c)).max())
    log('sqrt(p)=M^T in Z-span{X_j}: (residual, integrality)', inspan(MoA), ' (1+sqrt p)/2:', inspan((MoA + O(np.eye(n, dtype=np.int64))) // 2))
    cl = max(max(inspan(Xs[i].dot(Xs[j]))) for i in range(8) for j in range(8)); log('Z-span{X_j} closed under multiplication: max defect', cl)
    OUT['order_closed_defect'] = cl
    # trace form discriminant: Tr_{F/Q} trd(x ybar) = Tr_W(X Y^*)/(2r), Y^* = B Y^T B^-1 (row convention)
    Bf = BoA.astype(float); Binv = np.linalg.inv(Bf)
    TF = np.zeros((8, 8))
    for i in range(8):
        for j in range(8):
            TF[i, j] = np.trace(np.array(Xs[i], dtype=float) @ Bf @ np.array(Xs[j], dtype=float).T @ Binv) / (2 * r)
    TFi = np.rint(TF); log('trace form integrality defect', abs(TF - TFi).max())
    dT, _ = bareiss(TFi.astype(np.int64).tolist(), pivot=True)
    log('disc of order (det Tr_{F/Q} trd(x_i xbar_j)) =', dT, ' p^4 =', p ** 4, ' MAXIMAL (disc == p^4):', dT == p ** 4)
    OUT['order_disc'] = str(dT); OUT['order_maximal_by_disc'] = bool(dT == p ** 4)
    # ---- O.w for the violators, direct determinant
    for f, key in [('A/violator_p%d.json' % p, 'x')]:
        try: d = json.load(open(R + f))
        except Exception: continue
        w = O(d[key])
        Wr = O([list(w.dot(Xo)) for Xo in Xs]); G8 = Wr.dot(BoA).dot(Wr.T); d8, _ = bareiss(G8.tolist(), pivot=True)
        b = int(w.dot(BoA).dot(w)); e = int(w.dot(EoA).dot(w)); N = p * b * b - e * e
        form = Fraction(N, 16) ** 4
        log('REFEREE O.w: det_B(O w) =', d8, ' ((p/4)m+m-)^4 = (N/16)^4 =', form, ' equal:', Fraction(d8) == form,
            ' r^8D^4 =', r ** 8 * DA ** 4, ' R_U >=', float(Fraction(r ** 8 * DA ** 4, d8)))
        # stability of the Q-span under Delta: rank of {w X_i X_j}
        allv = np.array([np.array(w.dot(Xs[i]).dot(Xs[j]), dtype=float) for i in range(8) for j in range(8)])
        log('  rank of Delta*(O w) =', np.linalg.matrix_rank(allv), '(8 => Delta-line)')
        gg = 0
        for _t in range(40):
            cols = rng.choice(n, size=8, replace=False); gg = gcd(gg, abs(bareiss(Wr[:, cols].tolist(), pivot=True)[0]))
            if gg == 1: break
        log('  referee O.w rows: gcd of sampled 8x8 minors =', gg, '(1 => O.w = Lambda cap U, R_U exact)')
        if 'Ow_rows' in d:
            Wa = np.array(d['Ow_rows'], dtype=float); Wm = np.array(Wr, dtype=float)
            cc, *_ = np.linalg.lstsq(Wm.T, Wa.T, rcond=None)
            log('  builder Ow_rows in Z-span of referee O.w: residual', float(abs(Wm.T @ cc - Wa.T).max()), 'integrality', float(abs(cc - np.rint(cc)).max()), '|det| of transition', abs(round(float(np.linalg.det(cc)))))
        OUT['referee_Ow_saturated_gcd'] = gg
        OUT['referee_Ow_' + f] = dict(det=str(d8), formula=str(form), equal=bool(Fraction(d8) == form), RU_lower=float(Fraction(r ** 8 * DA ** 4, d8)))
    # same formula on v0 and a few random lattice vectors
    for trial in range(3):
        w = O(rng.integers(-2, 3, size=n))
        Wr = O([list(w.dot(Xo)) for Xo in Xs]); G8 = Wr.dot(BoA).dot(Wr.T); d8, _ = bareiss(G8.tolist(), pivot=True)
        b = int(w.dot(BoA).dot(w)); e = int(w.dot(EoA).dot(w)); N = p * b * b - e * e
        log('  random w: det_B(O w) == (N/16)^4 :', Fraction(d8) == Fraction(N, 16) ** 4, ' Gamma/theta =', N / THR)
except Exception as ex:
    traceback.print_exc()
json.dump(OUT, open(R + 'referee/ref_p%d.json' % p, 'w'), indent=1, default=str)
log('done')
