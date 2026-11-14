#!/usr/bin/env python3
"""
hankel_lib.py -- exact tools for Fauzan-type Hankel determinants for odd zeta values.

For odd s >= 3 the Q-linear functional mu_X (values in Q[X], affine in X) is defined on
rational functions of t with at most simple poles at t = -j^2 (j = 1,2,...) by

    mu_X(t^e)        = (-1)^e B_{2e+2} (2e+s)! / ((s-1)! (2e+2)!)
    mu_X(1/(t+j^2))  = j^(s-1) (X - H_j^(s)) - 1/(s-1) + 1/(2j).

At X = zeta(s) it is integration of R(y^2) against the positive weight
w_s(y) = 2 y^s f^{(s-1)}(y)/(s-1)!,  f(y) = 1/(e^{2 pi y}-1).

Family F1(s; K, N, a, h):   G(X) = [ mu_X( D_N(t)^a t^{i+j} / D_K(t) ) ]_{0<=i,j<h},
D_m(t) = prod_{j<=m} (t + j^2).   (Fauzan: s=5, a=6, K=40n, N=3n, h=37n.)

Everything here is exact rational arithmetic (python-flint).  No claims are made by this
file; it only computes.
"""
import math
from flint import fmpq, fmpz, fmpq_mat, fmpz_mat, fmpz_poly, fmpq_poly, arb, ctx


def primes_upto(n):
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


def vp_int(n, p):
    n = int(n)
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def vp(q, p):
    """p-adic valuation of an fmpq."""
    q = fmpq(q)
    if q == 0:
        return 10 ** 9
    return vp_int(q.p, p) - vp_int(q.q, p)


_fact = [1]


def fact(n):
    while len(_fact) <= n:
        _fact.append(_fact[-1] * len(_fact))
    return _fact[n]


def mu_monomial(e, s):
    """mu(t^e) for the weight of zeta(s)."""
    return fmpq((-1) ** e) * fmpq.bernoulli(2 * e + 2) * fmpq(fact(2 * e + s), fact(s - 1) * fact(2 * e + 2))


def harmonic_s(K, s):
    H = [fmpq(0)]
    for j in range(1, K + 1):
        H.append(H[-1] + fmpq(1, j ** s))
    return H


def D_poly(m):
    """D_m(t) = prod_{j<=m} (t+j^2) as fmpz_poly."""
    P = fmpz_poly([1])
    for j in range(1, m + 1):
        P = P * fmpz_poly([j * j, 1])
    return P


def moments_F1(s, K, N, a, h, extra_num=None):
    """Return lists (alpha_e, beta_e), e = 0..2h-2, with
       mu_X( W t^e / D_K ) = alpha_e + beta_e X,  W = D_N^a * extra_num.
    """
    DK = D_poly(K)
    W = D_poly(N) ** a
    if extra_num is not None:
        W = W * extra_num
    H = harmonic_s(K, s)
    # residues of W/D_K at t=-j^2
    DKp = DK.derivative()
    rho = {}
    for j in range(1, K + 1):
        wj = W(fmpz(-j * j))
        if wj != 0:
            rho[j] = fmpq(wj, DKp(fmpz(-j * j)))
    nmom = 2 * h - 1
    alpha, beta = [], []
    # polynomial parts
    maxdeg = W.degree() + nmom - 1 - K
    mu_mon = [mu_monomial(e, s) for e in range(max(0, maxdeg) + 1)]
    Wt = W
    for e in range(nmom):
        if Wt.degree() >= K:
            quo, rem = divmod(Wt, DK)
            cs = quo.coeffs()
            polyval = sum((fmpq(int(cs[k])) * mu_mon[k] for k in range(len(cs))), fmpq(0))
        else:
            polyval = fmpq(0)
        al = polyval
        be = fmpq(0)
        for j, r in rho.items():
            rj = r * fmpq((-j * j) ** e)
            be += rj * (j ** (s - 1))
            al += rj * (-(j ** (s - 1)) * H[j] - fmpq(1, s - 1) + fmpq(1, 2 * j))
        alpha.append(al)
        beta.append(be)
        Wt = Wt * fmpz_poly([0, 1])
    return alpha, beta


def hankel(mom, h):
    return fmpq_mat(h, h, [mom[i + j] for i in range(h) for j in range(h)])


def det_poly_affine(A, B):
    """det(A + X B) as a list of fmpq coefficients c_0..c_h (exact)."""
    h = A.nrows()
    if B.rank() == h:
        # det(A+XB) = det(B) det(X + B^{-1}A) = det(B) * charpoly(-B^{-1}A)(X)
        Binv = B.inv()
        M = Binv * A
        negM = fmpq_mat(h, h, [-x for x in M.entries()])
        cp = negM.charpoly()
        dB = B.det()
        cs = [dB * c for c in cp.coeffs()]
        cs += [fmpq(0)] * (h + 1 - len(cs))
        return cs
    # fallback: evaluation + Newton interpolation
    xs = list(range(h + 1))
    ys = []
    for x in xs:
        M = fmpq_mat(h, h, [A[i, j] + x * B[i, j] for i in range(h) for j in range(h)])
        ys.append(M.det())
    # Newton divided differences
    n = len(xs)
    dd = ys[:]
    for k in range(1, n):
        for i in range(n - 1, k - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) / fmpq(xs[i] - xs[i - k])
    poly = fmpq_poly([0])
    for k in range(n - 1, -1, -1):
        poly = poly * fmpq_poly([-xs[k], 1]) + fmpq_poly([dd[k]])
    cs = list(poly.coeffs())
    cs += [fmpq(0)] * (h + 1 - len(cs))
    return cs


def content(cs):
    """content c>0 of a rational polynomial: cs/c is a primitive integer polynomial."""
    nz = [c for c in cs if c != 0]
    g = fmpz(0)
    l = fmpz(1)
    for c in nz:
        g = g.gcd(c.p)
        l = l.lcm(c.q)
    return fmpq(g, l)


def gauss_val(cs, p):
    return min(vp(c, p) for c in cs if c != 0)


def log_abs_fmpq(q, prec=200):
    old = ctx.prec
    ctx.prec = prec
    q = fmpq(q)
    val = float(arb(abs(q.p)).log() - arb(q.q).log())
    ctx.prec = old
    return val


def eval_poly_at_zeta(cs, s, prec):
    """evaluate sum c_k zeta(s)^k with arb at the given binary precision; returns arb."""
    old = ctx.prec
    ctx.prec = prec
    z = arb(s).zeta()
    acc = arb(0)
    for c in reversed(cs):
        acc = acc * z + arb(c.p) / arb(c.q)
    ctx.prec = old
    return acc


def log_abs_arb(x, prec=200):
    old = ctx.prec
    ctx.prec = max(prec, 200)
    v = abs(x).log()
    ctx.prec = old
    return float(v.mid()), float(v.rad())


def analyse_F1(s, K, N, a, h, pmax=None, prec=None, verbose=False):
    """Compute det G(X) for F1(s;K,N,a,h); return dict with
       coefficients, content, log|det(zeta(s))|, valuations of content at primes <= pmax.
    """
    alpha, beta = moments_F1(s, K, N, a, h)
    A = hankel(alpha, h)
    B = hankel(beta, h)
    cs = det_poly_affine(A, B)
    deg = max(k for k, c in enumerate(cs) if c != 0)
    cont = content(cs)
    if prec is None:
        prec = 4000 + 60 * K * K
    while True:
        val = eval_poly_at_zeta(cs, s, prec)
        if val.rad() < abs(val.mid()) * 1e-20 and not val.contains(0):
            break
        prec *= 2
    logdet, _ = log_abs_arb(val)
    sign = 1 if val > 0 else -1
    logcont = log_abs_fmpq(cont)
    res = dict(s=s, K=K, N=N, a=a, h=h, deg=deg, cs=cs, content=cont,
               logdet=logdet, sign=sign, logcontent=logcont,
               logprim=logdet - logcont)
    if pmax is None:
        pmax = 8 * K + 2 * a * N + 4 * h
    vals = {}
    for p in primes_upto(pmax):
        v = gauss_val(cs, p)
        if v != 0:
            vals[p] = v
    res['vals'] = vals
    # is content fully explained by primes <= pmax ?
    rest = cont
    for p, v in vals.items():
        rest = rest / (fmpq(p) ** v)
    res['content_rest'] = rest
    return res


if __name__ == '__main__':
    import sys
    # tiny self-test: positivity and the smallest instance
    r = analyse_F1(5, 6, 0, 0, 5)
    print('F1(s=5,K=6,N=0,a=0,h=5): deg', r['deg'], 'log|prim(zeta)|', r['logprim'], 'sign', r['sign'])
    for s in (3, 5, 7):
        r = analyse_F1(s, 6, 0, 0, 5)
        print(s, r['logprim'], r['logprim'] / r['deg'], r['vals'], r['content_rest'])
