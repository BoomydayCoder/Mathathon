#!/usr/bin/env python3
"""
zeta7_lib.py -- general exact engine for determinants of the simple-pole positive functional of zeta(s).

Functional (odd s):  mu_X(R) = "int_0^oo R(y^2) w_s(y) dy" with zeta(s) replaced by X,
    w_s(y) = 2 y^s f^{(s-1)}(y)/(s-1)!,  f = 1/(e^{2 pi y}-1)  (positive),
    mu(t^e)          = (-1)^e B_{2e+2} (2e+s)!/((s-1)!(2e+2)!),
    mu_X(1/(t+a^2))  = M_s(a) = a^{s-1} zeta(s,a) - 1/(2a) - 1/(s-1)      (Hermite's formula),
with  zeta(s,j) = X - H_{j-1}^{(s)}  (a = j integer)  and
      zeta(s,j-1/2) = (2^s-1) X - 2^s sum_{v<j} (2v-1)^{-s}  (half-integers).
A dilated weight w_s(b y) with poles at -j^2 is the same thing as w_s with poles at -(b j)^2, so
"faster/slower kernels" e^{-2 pi b y} are covered by the choice of the node set.

General determinant:  det [ mu_X( W(t) t^{lam_i + lam'_k} / D_S(t) ) ],  D_S = prod_{a in S}(t+a^2),
W a polynomial (W>0 on t>0 for positivity), lam, lam' increasing exponent lists (minors of the moment
matrix; Hankel for lam=lam'=0..h-1).  All arithmetic exact (python-flint).  This file only computes.
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from flint import fmpq, fmpz, fmpq_mat, fmpq_poly, arb, ctx
from hankel_lib import (mu_monomial, det_poly_affine, content, log_abs_fmpq, eval_poly_at_zeta,
                        log_abs_arb, primes_upto, gauss_val, vp)


def node_value(a, s, cache={}):
    """(alpha, beta) with M_s(a) = alpha + beta*X, a = fmpq integer or half-integer > 0."""
    a = fmpq(a)
    key = (a.p, a.q, s)
    if key in cache:
        return cache[key]
    if a.q == 1:
        j = int(a.p)
        Hjm1 = sum((fmpq(1, v ** s) for v in range(1, j)), fmpq(0))
        zc, zb = -Hjm1, fmpq(1)
    elif a.q == 2:
        j = (int(a.p) + 1) // 2          # a = j - 1/2
        O = sum((fmpq(1, (2 * v - 1) ** s) for v in range(1, j)), fmpq(0))
        zc, zb = -fmpq(2 ** s) * O, fmpq(2 ** s - 1)
    else:
        raise ValueError
    pw = a ** (s - 1)
    al = pw * zc - 1 / (2 * a) - fmpq(1, s - 1)
    be = pw * zb
    cache[key] = (al, be)
    return al, be


def D_of(nodes):
    P = fmpq_poly([1])
    for a in nodes:
        a = fmpq(a)
        P = P * fmpq_poly([a * a, 1])
    return P


def moments(s, nodes, W, nmom):
    """alpha_e, beta_e (e<nmom) with mu_X(W t^e / D_S) = alpha_e + beta_e X."""
    nodes = [fmpq(a) for a in nodes]
    DS = D_of(nodes)
    W = fmpq_poly(W)
    DSp = DS.derivative()
    rho = []
    for a in nodes:
        t0 = -a * a
        wj = W(t0)
        if wj != 0:
            al, be = node_value(a, s)
            rho.append((t0, wj / DSp(t0), al, be))
    maxdeg = max(0, W.degree() + nmom - 1 - DS.degree())
    mu_mon = [mu_monomial(e, s) for e in range(maxdeg + 1)]
    alpha, beta = [], []
    Wt = W
    tpoly = fmpq_poly([0, 1])
    tp = [fmpq(1)] * len(rho)
    for e in range(nmom):
        if Wt.degree() >= DS.degree():
            quo, _ = divmod(Wt, DS)
            cs = quo.coeffs()
            al = sum((cs[k] * mu_mon[k] for k in range(len(cs))), fmpq(0))
        else:
            al = fmpq(0)
        be = fmpq(0)
        for idx, (t0, r, a_, b_) in enumerate(rho):
            rj = r * tp[idx]
            al += rj * a_
            be += rj * b_
            tp[idx] = tp[idx] * t0
        alpha.append(al)
        beta.append(be)
        Wt = Wt * tpoly
    return alpha, beta


def analyse(s, nodes, W=None, lam=None, lamp=None, h=None, pmax=None):
    """det [mu_X(W t^{lam_i+lam'_k}/D_S)]; returns dict(deg, cost, logdet, total, sign, vals)."""
    if W is None:
        W = fmpq_poly([1])
    if lam is None:
        lam = list(range(h))
    if lamp is None:
        lamp = lam
    h = len(lam)
    nmom = lam[-1] + lamp[-1] + 1
    al, be = moments(s, nodes, W, nmom)
    A = fmpq_mat(h, h, [al[i + k] for i in lam for k in lamp])
    B = fmpq_mat(h, h, [be[i + k] for i in lam for k in lamp])
    cs = det_poly_affine(A, B)
    nz = [k for k, c in enumerate(cs) if c != 0]
    if not nz:
        return dict(deg=-1)
    deg = max(nz)
    cont = content(cs)
    prec = 3000 + 80 * h * h
    while True:
        val = eval_poly_at_zeta(cs, s, prec)
        if not val.contains(0) and val.rad() < abs(val.mid()) * 1e-20:
            break
        prec *= 2
        if prec > 4e6:
            return dict(deg=deg, cost=-log_abs_fmpq(cont), logdet=float('nan'), total=float('nan'), sign=0, cs=cs)
    logdet, _ = log_abs_arb(val)
    cost = -log_abs_fmpq(cont)
    res = dict(deg=deg, h=h, cost=cost, logdet=logdet, total=cost + logdet, sign=(1 if val > 0 else -1), cs=cs)
    if pmax:
        res['vals'] = {p: gauss_val(cs, p) for p in primes_upto(pmax) if gauss_val(cs, p) != 0}
    return res


def D_int(js):
    return D_of([fmpq(j) for j in js])


def fmt(r, K):
    if r.get('deg', -1) < 0:
        return "det == 0"
    K2 = K * K
    return (f"h={r['h']:3d} deg={r['deg']:3d} sign={r['sign']:+d} | cost/K^2={r['cost'] / K2:8.4f} "
            f"logdet/K^2={r['logdet'] / K2:8.4f} total/K^2={r['total'] / K2:+8.4f} | total/deg^2={r['total'] / max(1, r['deg']) ** 2:+8.4f}")
