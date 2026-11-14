#!/usr/bin/env python3
"""
catalan2_lib.py -- exact tools for "mixed-node" Hankel determinants for Catalan's constant G.

Two positive weights on (0,inf) (both verified numerically to 30 digits in catalan2_checks.py):

  'bose':  w_B(y) = -2 y^2 d/dy [1/(e^{2 pi y}-1)]            mu(t^e) = |B_{2e+2}|
           mu(1/(t+a^2)) = a zeta(2,a) - 1/(2a) - 1           (Hermite, any real a>0)
           zeta(2,1/4) = pi^2+8G, zeta(2,3/4) = pi^2-8G, zeta(2,1/2) = pi^2/2, zeta(2,1) = pi^2/6
  'sinh':  w_I(y) = pi y^2 cosh(pi y)/sinh(pi y)^2            mu(t^e) = (2^{2e+2}-1)|B_{2e+2}|
           mu(1/(t+a^2)) = a eta(2,a) - 1/(2a),  eta(2,a) = sum_{m>=0} (-1)^m/(m+a)^2
           eta(2,1/2) = 4G, eta(2,1) = pi^2/12
  (mu(R) = int_0^inf R(y^2) w(y) dy.)   So for nodes a in (1/4)Z_{>0} ('bose') or (1/2)Z_{>0} ('sinh')
  mu(1/(t+a^2)) = c0 + cG*G + cP*pi^2 with explicit rationals; returned by pole_value().

A test function is  Psi(t) = sum_a r_a/(t+a^2) + P(t)  (dict of residues + polynomial part); moments
mu(t^e Psi) are triples (c0,cG,cP).  Hankel matrices [mu(t^{i+j} Psi)]_{i,j<h}; when all cP vanish the
determinant is a polynomial in X=G:  Delta(X) = det(A + X B).

Everything is exact rational arithmetic (python-flint).  This file only computes; it proves nothing.
"""
import math, sys
from flint import fmpq, fmpz, fmpq_mat, fmpq_poly, arb, ctx
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import primes_upto, vp, det_poly_affine, content, gauss_val, log_abs_fmpq, log_abs_arb

_PS = {}


def _partial(kernel, alpha, j):
    """sum_{v<j} eps^v (alpha+v)^-2, eps = 1 (bose) or -1 (sinh); cached."""
    key = (kernel, alpha)
    L = _PS.setdefault(key, [fmpq(0)])
    while len(L) <= j:
        v = len(L) - 1
        sg = 1 if kernel == 'bose' or v % 2 == 0 else -1
        L.append(L[-1] + fmpq(sg) / (alpha + v) ** 2)
    return L[j]


def pole_value(kernel, a):
    """(c0, cG, cP):  mu(1/(t+a^2)) = c0 + cG*G + cP*pi^2."""
    a = fmpq(a)
    fl = a.p // a.q
    alpha = a - fl
    if alpha == 0:
        alpha = fmpq(1)
        j = fl - 1
    else:
        j = fl
    S = _partial(kernel, alpha, j)
    if kernel == 'bose':
        zG, zP = {fmpq(1, 4): (8, fmpq(1)), fmpq(3, 4): (-8, fmpq(1)), fmpq(1, 2): (0, fmpq(1, 2)), fmpq(1): (0, fmpq(1, 6))}[alpha]
        return (-a * S - 1 / (2 * a) - 1, a * zG, a * zP)
    if kernel == 'sinh':
        zG, zP = {fmpq(1, 2): (4, fmpq(0)), fmpq(1): (0, fmpq(1, 12))}[alpha]
        sg = (-1) ** j
        return (-sg * a * S - 1 / (2 * a), sg * a * zG, sg * a * zP)
    raise ValueError(kernel)


def mu_mono(kernel, e):
    B = abs(fmpq.bernoulli(2 * e + 2))
    return B if kernel == 'bose' else (2 ** (2 * e + 2) - 1) * B


def moments(kernel, res, P, E):
    """triples for mu(t^e Psi), e = 0..E-1;  Psi = sum_a res[a]/(t+a^2) + P(t)  (P list of fmpq coefficients)."""
    P = list(P) if P is not None else []
    nodes = sorted(res)
    pv = {a: pole_value(kernel, a) for a in nodes}
    mm = [mu_mono(kernel, k) for k in range(E + len(P) + 1)]
    x = {a: -a * a for a in nodes}
    cur = {a: res[a] for a in nodes}          # r_a (-a^2)^e
    sig = []                                   # sigma_m = sum_a r_a (-a^2)^m
    out = []
    for e in range(E):
        c0 = sum((cur[a] * pv[a][0] for a in nodes), fmpq(0))
        cG = sum((cur[a] * pv[a][1] for a in nodes), fmpq(0))
        cP = sum((cur[a] * pv[a][2] for a in nodes), fmpq(0))
        c0 += sum((sig[e - 1 - k] * mm[k] for k in range(e)), fmpq(0))
        c0 += sum((P[i] * mm[i + e] for i in range(len(P))), fmpq(0))
        out.append((c0, cG, cP))
        sig.append(sum((cur[a] for a in nodes), fmpq(0)))
        for a in nodes:
            cur[a] = cur[a] * x[a]
    return out


def dprime(nodes):
    """{a: D'(-a^2)},  D(t) = prod_b (t+b^2)."""
    out = {}
    for a in nodes:
        pr = fmpq(1)
        for b in nodes:
            if b != a:
                pr *= (b * b - a * a)
        out[a] = pr
    return out


def psi_from_fraction(Ncoef, nodes):
    """Psi = N/D, D = prod (t+a^2): returns (residues, polynomial part coefficients)."""
    N = fmpq_poly(Ncoef)
    D = fmpq_poly([1])
    for a in nodes:
        D = D * fmpq_poly([a * a, 1])
    dp = dprime(nodes)
    res = {a: N(-a * a) / dp[a] for a in nodes}
    quo = N // D if N.degree() >= D.degree() else fmpq_poly([0])
    return res, list(quo.coeffs())


def eval_poly(cs, x):
    acc = arb(0)
    for c in reversed(cs):
        acc = acc * x + arb(c.p) / arb(c.q)
    return acc


def analyse_mom(mom, h, pmax, label='', want_vals=True, shift=0):
    """Hankel [mom[i+j+shift]] (entries must have cP = 0); exact det polynomial in X=G; cost/decay/total."""
    assert all(m[2] == 0 for m in mom[shift:shift + 2 * h - 1]), 'pi^2 part not cancelled'
    A = fmpq_mat(h, h, [mom[i + j + shift][0] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [mom[i + j + shift][1] for i in range(h) for j in range(h)])
    cs = det_poly_affine(A, B)
    nz = [k for k, c in enumerate(cs) if c != 0]
    if not nz:
        return dict(label=label, h=h, deg=-1)
    deg = max(nz)
    cont = content(cs)
    prec = 3000 + 40 * h * h
    while True:
        old = ctx.prec
        ctx.prec = prec
        val = eval_poly(cs, arb.const_catalan())
        ctx.prec = old
        if (not val.contains(0)) and val.rad() < abs(val.mid()) * 1e-20:
            break
        prec *= 2
        if prec > 4e6:
            return dict(label=label, h=h, deg=deg, zero=True)
    logdet, _ = log_abs_arb(val)
    sign = 1 if val > 0 else -1
    logcont = log_abs_fmpq(cont)
    r = dict(label=label, h=h, deg=deg, sign=sign, logdet=logdet, cost=-logcont, total=logdet - logcont,
             loglead=log_abs_fmpq(cs[deg]), cs=cs)
    if want_vals:
        # valuations of the monic polynomial's "denominator": content vs leading coefficient
        vals, valslc = {}, {}
        rest = cont
        for p in primes_upto(pmax):
            v = gauss_val(cs, p)
            vals[p] = v
            valslc[p] = vp(cs[deg], p) - v      # v_p(lc of primitive polynomial) >= 0
            rest = rest / fmpq(p) ** v
        r['vals'] = vals
        r['vlc'] = valslc
        r['rest'] = rest
    return r


def fmt(r, scale):
    if r.get('deg', -1) < 0 or r.get('zero'):
        return f"{r['label']:44s} h={r['h']:3d}  DEGENERATE"
    v2 = r['vlc'].get(2, 0) if 'vlc' in r else 0
    gam = (r['loglead'] + r['cost']) / scale      # log lc(P*)/scale
    return (f"{r['label']:44s} h={r['h']:3d} deg={r['deg']:3d} sign={r['sign']:+d} gamma={gam:7.4f} "
            f"delta={(r['loglead'] - r['logdet']) / scale:7.4f} total={r['total'] / scale:+8.4f} | v2(lc)={v2:5d} "
            f"({v2 * math.log(2) / scale:.4f})")
