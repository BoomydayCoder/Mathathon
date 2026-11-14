#!/usr/bin/env python3
"""
catalan_lib.py -- exact tools for Hankel determinants of "Hermite / Abel-Plana functionals" for
Hurwitz-type constants.  Generalises hankel_lib.py (which is the case kernel='bose', alpha=1, s odd).

Two kernels f, any integer s >= 1 (s >= 2 for 'bose'), any rational shift alpha in (0,1]:

  kernel 'bose':  f(y) = 1/(e^{2 pi y}-1),  weight  w(y) = 2 (-1)^(s-1) y^s f^{(s-1)}(y)/(s-1)!  > 0
      mu(t^e)             = (-1)^e B_{2e+2} (2e+s)! / ((s-1)! (2e+2)!)
      mu_X(1/(t+a^2))     = a^(s-1) (X - sum_{v<j} (v+alpha)^(-s)) - 1/(2a) - 1/(s-1),   a = alpha+j,
      X0 = zeta(s, alpha)  (Hurwitz).          [Hermite's formula]

  kernel 'sinh':  f(y) = 1/sinh(pi y),      weight  w(y) = (-1)^(s-1) y^s f^{(s-1)}(y)/(s-1)!  > 0
      mu(t^e)             = (2^(2e+2)-1) |B_{2e+2}| (2e+s)! / ((s-1)! (2e+2)!)
      mu_X(1/(t+a^2))     = (-1)^j a^(s-1) (X - sum_{v<j} (-1)^v (v+alpha)^(-s)) - 1/(2a),   a = alpha+j,
      X0 = eta(s, alpha) = sum_{n>=0} (-1)^n (n+alpha)^(-s).   [Abel-Plana for alternating series]
      eta(s,1/2) = 2^s beta(s)  (s=2: 4*Catalan; s=1: pi/2);  eta(s,1) = (1-2^(1-s)) zeta(s)  (s=1: log 2).

  kernel 'cosh' (s=2, alpha=1 only; "functional II" for Catalan):
      weight w(y) = (pi y/4) sinh(pi y)/cosh(pi y)^2 > 0   [= -(y/2) d/dy (1/(2 cosh(pi y)))]
      mu(t^e)             = (2e+1) |E_{2e}| / 2^(2e+3)            (Euler numbers)
      mu_X(1/(t+n^2))     = (-1)^n (X - sum_{k<2n} chi_4(k)/k^2),  n = 1,2,...;   X0 = G (Catalan).
      [Abel-Plana for sum chi_4(n) g(n) = int_0^inf Re g(iy) / (2 cosh(pi y/2)) dy, one integration by parts]

In all cases mu_{X0}(R) = int_0^infty R(y^2) w(y) dy for R with simple poles at t = -a^2 (a in alpha+N).
(Derivation in the report; numerical verification in catalan_functional_check.py.)

Hankel family:  nodes a_j = alpha + j, j = j0..j0+K-1,  D(t) = prod (t + a_j^2),
    G(X) = [ mu_X( t^{i+i'} / D(t) ) ]_{0<=i,i'<h},   Delta = det G(X) in Q[X].

Everything is exact rational arithmetic (python-flint).  This file only computes; it proves nothing.
"""
import math, sys
from flint import fmpq, fmpz, fmpq_mat, fmpq_poly, arb, ctx

sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import primes_upto, vp, fact, det_poly_affine, content, gauss_val, log_abs_fmpq, log_abs_arb


def mu_monomial(e, s, kernel):
    B = fmpq.bernoulli(2 * e + 2)
    c = fmpq(fact(2 * e + s), fact(s - 1) * fact(2 * e + 2))
    if kernel == 'bose':
        return fmpq((-1) ** e) * B * c
    if kernel == 'sinh':
        return fmpq(2 ** (2 * e + 2) - 1) * abs(B) * c
    if kernel == 'cosh':
        assert s == 2
        return fmpq(2 * e + 1) * fmpq(abs(int(fmpz.euler_number(2 * e))), 2 ** (2 * e + 3))
    raise ValueError(kernel)


def pole_values(s, kernel, alpha, J):
    """list of (const, coefX) with mu_X(1/(t+a_j^2)) = const + coefX*X, a_j = alpha+j, j=0..J-1."""
    alpha = fmpq(alpha)
    out = []
    S = fmpq(0)  # partial sum
    for j in range(J):
        a = alpha + j
        if kernel == 'cosh':
            assert s == 2 and alpha == 1
            n = j + 1          # node t = -n^2;  S = sum_{k<2n} chi_4(k)/k^2
            Sd = sum((fmpq((1 if k % 4 == 1 else -1), k * k) for k in range(1, 2 * n, 2)), fmpq(0))
            cX = fmpq((-1) ** n)
            c0 = -cX * Sd
        elif kernel == 'bose':
            cX = a ** (s - 1)
            c0 = -cX * S - 1 / (2 * a) - fmpq(1, s - 1)
            S += 1 / a ** s
        else:
            sg = (-1) ** j
            cX = sg * a ** (s - 1)
            c0 = -cX * S - 1 / (2 * a)
            S += sg / a ** s
        out.append((c0, cX))
    return out


def moments(s, kernel, alpha, K, h=None, j0=0, tscale=1):
    """(alpha_e, beta_e), e=0..2h-2: mu_X(t^e/D) = alpha_e + beta_e X, D = prod_{j=j0}^{j0+K-1}(t+a_j^2).
    tscale: use the variable T = tscale * t (nodes -tscale*a_j^2) -- changes Delta only by tscale^(O(K))."""
    if h is None:
        h = K
    alpha = fmpq(alpha)
    nodes = [-(alpha + j) ** 2 * tscale for j in range(j0, j0 + K)]
    pv = pole_values(s, kernel, alpha, j0 + K)[j0:]
    # mu_X(1/(T - T_j)) with T = tscale*t:  1/(tscale*(t+a^2)) -> pv/tscale
    D = fmpq_poly([1])
    for tj in nodes:
        D = D * fmpq_poly([-tj, 1])
    Dp = D.derivative()
    res = [1 / Dp(tj) for tj in nodes]
    nmom = 2 * h - 1
    maxdeg = max(0, nmom - 1 - K)
    mm = [mu_monomial(e, s, kernel) * fmpq(tscale) ** e for e in range(maxdeg + 1)]
    al, be = [], []
    Tpow = fmpq_poly([1])
    X = fmpq_poly([0, 1])
    pw = [fmpq(1)] * K
    for e in range(nmom):
        a0 = fmpq(0)
        b0 = fmpq(0)
        if e >= K:
            quo = Tpow // D
            cs = quo.coeffs()
            a0 += sum((cs[k] * mm[k] for k in range(len(cs))), fmpq(0))
        for j in range(K):
            r = res[j] * pw[j] / tscale
            a0 += r * pv[j][0]
            b0 += r * pv[j][1]
            pw[j] = pw[j] * nodes[j]
        al.append(a0)
        be.append(b0)
        Tpow = Tpow * X
    return al, be


def hankel(mom, h):
    return fmpq_mat(h, h, [mom[i + j] for i in range(h) for j in range(h)])


def X0_arb(s, kernel, alpha, prec):
    """the constant X0 as an arb at binary precision prec."""
    old = ctx.prec
    ctx.prec = prec
    alpha = fmpq(alpha)
    try:
        from flint import acb
        a = acb(arb(alpha.p) / arb(alpha.q))
        if kernel == 'cosh':
            val = arb.const_catalan()
        elif kernel == 'bose':
            val = acb(s).zeta(a).real
        else:
            # eta(s,alpha) = 2^-s (zeta(s,alpha/2) - zeta(s,(alpha+1)/2)); for s=1 use digamma
            if s == 1:
                val = ((a + 1) / 2).digamma().real / 2 - (a / 2).digamma().real / 2
            else:
                val = (acb(s).zeta(a / 2).real - acb(s).zeta((a + 1) / 2).real) / arb(2) ** s
    finally:
        ctx.prec = old
    return val


def eval_poly(cs, x):
    acc = arb(0)
    for c in reversed(cs):
        acc = acc * x + arb(c.p) / arb(c.q)
    return acc


def analyse(s, kernel, alpha, K, h=None, j0=0, tscale=1, want_vals=False):
    if h is None:
        h = K
    al, be = moments(s, kernel, alpha, K, h, j0, tscale)
    A, B = hankel(al, h), hankel(be, h)
    cs = det_poly_affine(A, B)
    deg = max(k for k, c in enumerate(cs) if c != 0)
    cont = content(cs)
    prec = 4000 + 80 * K * K
    while True:
        old = ctx.prec
        ctx.prec = prec
        x0 = X0_arb(s, kernel, alpha, prec)
        val = eval_poly(cs, x0)
        ctx.prec = old
        if (not val.contains(0)) and val.rad() < abs(val.mid()) * 1e-20:
            break
        prec *= 2
    logdet, _ = log_abs_arb(val)
    sign = 1 if val > 0 else -1
    logcont = log_abs_fmpq(cont)
    res = dict(s=s, kernel=kernel, alpha=alpha, K=K, h=h, deg=deg, cs=cs, content=cont, sign=sign,
               logdet=logdet, cost=-logcont, total=logdet - logcont, lead=cs[deg])
    if want_vals:
        vals = {}
        rest = cont
        for p in primes_upto(8 * K + 50):
            v = gauss_val(cs, p)
            if v != 0:
                vals[p] = v
                rest = rest / fmpq(p) ** v
        res['vals'] = vals
        res['content_rest'] = rest
    return res


if __name__ == '__main__':
    # self-test against hankel_lib (bose, alpha=1, s=3)
    import hankel_lib as HL
    al, be = HL.moments_F1(3, 8, 0, 0, 8)
    cs_old = HL.det_poly_affine(HL.hankel(al, 8), HL.hankel(be, 8))
    r = analyse(3, 'bose', 1, 8)
    print('agreement with hankel_lib (s=3,K=8):', all(a == b for a, b in zip(cs_old, r['cs'])))
    print('zeta(3): total/K^2 =', r['total'] / 64, 'sign', r['sign'], 'deg', r['deg'])
    r = analyse(2, 'sinh', fmpq(1, 2), 8, want_vals=True)
    print('Catalan (X=4G): deg', r['deg'], 'sign', r['sign'], 'cost/K^2', r['cost'] / 64, 'logdet/K^2', r['logdet'] / 64,
          'total/K^2', r['total'] / 64, r['vals'], r['content_rest'])
