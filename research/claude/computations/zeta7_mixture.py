#!/usr/bin/env python3
"""
zeta7_mixture.py -- direction 1 (mixtures): the positive weight  w_7(y) - y^2 w_5(y) > 0
(positive because (2 pi)^2 > 30) gives the functional  Lambda(R) = mu7_X7(R) - mu5_X5(t R),  whose value
at 1/(t+j^2) is  j^6 [ (X7+X5) - (H_j^(7)+H_j^(5)) ] + rational.   So all determinants are polynomials
in the ONE variable Z = zeta(7)+zeta(5)  (and w_7 + y^2 w_5 gives Z = zeta(7)-zeta(5)).
Exact computation of cost/decay/total for det[Lambda(W t^{i+k}/D_S)], evaluated at Z.
"""
import sys
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from zeta7_lib import *
from flint import arb, ctx


def run(K, N, a, h, sign=-1):
    nodes = list(range(N + 1, K + 1))
    W = D_int(range(1, N + 1)) ** a
    nm = 2 * h - 1
    a7, b7 = moments(7, nodes, W, nm)
    a5, b5 = moments(5, nodes, W * fmpq_poly([0, 1]), nm)
    assert all(b7[e] == -b5[e] for e in range(nm))
    # sign=-1: weight w7 - y^2 w5, variable Z = X7+X5 ; sign=+1: weight w7 + y^2 w5, Z = X7-X5
    A = fmpq_mat(h, h, [a7[i + k] + sign * a5[i + k] for i in range(h) for k in range(h)])
    B = fmpq_mat(h, h, [b7[i + k] for i in range(h) for k in range(h)])
    cs = det_poly_affine(A, B)
    deg = max(k for k, c in enumerate(cs) if c != 0)
    cost = -log_abs_fmpq(content(cs))
    prec = 4000 + 80 * h * h
    old = ctx.prec
    ctx.prec = prec
    z = arb(7).zeta() - sign * arb(5).zeta()
    acc = arb(0)
    for c in reversed(cs):
        acc = acc * z + arb(c.p) / arb(c.q)
    ok = not acc.contains(0)
    ld = float(abs(acc).log().mid())
    sg = 1 if acc > 0 else -1
    ctx.prec = old
    name = "zeta(7)+zeta(5)" if sign < 0 else "zeta(7)-zeta(5)"
    print(f"Z={name}: K={K} N={N} W=D_N^{a} h={h} deg={deg} sign={sg:+d} ok={ok} | cost/K^2={cost / K ** 2:.4f} logdet/K^2={ld / K ** 2:.4f} total/K^2={(cost + ld) / K ** 2:+.4f}")
    sys.stdout.flush()


if __name__ == '__main__':
    for K in (20, 40, 60):
        run(K, 0, 0, K, -1)
        run(K, 0, 0, K, +1)
    run(40, 3, 5, 37, -1)
    run(80, 6, 5, 74, -1)
