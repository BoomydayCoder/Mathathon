#!/usr/bin/env python3
"""
hankel_Htheta_interval.py -- rigorous (ball arithmetic, python-flint/arb) verification of

  (H_theta, delta):   2 phi(r) + delta*theta*r <= 2 pi theta r + log(1+theta^2 r^2) + 2 theta r arctan(1/(theta r))
                      for all r >= 1,     phi(r) = (r+1) log(r+1) - (r-1) log(r-1)   (phi(1) = 2 log 2).

Method: [1, 1+1/64] by a crude monotone bound; [1+1/64, 40] by subdivision into balls;
r >= 40: phi(r) <= 2 log(r+1) + 2 and RHS >= 2 pi theta r  (checked at r = 40 + monotonicity of
2 pi theta r - delta theta r - 4 log(r+1) - 4 for r >= 40).
Also checks  16 * exp(-2 pi theta - I(theta)) < 1,  I(theta) = log(1+theta^2) + 2 theta arctan(1/theta).
"""
import sys
from flint import arb, ctx, fmpq
ctx.prec = 80

def check(theta_q, delta_q, nsub=40000, R=40):
    th = arb(theta_q.p) / arb(theta_q.q); de = arb(delta_q.p) / arb(delta_q.q)
    pi = arb.pi()
    def rhs_minus_lhs(rball):
        r = rball
        phi = (r + 1) * (r + 1).log() - (r - 1) * (r - 1).log()
        rhs = 2 * pi * th * r + (1 + th * th * r * r).log() + 2 * th * r * (1 / (th * r)).atan()
        return rhs - 2 * phi - de * th * r
    # piece 1: r in [1, 1+1/64]:  -(r-1)log(r-1) <= (1/64) log 64,  (r+1)log(r+1) <= (2+1/64) log(2+1/64)
    e = arb(1) / 64
    lhs_max = 2 * ((2 + e) * (2 + e).log() + e * arb(64).log()) + de * th * (1 + e)
    rhs_min = 2 * pi * th + (1 + th * th).log() + 2 * th * (1 / (th * (1 + e))).atan()  # each term's lower bound on the piece
    ok1 = bool(rhs_min - lhs_max > 0)
    # piece 2: subdivision of [1+1/64, R]
    a = 1 + 1.0 / 64
    worst = None
    ok2 = True
    for k in range(nsub):
        lo = a + (R - a) * k / nsub; hi = a + (R - a) * (k + 1) / nsub
        ball = arb((lo + hi) / 2, (hi - lo) / 2 * 1.0000001 + 1e-15)
        v = rhs_minus_lhs(ball)
        if not (v > 0):
            ok2 = False
            print("   FAILED on", lo, hi, v)
            break
        low = float((v.mid() - v.rad()).mid()) if hasattr(v, 'mid') else None
        low = float(v.mid()) - float(v.rad())
        if worst is None or low < worst[0]:
            worst = (low, lo)
    # piece 3: r >= R
    r = arb(R)
    g = 2 * pi * th * r - de * th * r - 4 * (r + 1).log() - 4
    gprime_pos = bool(2 * pi * th - de * th - arb(4) / (R + 1) > 0)
    ok3 = bool(g > 0) and gprime_pos
    I = (1 + th * th).log() + 2 * th * (1 / th).atan()
    q = 16 * (-2 * pi * th - I).exp()
    return ok1, ok2, ok3, worst, q

if __name__ == '__main__':
    for (t, d) in [(fmpq(2, 5), fmpq(1, 20)), (fmpq(21, 50), fmpq(2, 5)), (fmpq(9, 20), fmpq(9, 10)), (fmpq(1, 2), fmpq(3, 2))]:
        ok1, ok2, ok3, worst, q = check(t, d)
        print(f"theta={t} delta={d}: [1,1+1/64]: {ok1}; [1+1/64,40]: {ok2} (min certified margin {worst[0]:.5f} near r={worst[1]:.3f}); r>=40: {ok3};  16 exp(-2 pi theta - I(theta)) = {float(q.mid()):.4f} < 1: {bool(q < 1)}")
