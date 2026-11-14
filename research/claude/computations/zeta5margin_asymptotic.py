#!/usr/bin/env python3
"""
zeta5margin_asymptotic.py -- exact asymptotic constants for the zeta(5) normalisation when the outer-range
bound (4.14) is replaced by  max( (4.14), water filling )  [Theorem W of the report].

Continuum water filling (limit of gamma_p^wf / p, x = K/p), family W = D_N^a, weight of zeta(s):
  classes z in (0,1/2):  l_K = [x-z]+[x+z]+1,  l_N = [alpha x - z]+[alpha x + z]+1,
  doubled base B(z) = a l_N - l_K - (s-1); doubled weights B + 2i (i>=0); take the lambda*x smallest (measure).
  Gamma_wf(x) = sum of the selected doubled weights.
All arithmetic is exact (fractions).  Piecewise linearity is established from an explicit superset of
breakpoints and re-checked by exact collinearity tests on every cell.
Nothing is claimed by this file; it only computes.
"""
import sys, math
from fractions import Fraction as Fr
from math import floor

def frac(u): return u - floor(u)

class Model:
    def __init__(self, s, a, alpha, lam):
        self.s, self.a, self.alpha, self.lam = s, a, Fr(alpha), Fr(lam)

    def types(self, x):
        al = self.alpha
        cuts = {Fr(0), Fr(1, 2)}
        for c in (frac(x), 1 - frac(x), frac(al * x), 1 - frac(al * x)):
            if 0 < c < Fr(1, 2): cuts.add(c)
        cuts = sorted(cuts)
        out = []
        for u, v in zip(cuts, cuts[1:]):
            z = (u + v) / 2
            lK = floor(x - z) + floor(x + z) + 1
            lN = floor(al * x - z) + floor(al * x + z) + 1
            out.append((v - u, self.a * lN - lK - (self.s - 1)))
        return out

    def wf(self, x):
        """(Gamma_wf(x), doubled level)"""
        T = self.types(x)
        r = self.lam * x
        w = min(B for mu, B in T)
        tot = Fr(0)
        while True:
            avail = sum(mu for mu, B in T if B <= w and (w - B) % 2 == 0)
            take = min(avail, r)
            tot += take * w
            r -= take
            if r == 0: return tot, w
            w += 1

    def scalN(self, x):
        """N(x): v_p(S_K) = p N(K/p) + O(x),  S_K = (K!)^{2h} 4^{h-1} / ((N!)^{2ah} prod_{i<h} (2i)!^2)."""
        lam, al = self.lam, self.alpha
        m = floor(2 * lam * x)
        J = m * lam * x - Fr(m * (m + 1), 4)
        return 2 * lam * x * floor(x) - 2 * self.a * lam * x * floor(al * x) - 2 * J

    def R(self, x):
        return -self.wf(x)[0] - self.scalN(x)

    def xbreaks(self, x0, x1):
        """superset of the breakpoints of Gamma_wf and N on [x0,x1]."""
        al, lam = self.alpha, self.lam
        S = {x0, x1}
        for coef in (Fr(2), 2 * al, 1 - al, 1 + al, Fr(1), al, 2 * lam):   # coef*x integer
            k = math.ceil(coef * x0)
            while Fr(k) / coef <= x1:
                S.add(Fr(k) / coef); k += 1
        S = sorted(S)
        extra = set()
        for u, v in zip(S, S[1:]):
            xa, xb = u + (v - u) / 3, u + 2 * (v - u) / 3
            Ta, Tb = self.types(xa), self.types(xb)
            assert [B for _, B in Ta] == [B for _, B in Tb]
            wlo = min(B for _, B in Ta)
            whi = max(self.wf(xa)[1], self.wf(xb)[1]) + 8
            for w in range(wlo, whi + 1):
                cnt = [max(0, (w - B) // 2 + 1) for _, B in Ta]
                fa = sum(m * c for (m, _), c in zip(Ta, cnt)) - lam * xa
                fb = sum(m * c for (m, _), c in zip(Tb, cnt)) - lam * xb
                if fa != fb:
                    xs = xa - fa * (xb - xa) / (fb - fa)
                    if u < xs < v: extra.add(xs)
        return sorted(set(S) | extra)

def collinear(f, u, v):
    p = [u + (v - u) * Fr(k, 8) for k in (1, 3, 4, 6, 7)]
    val = [f(t) for t in p]
    sl = (val[-1] - val[0]) / (p[-1] - p[0])
    assert all(val[i] == val[0] + sl * (p[i] - p[0]) for i in range(5)), (float(u), float(v))
    return sl, val[0] - sl * p[0]            # f(t) = sl*t + ic on (u,v)

def int_R_over_x3(mod, x0, x1):
    B = mod.xbreaks(x0, x1)
    tot = Fr(0)
    for u, v in zip(B, B[1:]):
        c1, c0 = collinear(mod.R, u, v)
        tot += c1 * (1 / u - 1 / v) + c0 * (1 / u ** 2 - 1 / v ** 2) / 2
    return tot, len(B) - 1

# ---------------------------------------------------------------- Fauzan's outer functions (s=5, a=6, alpha=3/40)
AL, LA = Fr(3, 40), Fr(37, 40)
pos = lambda t: max(t, Fr(0))
def R0(y):
    if y < Fr(1, 2): return 8 - 9 * y - 8 * AL - 5 * min(AL, 1 - 2 * y) - 5 * pos(1 + AL - 3 * y)
    if y < 1: return 7 * (1 - y) - 6 * min(AL, 1 - y) - 6 * pos(1 + AL - 2 * y) + pos(1 + 4 * AL - 2 * y)
    return Fr(0)
def dd(y):
    return pos(1 + 4 * AL - 3 * y - pos(1 + AL - 3 * y)) if y < Fr(1, 2) else Fr(0)
def Cout(y): return R0(y) - dd(y)
def scal(y): return -2 * LA * floor(1 / y) + sum(pos(2 * LA - j * y) for j in range(1, 6))

if __name__ == "__main__":
    F = Model(5, 6, AL, LA)
    # ---- check 1: reproduce (5.18) from the water-filling continuum function
    I318, npc = int_R_over_x3(F, Fr(3), Fr(20))
    ref = Fr(322437603634266857629, 7535670527041937280000)
    print(f"[check] int_3^20 R/x^3 with R = -Gamma_wf - N : {I318} = {float(I318):.12f}  ({npc} cells);  equals Fauzan (5.18): {I318 == ref}")
    # ---- outer range in y = p/K
    Cwf = lambda y: -y * F.wf(1 / y)[0]
    ys = {1 / x for x in F.xbreaks(Fr(1), Fr(3))} | {Fr(1, 3), Fr(1, 2), Fr(1), 2 * LA}
    ys |= {(1 - AL) / 2, (1 + AL) / 3, (1 + 4 * AL) / 3, 1 - AL, (1 + AL) / 2, (1 + 4 * AL) / 2, Fr(1, 3) + AL}  # (5.8),(5.9)
    ys |= {2 * LA / j for j in range(1, 6)} | {Fr(1, k) for k in (1, 2, 3)}
    ys = sorted(t for t in ys if Fr(1, 3) <= t <= 2 * LA)
    Iold = Inew = Iwf = Fr(0)
    pieces = []
    for u, v in zip(ys, ys[1:]):
        o1, o0 = collinear(lambda t: Cout(t), u, v)
        s1, s0 = collinear(lambda t: scal(t), u, v)
        lin = lambda c1, c0, a_, b_: c1 * (b_ * b_ - a_ * a_) / 2 + c0 * (b_ - a_)
        Iold += lin(o1 + s1, o0 + s0, u, v)
        if v <= 1:
            w1, w0 = collinear(Cwf, u, v)
            Iwf += lin(w1 + s1, w0 + s0, u, v)
            cuts = [u, v]
            if w1 != o1:
                t = (o0 - w0) / (w1 - o1)
                if u < t < v: cuts = [u, t, v]
            for a_, b_ in zip(cuts, cuts[1:]):
                mid = (a_ + b_) / 2
                use_wf = (w1 * mid + w0) < (o1 * mid + o0)
                c1, c0 = (w1, w0) if use_wf else (o1, o0)
                Inew += lin(c1 + s1, c0 + s0, a_, b_)
                pieces.append((a_, b_, 'wf ' if use_wf else 'out', o1, o0, w1, w0))
        else:
            Inew += lin(o1 + s1, o0 + s0, u, v)
    print(f"[check] I_out from (5.8)-(5.10): {Iold}  equals 127751/96000: {Iold == Fr(127751, 96000)}")
    gain = Iold - Inew
    print(f"I_out^new = int min(C_out, C_wf) + scalar = {Inew} = {float(Inew):.9f}")
    print(f"gain  g = I_out - I_out^new = {gain} = {float(gain):.9f}")
    print(f"(all-water-filling on 1/3<y<1 would give {float(Iwf + Iold - Iold):.6f} for the y<1 part; not used)")
    print("\npieces of the outer range (y = p/K; cost per K: C(y) = c1*y + c0):")
    print("    y-interval            x=K/p interval     used | C_out = c1 y + c0        | C_wf = c1 y + c0")
    for a_, b_, wh, o1, o0, w1, w0 in pieces:
        print(f"  [{float(a_):.6f},{float(b_):.6f}]  [{float(1/b_):.5f},{float(1/a_):.5f}]  {wh} | {str(o1):>8s} y + {str(o0):>9s} | {str(w1):>10s} y + {str(w0):>9s}   [{a_}, {b_}]")
    wfp = [(a_, b_) for a_, b_, wh, *_ in pieces if wh == 'wf ']
    if wfp:
        print(f"water filling is used exactly for y in ({min(a for a, b in wfp)}, {max(b for a, b in wfp)}), i.e. "
              f"{float(1/max(b for a, b in wfp)):.6f} < K/p < {float(1/min(a for a, b in wfp)):.6f}")
    # ---- constants
    Astar = Fr(9928298118277006344769, 7535670527041937280000)
    AM = lambda M: Astar + 7 * LA / M - (Fr(2923, 240) - Fr(1, 4)) / M ** 2 + Fr(32, M ** 3)
    print("\nconstants (Fauzan (5.19),(5.20));  new = old - g")
    for nm, val in (("A*", Astar), ("A_200", AM(200)), ("A_400", AM(400)), ("A_100000", AM(100000))):
        print(f"  {nm:9s} old {float(val):.6f}   new {float(val - gain):.6f}   margin vs U=-1.36700: old {float(val) - 1.36700:+.6f} new {float(val - gain) - 1.36700:+.6f}"
              f"   vs -1.3855: old {float(val) - 1.3855:+.5f} new {float(val - gain) - 1.3855:+.5f}")
    print(f"  exact new A* = {Astar - gain}")
    # ---- numerical tail of int_3^X R/x^3 (water filling everywhere): exact on [3,X], then Fauzan's bounds
    if len(sys.argv) > 1:
        X = int(sys.argv[1])
        IX, npc = int_R_over_x3(F, Fr(3), Fr(X))
        P_X = 74 * frac(AL * X) * (1 - frac(AL * X)) - LA * frac(Fr(X)) * (1 - frac(Fr(X)))
        print(f"\n int_3^{X} R/x^3 = {float(IX):.9f} ({npc} cells);  -lam/X = {float(-LA/X):.6f}")
        print(f" => int_3^inf R/x^3 ~ {float(IX - LA/X):.6f} + O(X^-2)   (Fauzan uses (5.18)+(5.16) = {float(ref - Fr(2689,48000)):.6f} as upper bound)")
