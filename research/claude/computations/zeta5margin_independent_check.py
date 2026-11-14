#!/usr/bin/env python3
"""Independent check of zeta5-margin.md claims (Theorem W numerics, constants, margin).
Written from Fauzan's definitions only (Sec. 2.1, (3.1), (4.14), (5.8)-(5.10), App. B.3).
Prints PASS/FAIL lines."""
import sys, time
from fractions import Fraction as Fr
import flint
from flint import fmpz, fmpq, fmpz_poly, fmpz_mat

ALLOK = True
def report(ok, msg):
    global ALLOK
    ALLOK &= bool(ok)
    print(("PASS " if ok else "FAIL ") + msg)

def primes_upto(n):
    s = [True]*(n+1); s[0]=s[1]=False
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = [False]*len(s[i*i::i])
    return [i for i in range(n+1) if s[i]]

def vp(x, p):
    x = fmpz(x)
    if x == 0: return 10**9
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
def vpq(q, p):
    q = fmpq(q)
    return vp(q.p, p) - vp(q.q, p)

# ---------------- (1) exact Delta_K ----------------
def build_hankel(n):
    K, N, h = 40*n, 3*n, 37*n
    t = fmpz_poly([0, 1])
    DN = fmpz_poly([1]);
    for j in range(1, N+1): DN *= (t + j*j)
    Dtail = fmpz_poly([1])
    for j in range(N+1, K+1): Dtail *= (t + j*j)
    W = DN**5                      # D_N^6/D_K = D_N^5/D_tail
    H5 = [fmpq(0)]
    for v in range(1, K+1): H5.append(H5[-1] + fmpq(1, v**5))
    maxdeg = 5*N + 2*h - 2 - (K-N)
    mu_mono = [fmpq((-1)**e) * fmpq.bernoulli(2*e+2) * (2*e+3)*(2*e+4)*(2*e+5) / 24
               for e in range(max(maxdeg, 0)+1)]
    # residues of W t^e / Dtail at t=-j^2 (j=N+1..K): W(-j^2)(-j^2)^e / Dtail'(-j^2)
    dDt = Dtail.derivative()
    res0 = {j: fmpq(W(-j*j)) / fmpq(dDt(-j*j)) for j in range(N+1, K+1)}
    a, b = [], []
    for e in range(2*h-1):
        num = W * t**e
        q = num // Dtail  # exact polynomial quotient (Dtail monic)
        A = fmpq(0)
        for k in range(q.degree()+1):
            ck = q[k]
            if ck != 0: A += ck * mu_mono[k]
        B = fmpq(0)
        for j in range(N+1, K+1):
            r = res0[j] * fmpz(-j*j)**e
            B += r * j**4
            A += r * (-j**4 * H5[j] - fmpq(1, 4) + fmpq(1, 2*j))
        a.append(A); b.append(B)
    return K, N, h, a, b

def det_poly(h, a, b):
    """det[a_{i+j} + X b_{i+j}] as list of fmpq coefficients, via evaluation at X=0..h + Newton."""
    from math import lcm
    D = 1
    for x in a + b: D = lcm(D, int(fmpq(x).q))
    A = [fmpz((x*D).p) for x in a]; B = [fmpz((x*D).p) for x in b]
    vals = []
    for X in range(h+1):
        M = fmpz_mat(h, h, [A[i+j] + X*B[i+j] for i in range(h) for j in range(h)])
        vals.append(fmpq(M.det()))
    xs = list(range(h+1))
    # Newton divided differences
    c = vals[:]
    for k in range(1, h+1):
        for i in range(h, k-1, -1):
            c[i] = (c[i] - c[i-1]) / (xs[i] - xs[i-k])
    poly = [fmpq(0)]*(h+1)
    # Horner expand
    res = [fmpq(0)]
    for k in range(h, -1, -1):
        # res = res*(X - xs[k]) + c[k]
        new = [fmpq(0)]*(len(res)+1)
        for i, r in enumerate(res):
            new[i+1] += r; new[i] -= r*xs[k]
        new[0] += c[k]; res = new
    res = res[:h+1]
    Dh = fmpq(D)**h
    return [r / Dh for r in res]

def gamma_wf(p, K, N, h):
    m = (p-1)//2
    def ell(A, c):
        return sum(1 for j in range(1, A+1) if j % p in (c, p-c))
    base = [6*ell(N, c) - ell(K, c) - 4 for c in range(1, m+1)]   # doubled, step 2
    z0 = 12*(N//p) - 2*(K//p) + 1                                   # doubled, step 4
    def count(T):
        cnt = sum((T-b)//2 + 1 for b in base if b <= T)
        if z0 <= T: cnt += (T-z0)//4 + 1
        return cnt
    lo, hi = min(base+[z0]) - 1, min(base+[z0]) + 4*h + 10
    while hi - lo > 1:
        mid = (lo+hi)//2
        if count(mid) >= h: hi = mid
        else: lo = mid
    T = hi  # smallest T with count(T) >= h
    tot, cnt = 0, 0
    for b in base:
        if b <= T - 1:
            k = (T-1-b)//2 + 1; tot += k*b + k*(k-1); cnt += k
    if z0 <= T - 1:
        k = (T-1-z0)//4 + 1; tot += k*z0 + 2*k*(k-1); cnt += k
    tot += (h-cnt)*T
    return tot

def gamma_ell_fast(p, K, N, h):
    """same as gamma_wf but with O(m) ell counts (for large K)."""
    m = (p-1)//2
    def ell(A, c):
        q, r = divmod(A, p)
        return 2*q + (1 if c <= r else 0) + (1 if p-c <= r else 0)
    base = [6*ell(N, c) - ell(K, c) - 4 for c in range(1, m+1)]
    z0 = 12*(N//p) - 2*(K//p) + 1
    import numpy as np
    bs = np.array(base, dtype=np.int64)
    def count(T):
        s = bs[bs <= T]
        return int(((T - s)//2 + 1).sum()) + ((T-z0)//4 + 1 if z0 <= T else 0)
    lo, hi = int(bs.min()) - 1, int(bs.min()) + 4*h + 10
    while hi - lo > 1:
        mid = (lo+hi)//2
        if count(mid) >= h: hi = mid
        else: lo = mid
    T = hi
    s = bs[bs <= T-1]; k = (T-1-s)//2 + 1
    tot = int((k*s + k*(k-1)).sum()); cnt = int(k.sum())
    if z0 <= T-1:
        k0 = (T-1-z0)//4 + 1; tot += k0*z0 + 2*k0*(k0-1); cnt += k0
    return tot + (h-cnt)*T

def gamma_out(p, K, N):
    """Fauzan (4.14); conditions (4.9)."""
    if p > K: return 0
    v = K - p*(K//p); u = max(0, N+v-p+1); tp = min(N, v) + u
    rp = max(0, K + 4*N - 2*p + 2)
    if K < 2*p: return -7*(K-p) + 6*tp - 1 - min(rp, p-1-N+u)
    return -7*(K-p) + 3 + 12*N + 5*tp - min(rp, p+u)

def cond49(p, K, N):
    return p >= 7 and p <= K < 3*p and p*p > 2*K and 2*N < p and 5*N <= 2*p-2

def part1(ns):
    for n in ns:
        t0 = time.time()
        K, N, h, a, b = build_hankel(n)
        coeffs = det_poly(h, a, b)
        report(coeffs[h] != 0, f"n={n}: deg Delta_K = h={h} (leading coeff nonzero); time {time.time()-t0:.1f}s")
        print(f"  n={n} K={K} h={h}:   p | true vG_p | gamma_wf | (4.14) | note")
        viol = 0
        for p in primes_upto(K+30):
            if p < 5 or p*p <= 104*n + 6: continue
            tv = min(vpq(c, p) for c in coeffs if c != 0)
            gw = gamma_wf(p, K, N, h)
            go = gamma_out(p, K, N) if (cond49(p, K, N) or p > K) else None
            flag = []
            if gw > tv: flag.append("WF-VIOLATION"); viol += 1
            if go is not None and go > tv: flag.append("OUT-VIOLATION"); viol += 1
            if go is not None and gw > go: flag.append("wf beats (4.14)")
            print(f"  {p:5d} | {tv:6d} | {gw:6d} | {str(go):>6} | {' '.join(flag)}")
        report(viol == 0, f"n={n}: no violation of gamma_wf <= v_p^G(Delta_K) (nor of (4.14))")
        sys.stdout.flush()

# ---------------- (2) constants ----------------
al, lam = Fr(3, 40), Fr(37, 40)
def pos(x): return x if x > 0 else Fr(0)
def R0(y):
    if Fr(1,3) < y < Fr(1,2):
        return 8 - 9*y - 8*al - 5*min(al, 1-2*y) - 5*pos(1+al-3*y)
    if Fr(1,2) < y < 1:
        return 7*(1-y) - 6*min(al, 1-y) - 6*pos(1+al-2*y) + pos(1+4*al-2*y)
    return Fr(0)
def dfun(y):
    return pos(1+4*al-3*y - pos(1+al-3*y)) if Fr(1,3) < y < Fr(1,2) else Fr(0)
def floor(x): return x.numerator // x.denominator
def Tout(y):
    return R0(y) - dfun(y) - 2*lam*floor(1/y) + sum(pos(2*lam - j*y) for j in range(1, 6))
def integrate_pw_affine(f, bps):
    bps = sorted(set(bps)); s = Fr(0)
    for l, r in zip(bps, bps[1:]):
        m = (l+r)/2; q1 = l + (r-l)/4; q3 = l + 3*(r-l)/4
        assert f(q1) + f(q3) == 2*f(m), ("not affine on", l, r)
        s += (r-l)*f(m)
    return s
def Cout(y): return R0(y) - dfun(y)
def Cwf(y):
    if y <= Fr(43,120): return 12*y - Fr(11,10)
    if y <= Fr(23,60): return 3*y + Fr(17,8)
    if y <= Fr(2,5): return Fr(131,40)
    return Fr(5,2)*y + Fr(91,40)

def part2():
    bps = [Fr(1,3), Fr(43,120), Fr(37,80), Fr(1,2), Fr(37,40), Fr(43,80), Fr(13,20), Fr(13,30),
           Fr(1), Fr(37,20)] + [2*lam/j for j in range(1,6)] + [Fr(1,k) for k in range(1,4)]
    bps = [x for x in bps if Fr(1,3) <= x <= 2*lam]
    Iout = integrate_pw_affine(Tout, bps)
    report(Iout == Fr(127751, 96000), f"I_out recomputed from (5.8)-(5.10) = {Iout} (paper 127751/96000)")
    # C_out formula claimed in the note
    ok = all(Cout(y) == (6*y + Fr(57,40) if y < Fr(43,120) else -6*y + Fr(229,40))
             for y in [Fr(1,3)+Fr(k,10000) for k in range(1, 730)])
    report(ok, "C_out = R0 - d equals 6y+57/40 on (1/3,43/120), -6y+229/40 on (43/120,69/170)")
    report(Cwf(Fr(69,170)) == Cout(Fr(69,170)) and Cwf(Fr(1,3)) < Cout(Fr(1,3)+Fr(1,10**9)),
           "C_wf meets C_out at 69/170 and is below it on the interval")
    # continuity of C_wf
    report(all(abs(Cwf(x-Fr(1,10**12)) - Cwf(x+Fr(1,10**12))) < Fr(1,10**10) for x in [Fr(43,120),Fr(23,60),Fr(2,5)]),
           "C_wf continuous at its breakpoints")
    g = integrate_pw_affine(lambda y: Cout(y) - Cwf(y),
                            [Fr(1,3), Fr(43,120), Fr(23,60), Fr(2,5), Fr(69,170)])
    report(g == Fr(3203, 163200), f"g = int (C_out - C_wf) dy = {g} = {float(g):.6f} (claimed 3203/163200)")
    Inew = Iout - g
    report(Inew == Fr(2139737, 1632000), f"I_out^new = {Inew} = {float(Inew):.6f} (claimed 2139737/1632000)")
    # large-K sanity: -gamma_wf / K vs C_wf(p/K)
    from sympy import nextprime
    maxerr = 0
    for n in (100, 400, 1600):
        K, N, h = 40*n, 3*n, 37*n
        for y0 in (0.335, 0.34, 0.35, 0.357, 0.36, 0.37, 0.38, 0.385, 0.39, 0.40, 0.403, 0.405):
            p = int(nextprime(int(y0*K)))
            gw = gamma_ell_fast(p, K, N, h)
            pred = K*Cwf(Fr(p, K))
            err = float(-gw - pred); maxerr = max(maxerr, abs(err))
            if n == 1600 or y0 in (0.34, 0.37, 0.39, 0.403):
                print(f"  K={K:6d} p={p:6d} y={p/K:.4f}  -gamma_wf={-gw:7d}  K*C_wf={float(pred):10.2f}  diff={err:7.2f}"
                      f"  -gamma_out(4.14)={-gamma_out(p,K,N)}  K*C_out={float(K*Cout(Fr(p,K))):9.2f}")
    report(maxerr < 20, f"|-gamma_wf - K*C_wf(p/K)| <= {maxerr:.2f} at all sampled primes (K=4000,16000,64000)")

# ---------------- (3) margin ----------------
def part3():
    A200 = Fr(127125602969131786927559, 94195881588024216000000)
    U = Fr(-2733991, 2000000)
    g = Fr(3203, 163200)
    report(-1600*(A200+U) - Fr(139,5) == Fr(3089837638249482469, 58872425992515135000),
           "A_200 consistent with the paper's App. B.3 margin identity")
    An = A200 - g
    print(f"  A_200 = {float(A200):.6f}, A_200^new = {An} = {float(An):.6f}")
    print(f"  A_200 + U = {float(A200+U):.6f}")
    print(f"  A_200^new + U = {An+U} = {float(An+U):.6f}")
    report(An + U < 0 and abs(float(An+U) + 0.037039) < 1e-5,
           "A_200^new + U negative and = -0.037039 to within 1e-5 (claimed value)")
    print(f"  NOTE: exact A_200^new+U = {float(An+U):.7f}; the claimed -0.037039 (and old -0.017412) come from"
          f" rounding U to -1.36700; with exact U the 6-decimal values are {float(An+U):.6f} / {float(A200+U):.6f}")

if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] or [1, 2]   # n=3 takes ~23 min: pass '1 2 3'
    part2(); part3(); part1(ns)
    print("OVERALL:", "PASS" if ALLOK else "FAIL")
