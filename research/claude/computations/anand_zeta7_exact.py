#!/usr/bin/env python3
"""
anand_zeta7_exact.py -- exact test of the determinant of P. Anand, "zeta(7) is irrational"
(Zenodo 22920911), at K = 40n, N = 3n, h = 37n for small n.

Construction (preprint Def. 2.1, 2.2; eq. (5)-(8)):
    mu_X(t^e)       = (-1)^e B_{2e+2} (2e+3)(2e+4)(2e+5)(2e+6)(2e+7)/720
    mu_X(1/(t+j^2)) = j^6 (X - H_j^(7)) + 1/(2j) - 1/6
    G_K(X) = [ mu_X( D_N(t)^8 t^{i+j} / D_K(t) ) ]_{0<=i,j<h},   Delta_K = det G_K
    S_K = (K!)^{2h} 4^{h-1} / ( (N!)^{16h} prod_{i=1}^{h-1} ((2i)!)^2 ),   F_K = S_K Delta_K.
(hankel_lib.moments_F1 with s=7, a=8 is exactly this functional and matrix.)

For each n the script computes, EXACTLY (python-flint rationals):
  * Delta_K(X) in Q[X];  checks deg = h and the leading coefficient formula (10) of Lemma 2.3;
  * the content of Delta_K, hence the primitive integer polynomial P_K = Delta_K / cont(Delta_K);
  * Gauss valuations v_p(Delta_K) for every prime p (all primes dividing the content are listed);
  * v_p(S_K), v_p(F_K), the preprint's gamma_p^out (59) wherever its hypotheses (52) hold, the crude
    bound (41), and the claim "v_p(Delta_K) >= 0 for p > K" of Prop. 4.6;
and, in ball arithmetic (arb, precision raised until the radius is negligible):
  * log F_K(zeta(7)), log P_K(zeta(7)), sign.
Every Q in Z[X] proportional to Delta_K is an integer multiple of P_K, so |Q(zeta(7))| >= |P_K(zeta(7))|.

No claims are made by this file; it only computes.   Usage: python3 anand_zeta7_exact.py n [n ...]
"""
import sys, math, time, pickle, os
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import (primes_upto, moments_F1, hankel, det_poly_affine, content, fact,
                        eval_poly_at_zeta, log_abs_arb, log_abs_fmpq)
import gmpy2
from flint import fmpq, fmpz, arb, ctx
from fractions import Fraction as Fr
from anand_zeta7_constants import Gamma

S_ZETA, A_EXP = 7, 8


def vq(q, p):
    num = gmpy2.mpz(int(q.p))
    a = gmpy2.remove(num, p)[1]
    if a:
        return int(a)
    return -int(gmpy2.remove(gmpy2.mpz(int(q.q)), p)[1])


def legendre(m, p):
    v = 0
    while m:
        m //= p
        v += m
    return v


def vp_SK(K, N, h, p):
    v = 2 * h * legendre(K, p) - 16 * h * legendre(N, p) - 2 * sum(legendre(2 * i, p) for i in range(1, h))
    if p == 2:
        v += 2 * (h - 1)
    return v


def log_SK(K, N, h):
    return (2 * h * math.lgamma(K + 1) + (h - 1) * math.log(4) - 16 * h * math.lgamma(N + 1)
            - 2 * sum(math.lgamma(2 * i + 1) for i in range(1, h)))


def gamma_out(K, N, h, p):
    """preprint (59), with r_p from (54); None unless hypotheses (52) hold."""
    if not (p >= 11 and 3 * p > K and p <= K and p * p > 2 * K and 7 * N <= 2 * p - 2):
        return None
    v = K - p * (K // p)
    u = max(0, N + v - p + 1)
    tp = min(N, v) + u
    rp = min(h, max(0, K + 6 * N - 3 * p + 4))
    if K < 2 * p:
        return -9 * (K - p) + 8 * tp - 1 - min(rp, p - 1 - N + u)
    return -9 * (K - p) + 3 + 16 * N + 7 * tp - min(rp, p + u)


def ilog(p, x):
    """floor(log_p x) for integer x >= 1"""
    k = 0
    while p ** (k + 1) <= x:
        k += 1
    return k


def crude41(K, h, p):
    """preprint (41): v_p(F_K) >= -8h floor(log_p(5K)) - h v_p(720)"""
    v720 = {2: 4, 3: 2, 5: 1}.get(p, 0)
    return -8 * h * ilog(p, 5 * K) - h * v720


def run(n, cache_dir):
    K, N, h = 40 * n, 3 * n, 37 * n
    t0 = time.time()
    cf = os.path.join(cache_dir, f"delta_n{n}.pkl")
    if os.path.exists(cf):
        with open(cf, 'rb') as f:
            raw = pickle.load(f)
        cs = [fmpq(fmpz(a), fmpz(b)) for a, b in raw]
        how = "exact Delta_K(X) loaded from cache of an earlier run of this script"
    else:
        al, be = moments_F1(S_ZETA, K, N, A_EXP, h)
        cs = det_poly_affine(hankel(al, h), hankel(be, h))
        with open(cf, 'wb') as f:
            pickle.dump([(int(c.p), int(c.q)) for c in cs], f)
        how = f"Delta_K(X) computed exactly in {time.time()-t0:.0f}s"
    print(f"=== n={n}: K={K} N={N} h={h}   ({how}; first runs took 0s, 5s, 34s, 157s for n=1..4)")
    deg = max(k for k, c in enumerate(cs) if c != 0)
    # Lemma 2.3 leading coefficient
    lead = fmpz(1)
    for j in range(N + 1, K + 1):
        DN = fmpz(1)
        for m in range(1, N + 1):
            DN *= (m * m - j * j)
        lead *= fmpz(j) ** 6 * DN ** 7
    if (h * (h - 1) // 2) % 2:
        lead = -lead
    print(f"  degree = {deg} (h = {h});  leading coefficient equals formula (10): {cs[h] == fmpq(lead)}")
    cont = content(cs)
    logcont = log_abs_fmpq(cont)
    lS = log_SK(K, N, h)
    prec = 4000 + 80 * K * K
    while True:
        val = eval_poly_at_zeta(cs, S_ZETA, prec)
        if not val.contains(0) and val.rad() < abs(val.mid()) * 1e-25:
            break
        prec *= 2
    logdet, rad = log_abs_arb(val)
    sign = '+' if val > 0 else '-'
    K2 = K * K
    logF = lS + logdet
    logP = logdet - logcont
    loglead = float(arb(abs(lead)).log()) if True else 0
    print(f"  sign Delta_K(zeta(7)) = {sign}   (arb precision {prec} bits, log-radius {rad:.1e})")
    print(f"  log Delta_K(zeta7)        = {logdet:.6f}")
    print(f"  log S_K                   = {lS:.6f}")
    print(f"  log F_K(zeta7)      /K^2  = {logF / K2:+.6f}      [preprint (108): <= U = +1.76 (+ lower order)]")
    print(f"  log cont(F_K)       /K^2  = {(lS + logcont) / K2:+.6f}      [= -log m_true/K^2; preprint claims m_K,M <= exp(A_M K^2), A_200 = -13.7529]")
    print(f"  log |lead F_K|      /K^2  = {(lS + loglead) / K2:+.6f}      [m * lead(F_K) must be a nonzero integer => log m >= -this]")
    print(f"  log P_K(zeta7)      /K^2  = {logP / K2:+.6f}      [P_K primitive; any Q in Z[X] prop. to Delta_K has |Q(zeta7)| >= P_K(zeta7)]")
    print(f"  log P_K(zeta7)            = {logP:+.3f}    vs  claimed bound -19000 n^2 = {-19000 * n * n}   (= {-19000 / 1600:.4f} K^2)")
    # elementary lower bound (see anand_zeta7_checks.py part 3), (a,b) = (3/40, 1/5)
    a_, b_, al_, lam_ = 3 / 40, 1 / 5, 3 / 40, 37 / 40
    cG = lam_ * (-8 * al_ * math.log(a_) + math.log(1 + b_) + 2 * math.pi * math.sqrt(b_)) - lam_ ** 2 * math.log((b_ - a_) / 4)
    lb = -cG * K2 + h * (math.log((2 * math.pi) ** 6 / 720) + 3 * math.log(a_) + 2 * math.log(2) - math.log(2 * h))
    print(f"  log(Delta_K(zeta7)/|lead Delta_K|) = {logdet - loglead:+.3f}  >=  elementary lower bound {lb:+.3f}  >  -19000 n^2 : {logdet - loglead >= lb > -19000 * n * n}")
    # prime table
    rest = cont
    rows = []
    P = primes_upto(max(2 * h + 60, K + 6 * N + 60))
    viol_out, viol_41, viol_big = [], [], []
    for p in P:
        v = vq(cont, p)
        rest = rest / fmpq(p) ** v
        vS = vp_SK(K, N, h, p)
        go = gamma_out(K, N, h, p)
        c41 = crude41(K, h, p)
        rows.append((p, v, vS, v + vS, go, c41))
        if go is not None and v < go:
            viol_out.append((p, v, go))
        if v + vS < c41:
            viol_41.append((p, v + vS, c41))
        if p > K and v < 0:
            viol_big.append((p, v))
    print(f"  content of Delta_K fully explained by primes <= {P[-1]}: {rest == 1 or rest == -1}")
    print("     p   range | v_p(Delta) true | gamma_out (59) | true-gamma | v_p(S_K) | v_p(F_K) true | crude (41) | p*Gamma(K/p) [continuum form (66),(70) of gamma_in; heuristic: K is outside the range K>=200M^2 of Prop. 4.4]")
    rng = lambda p: 'p^2<=5K' if p * p <= 5 * K else ('<=K/3' if 3 * p <= K else ('outer' if p <= K else ('K<p<=2h' if p <= 2 * h else '>2h')))
    for (p, v, vS, vF, go, c41) in rows:
        if v == 0 and vS == 0 and go is None:
            continue
        g = '      .' if go is None else f"{go:7d}"
        d = '      .' if go is None else f"{v - go:7d}"
        gi = ''
        if p * p > 5 * K and 3 * p <= K:
            gi = f"   | {float(p * Gamma(Fr(K, p))):9.1f}"
        print(f"  {p:5d} {rng(p):8s}| {v:8d}        | {g}        | {d}    | {vS:7d}  | {vF:8d}      | {c41:7d}{gi}")
    print(f"  Prop 4.6 formula (59) violated (true v_p(Delta_K) < gamma_out) at primes satisfying (52): {viol_out if viol_out else 'none'}")
    print(f"  (41) violated (true v_p(F_K) < crude bound): {viol_41 if viol_41 else 'none'}")
    print(f"  'v_p(Delta_K) >= 0 for p > K' violated: {viol_big if viol_big else 'none'}")
    # cost split
    outer = [(p, v, vS, go) for (p, v, vS, vF, go, c41) in rows if go is not None]
    s_true = sum((v + vS) * math.log(p) for (p, v, vS, go) in outer) / K2
    s_claim = sum((go + vS) * math.log(p) for (p, v, vS, go) in outer) / K2
    s_big = sum(vS * math.log(p) for (p, v, vS, vF, go, c41) in rows if p > K) / K2
    s_all = sum(vF * math.log(p) for (p, v, vS, vF, go, c41) in rows) / K2
    print(f"  sum_p v_p(F_K) log p / K^2 over ALL primes (= log cont F_K / K^2)      : {s_all:+.6f}")
    print(f"  outer primes K/3<p<=K with (52): sum (v_p(S_K)+true v_p(Delta)) log p/K^2: {s_true:+.6f}")
    print(f"  outer primes, preprint's exponents L_p = v_p(S_K)+gamma_out : sum/K^2     : {s_claim:+.6f}")
    print(f"  primes K<p<=2h: sum v_p(S_K) log p / K^2  (gamma_out := 0 there)          : {s_big:+.6f}")
    print(f"  # time {time.time()-t0:.0f}s", flush=True)


if __name__ == '__main__':
    cache = '/tmp/claude-0/-home-user-Mathathon/a9b0d4a8-9889-5879-bd9f-5fba824e9410/scratchpad/anand'
    os.makedirs(cache, exist_ok=True)
    for a in sys.argv[1:]:
        run(int(a), cache)
