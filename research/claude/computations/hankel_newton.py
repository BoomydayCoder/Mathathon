#!/usr/bin/env python3
"""
hankel_newton.py -- the "two-sided Newton basis" for the Hankel determinant

    G_K(X) = [ mu_X( t^{i+j} / D_K(t) ) ]_{0<=i,j<K},   D_K(t) = prod_{j<=K}(t+j^2),

(W = 1, h = K), for odd s.  With the unitriangular bases
    f_i(t)  = prod_{j=1}^{i} (t+j^2)            (rows,    i  = 0..K-1)
    g_i'(t) = prod_{j=K-i'+1}^{K} (t+j^2)       (columns, i' = 0..K-1)
one has det G_K = det C_K,  C_K = [c_{i,i'}],  c_{i,i'} = mu_X(f_i g_i' / D_K):

    c_{i,i'} = mu_X( 1 / prod_{j=i+1}^{K-i'} (t+j^2) )      if i+i' <= K-1   ("pole entry",
               a divided difference of m_j(X) = mu_X(1/(t+j^2)) over the window [i+1, K-i'])
    c_{i,i'} = mu( prod_{j=K-i'+1}^{i} (t+j^2) )            if i+i' >= K     ("polynomial entry").

This script
  (1) checks det C_K == det G_K exactly (small K),
  (2) checks the entrywise p-adic lower bounds used in the report (s = 3):
        (B1)  sqrt(2K)+3 < p :  v_p(c_{i,i'}) >= -3 - 2 N_p(i,i'),
              N_p(i,i') = floor((K-i'-1)/p) - floor((i+1)/p),
        (B2)  all p:  v_p(c_{i,i'}) >= -7 floor(log_p(3K)) - 2 sum_k N_{p^k}(i,i'),
        (B3)  K/2 < p <= K: pole entries are p-integral if K-i' < p and have v_p >= -5 otherwise;
              polynomial entries are p-integral if i+i'-K <= p-3 and have v_p >= -1 otherwise,
        (B4)  p > K: all entries p-integral,
  (3) compares  true v_p(det G_K)  with  the tropical determinant of the valuation matrix
      and with the closed bounds  -3K + 2 floor(K/p)  and  -6(K-p+1).
Usage:  python3 hankel_newton.py [s] [K1,K2,...]

-----------------------------------------------------------------------------------------------
WHY THE BOUNDS HOLD (s = 3; proof notes for the statements this script tests numerically).
Notation: m_j = j^2 (X - H_j) - 1/2 + 1/(2j), H_j = sum_{v<=j} v^-3; R(A,B) = mu_X(1/prod_{j=A+1}^{B}(t+j^2))
= sum_{j in W} rho_j m_j, W = [A+1,B], rho_j = 1/prod_{m in W, m != j} (m^2-j^2); v_p = Gauss valuation.
tau(Q) = L(Q')/2, L(x^k) = B_k (B_1 = -1/2); mu(P) = tau(x^3 P(-x^2)).

(L1) reflection congruence: for an odd prime p and 0<a<p:  m_a == m_{p-a} (mod p).
     [ (p-a)^2 == a^2;  H_{p-a} == H_{p-1} + H_{a-1} == H_a - a^-3 (mod p) as H_{p-1} == 0;  1/(2(p-a)) == -1/(2a). ]
(L2) windows below p: B <= p-1  =>  R(A,B) in Z_(p)[X].
     [ p | m^2-j^2 (m != j in W) iff m = p-j, exactly once; unpaired terms are integral; a pair {a,p-a}
       contributes (m_a/Pi' - m_{p-a}/Pi'')/((p-a)^2-a^2) with Pi' == Pi'' units, numerator == 0 by (L1). ]
(L3) = (B1) for pole entries: p^2-p > 2K, N = #{multiples of p in [A+2,B-1]}  =>  v_p(R(A,B)) >= -3-2N.
     [ For m != j in W, v_p(m^2-j^2) = [p|m-j] + [p|m+j].  A class {m == +-j} (p not| j) has l <= 2N+2 elements
       (2 per full period between interior multiples, <= 2 in each end piece); v_p(m_j) >= -3; so the term
       has v_p >= -3-(l-1) >= -3-2N unless l = 2N+2.  In that case write j = |a+pq|; the index set Q is stable
       under q -> -1-q, m^2-j^2 = p(q'-q)(2a+p(q+q')), p^3 m_j == -a^2 H_{d(q)} (mod p) with d(q) = floor(j/p)
       = q (q>=0), -q-1 (q<0), and the leading term of the class sum is a unit times the divided difference
       sum_q H_{d(q)} / prod_{q' != q}(q'-q), which vanishes because the function is invariant and the
       node set is symmetric with an EVEN number of nodes.  Multiples of p: v_p(m_j) >= -1 and
       v_p(rho_j) >= -2(N+1). ]
(L4) = (B1) for polynomial entries: p > sqrt(2K)+3, P = prod_{j=B+1}^{A}(t+j^2), n = #{multiples of p in [B,A+1]}
     =>  v_p(mu(P)) >= 2n-3.
     [ n <= 1: von Staudt (v_p(mu(t^e)) >= -1).  n >= 2: Raabe's multiplication theorem gives
       tau(Q) = p^-2 sum_{a mod p} tau(Q(a+px)), Q = x^3 P(-x^2).  Q(a+px) = +- p^{l_a} T_a(x) U_a(x),
       T_a = prod (x - q_r) over the roots r == a, U_a in Z[px] with unit constant term, l_a >= 2(n-1);
       tau(x^m) = m B_{m-1}/2 is p-integral for m <= 2p-2 and has v_p >= -1 always, so v_p(tau(Q(a+px))) >= l_a;
       if l_a = 2(n-1) the roots q_r are symmetric under q -> -1-q, T_a(-1-x) = T_a(x), and
       tau(g(-1-x)) = -tau(g(x)) forces tau(T_a) = 0, which gains one more p.  a = 0 gains p^3 from x^3. ]
(L5) weights: v_p(c_{i,i'}) >= u_i + u'_i' for all i,i'  =>  v_p(det) >= sum u_i + sum u'_i'.
     With u_i = 2 floor((i+1)/p), u'_i' = -3 - 2 floor((K-i'-1)/p):  v_p(det G_K) >= -3K + 2 floor(K/p).
(L6) K/2 < p <= K: columns i' <= K-p get weight -5 [(L3) with N <= 1], rows i >= p-1 weight -1
     [mu(t^e) is p-integral for e <= p-3], all other entries are p-integral [(L2)]:  v_p(det G_K) >= -6(K-p+1).
(L7) p > K: every entry is p-integral, v_p(det G_K) >= 0.
(L8) = (B2), all p (crude): termwise for pole entries; for polynomial entries tau(binom(x,k)) = (-1)^{k-1} H_k^{(1)}/(k+1)
     and p^M | Q(x) for all integers x with M >= 2 sum_k n_{p^k} - 4 floor(log_p K).  =>  v_p(det G_K) >= -7K floor(log_p 3K).
Consequently  Phi_K det G_K in Z[X]  with  log Phi_K = 6 int_{K/2}^{K} theta(u) du + O(K^{3/2})   (theta = Chebyshev).
-----------------------------------------------------------------------------------------------
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import *


def newton_matrix(s, K):
    """entries as pairs (alpha, beta): c = alpha + beta*X ; returns (Al, Be) lists of lists of fmpq."""
    H = harmonic_s(K, s)
    # R[A][B] = mu_X(1/prod_{j=A+1}^{B}(t+j^2)), 0 <= A < B <= K, by the divided-difference recursion
    Ra = [[None] * (K + 2) for _ in range(K + 1)]
    Rb = [[None] * (K + 2) for _ in range(K + 1)]
    for j in range(1, K + 1):
        Ra[j - 1][j] = -(fmpq(j) ** (s - 1)) * H[j] - fmpq(1, s - 1) + fmpq(1, 2 * j)
        Rb[j - 1][j] = fmpq(j) ** (s - 1)
    for L in range(2, K + 1):
        for A in range(0, K - L + 1):
            B = A + L
            den = fmpq(B * B - (A + 1) * (A + 1))
            Ra[A][B] = (Ra[A][B - 1] - Ra[A + 1][B]) / den
            Rb[A][B] = (Rb[A][B - 1] - Rb[A + 1][B]) / den
    mu_mon = [mu_monomial(e, s) for e in range(0, K + 1)]
    Al = [[fmpq(0)] * K for _ in range(K)]
    Be = [[fmpq(0)] * K for _ in range(K)]
    # polynomial entries: T[B][A] = mu(prod_{j=B+1}^{A}(t+j^2)), B <= A
    T = {}
    for B in range(0, K + 1):
        P = fmpz_poly([1])
        T[(B, B)] = mu_mon[0]
        for A in range(B + 1, K):
            P = P * fmpz_poly([A * A, 1])
            cs = P.coeffs()
            T[(B, A)] = sum((fmpq(int(cs[k])) * mu_mon[k] for k in range(len(cs))), fmpq(0))
    for i in range(K):
        for ip in range(K):
            A, B = i, K - ip
            if A + 1 <= B:
                Al[i][ip] = Ra[A][B]
                Be[i][ip] = Rb[A][B]
            else:
                Al[i][ip] = T[(B, A)]
    return Al, Be


def val_entry(a, b, p):
    v = 10 ** 6
    if a != 0:
        v = min(v, vp(a, p))
    if b != 0:
        v = min(v, vp(b, p))
    return v


def N_p(i, ip, K, q):
    return (K - ip - 1) // q - (i + 1) // q


def check(s, K, do_det=True, verbose=True):
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    Al, Be = newton_matrix(s, K)
    out = {}
    if do_det:
        A = fmpq_mat(K, K, [Al[i][j] for i in range(K) for j in range(K)])
        B = fmpq_mat(K, K, [Be[i][j] for i in range(K) for j in range(K)])
        cs = det_poly_affine(A, B)
        al, be = moments_F1(s, K, 0, 0, K)
        cs0 = det_poly_affine(hankel(al, K), hankel(be, K))
        out['det_equal'] = all(x == y for x, y in zip(cs, cs0))
        out['cs'] = cs
    viol = {'B1': 0, 'B2': 0, 'B3': 0, 'B4': 0}
    slack = {'B1': 10 ** 6, 'B2': 10 ** 6}
    rows = []
    for p in primes_upto(3 * K):
        V = np.zeros((K, K))
        for i in range(K):
            for ip in range(K):
                a, b = Al[i][ip], Be[i][ip]
                v = val_entry(a, b, p)
                V[i, ip] = min(v, 10 ** 4)
                A_, B_ = i, K - ip
                if a == 0 and b == 0:
                    continue
                # (B2)
                bound2 = -7 * int(math.floor(math.log(3 * K) / math.log(p) + 1e-12))
                q = p
                while q <= 3 * K:
                    bound2 -= 2 * N_p(i, ip, K, q)
                    q *= p
                if v < bound2:
                    viol['B2'] += 1
                slack['B2'] = min(slack['B2'], v - bound2)
                if p > math.sqrt(2 * K) + 3 and p <= K:
                    bound1 = -s - 2 * N_p(i, ip, K, p)
                    if v < bound1:
                        viol['B1'] += 1
                    slack['B1'] = min(slack['B1'], v - bound1)
                if 2 * p > K and p <= K:
                    if A_ + 1 <= B_:
                        bd = 0 if B_ < p else -(s + 2)
                    else:
                        bd = 0 if (A_ - B_) <= p - s else -1
                    if v < bd:
                        viol['B3'] += 1
                if p > K and v < 0:
                    viol['B4'] += 1
        if do_det and p <= K:
            r, c = linear_sum_assignment(V)
            trop = int(V[r, c].sum())
            true = gauss_val(out['cs'], p)
            b1 = -s * K + 2 * (K // p) if p > math.sqrt(2 * K) + 3 else None
            b3 = -(s + 3) * (K - p + 1) if 2 * p > K else None
            rows.append((p, true, trop, b1, b3))
    out['viol'] = viol
    out['slack'] = slack
    out['rows'] = rows
    return out


if __name__ == '__main__':
    s = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    Ks = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [12, 20, 30, 41, 60]
    for K in Ks:
        o = check(s, K)
        print(f"s={s} K={K}: det(Newton)==det(monomial): {o['det_equal']};  violations of (B1),(B2),(B3),(B4): {o['viol']};  min slack (B1),(B2): {o['slack']}")
        print("   p | true v_p(det) | tropical det | -sK+2[K/p] (p>sqrt(2K)+3) | -(s+3)(K-p+1) (p>K/2)")
        for (p, true, trop, b1, b3) in o['rows']:
            flag = ''
            if b1 is not None and true < b1: flag += ' !!B1-det'
            if b3 is not None and true < b3: flag += ' !!B3-det'
            print(f"   {p:3d} | {true:6d} | {trop:6d} | {str(b1):>6} | {str(b3):>6}{flag}")
        sys.stdout.flush()
