#!/usr/bin/env python3
"""nogo7_wild.py -- exact data for the Schur-complement step of Theorem A (zeta7-status.md section 5) with the
hypothesis (H3) replaced by what the proof really uses, and the 'descent' congruence for the wild block.
Basis f_j=D_K/(t+j^2): G(Y)=P+diag(delta_j(Y-h_j)), X=Y/p^s, W={p,..,mp}, R=rest, omega=det G_WW(Y), g=v_Gauss(omega).
The proof of Theorem A works verbatim whenever  (C1) theta''=min_{j in R,j'} v(P_jj')-v(delta_j) >= 1  and
(C2) theta'=min_{j in R} v_Gauss(P_jW adj G_WW)-v(delta_j) > g   [then reduction of F_Y = const*omega_bar*prod_R(Y-h_j)].
Printed per prime p>sqrt(2K), p>s, m=[K/p]>=2: min v(P), g, a-priori m(2m-3), predicted g*=min_k 2k(k-1)+(2m-3)(m-k),
theta', theta'', C1&C2, and for m=2 the congruence P_ww' = mu^(s-4)(D_r) mod p (r=K-2p)  [EXACT COMPUTATION].
Usage: python3 nogo7_wild.py s K [s K ...]"""
import sys
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from flint import fmpq, fmpq_mat, fmpq_poly
from hankel_lib import mu_monomial, vp, primes_upto
def build(s, K):
    H = [fmpq(0)]
    for j in range(1, K + 1): H.append(H[-1] + fmpq(1, j ** s))
    DK = fmpq_poly([1])
    for j in range(1, K + 1): DK = DK * fmpq_poly([j * j, 1])
    mom = [mu_monomial(e, s) for e in range(K)]
    mu = lambda q: sum((c * mom[k] for k, c in enumerate(q.coeffs())), fmpq(0))
    f = [None]; c = [None]; phi = [None]; Pd = [None]
    for j in range(1, K + 1):
        fj = divmod(DK, fmpq_poly([j * j, 1]))[0]
        cj = fj(fmpq(-j * j))
        f.append(fj); c.append(cj); phi.append(mu(fj))
        Pd.append(mu(divmod(fj - fmpq_poly([cj]), fmpq_poly([j * j, 1]))[0]))
    P = [[None] * (K + 1) for _ in range(K + 1)]
    for i in range(1, K + 1):
        P[i][i] = Pd[i]
        for j in range(i + 1, K + 1):
            P[i][j] = P[j][i] = (phi[i] - phi[j]) / fmpq(j * j - i * i)
    a = [None] + [-H[j] + (fmpq(1, 2 * j) - fmpq(1, s - 1)) / fmpq(j ** (s - 1)) for j in range(1, K + 1)]
    return a, c, P, H, mom
def pdet(M):
    n = len(M)
    if n == 0: return fmpq_poly([1])
    if n == 1: return M[0][0]
    tot = fmpq_poly([0])
    for k in range(n):
        if M[0][k] == 0: continue
        sub = [row[:k] + row[k + 1:] for row in M[1:]]
        tot += (M[0][k] * pdet(sub)) * (1 if k % 2 == 0 else -1)
    return tot
def gv(q, p):
    cs = [x for x in q.coeffs() if x != 0]
    return min(vp(x, p) for x in cs) if cs else 10**9
def modp(q, p):
    q = fmpq(q); return (int(q.p) * pow(int(q.q), -1, p)) % p
def run(s, K):
    a, c, P, H, mom = build(s, K)
    for p in primes_upto(K):
        if p * p <= 2 * K or p <= s: continue
        m = K // p; r = K - m * p
        if m < 2: continue
        delta = [None] + [c[j] * fmpq(j ** (s - 1)) / fmpq(p ** s) for j in range(1, K + 1)]
        h = [None] + [-fmpq(p ** s) * a[j] for j in range(1, K + 1)]
        W = [j for j in range(1, K + 1) if j % p == 0]; R = [j for j in range(1, K + 1) if j % p]
        Pint = min(vp(P[i][j], p) for i in range(1, K + 1) for j in range(1, K + 1))
        G = [[fmpq_poly([P[i][j]]) + (fmpq_poly([-delta[i] * h[i], delta[i]]) if i == j else fmpq_poly([0])) for j in W] for i in W]
        om = pdet(G); g = gv(om, p)
        adj = [[pdet([[G[rr][cc] for cc in range(m) if cc != i] for rr in range(m) if rr != k]) * (1 if (i + k) % 2 == 0 else -1) for k in range(m)] for i in range(m)]
        th1 = th2 = 10 ** 9
        for j in R:
            vd = vp(delta[j], p)
            for k in range(m):
                u = fmpq_poly([0])
                for i in range(m): u += adj[i][k] * P[j][W[i]]
                th1 = min(th1, gv(u, p) - vd)
            th2 = min(th2, min(vp(P[j][jj], p) for jj in R) - vd)
        gstar = min(2 * k * (k - 1) + (2 * m - 3) * (m - k) for k in range(m + 1))
        extra = ''
        if m == 2 and s >= 7:
            Dr = fmpq_poly([1])
            for k in range(1, r + 1): Dr = Dr * fmpq_poly([k * k, 1])
            cr = sum((cc * mu_monomial(k, s - 4) for k, cc in enumerate(Dr.coeffs())), fmpq(0))
            okd = Pint >= 0 and vp(cr, p) >= 0 and all(modp(P[x][y], p) == modp(cr, p) for x in W for y in W)
            extra = f" | descent P_WW=mu^({s-4})(D_{r})*J mod p: {okd}, v_p(mu^({s-4})(D_r))={vp(cr, p)}, p|2^s+1: {(2**s+1) % p == 0}"
        print(f"s={s:2d} K={K:3d} p={p:3d} m={m} r={r:2d} unsat(2m<=s-3)={2*m<=s-3!s:5}: minv(P)={Pint:2d} g={g:2d} (H3-bound {m*(2*m-3)}, g*={gstar}) theta'={th1:2d} theta''={th2:2d} C1&C2={th1>g and th2>=1}{extra}", flush=True)
args = sys.argv[1:]
for i in range(0, len(args), 2): run(int(args[i]), int(args[i + 1]))
