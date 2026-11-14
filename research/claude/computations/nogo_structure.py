#!/usr/bin/env python3
"""nogo_structure.py -- exact check of every intermediate claim in the proof of Theorem A (see nogo_newton.py).
Gram matrix in the basis f_j=D_K/(t+j^2):  G(X)=P+diag(c_j m_j(X)),  m_j(X)=j^(s-1)(X+a_j);  X=Y/p^s:
c_j m_j = delta_j (Y-h_j), delta_j=c_j j^(s-1)/p^s, h_j=-p^s a_j.  Checks: P p-integral; v(delta_j); h_j p-integral and
h_j = H_[j/p] mod p (p∤j), = H_d - 1/(2d^s) (j=dp);  theta=min_{j in R} v(P_j*/delta_j) >= s-2m-1;
v_G(omega) <= m(2m-3), omega(Y)=det(P_WW+diag(delta_i(Y-h_i)));  Delta_K = V^-2 det G, lead = +-K!^(s-1).
Usage: python3 nogo_structure.py s K p [s K p ...]"""
import sys
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from flint import fmpq, fmpq_mat, fmpq_poly
from hankel_lib import mu_monomial, vp
def build(s, K):
    H = [fmpq(0)]
    for j in range(1, K + 1): H.append(H[-1] + fmpq(1, j ** s))
    DK = fmpq_poly([1])
    for j in range(1, K + 1): DK = DK * fmpq_poly([j * j, 1])
    mom = [mu_monomial(e, s) for e in range(K)]
    mu = lambda q: sum((c * mom[k] for k, c in enumerate(q.coeffs())), fmpq(0))
    f = [None] + [divmod(DK, fmpq_poly([j * j, 1]))[0] for j in range(1, K + 1)]
    c = [None] + [f[j](fmpq(-j * j)) for j in range(1, K + 1)]
    P = [[None] * (K + 1) for _ in range(K + 1)]
    for i in range(1, K + 1):
        for j in range(i, K + 1):
            q = divmod(f[j] - fmpq_poly([c[j]]), fmpq_poly([j * j, 1]))[0] if i == j else divmod(f[i], fmpq_poly([j * j, 1]))[0]
            P[i][j] = P[j][i] = mu(q)
    a = [None] + [-H[j] + (fmpq(1, 2 * j) - fmpq(1, s - 1)) / fmpq(j ** (s - 1)) for j in range(1, K + 1)]
    return a, c, P, H, mom
def run(s, K, p):
    a, c, P, H, mom = build(s, K); m = K // p; ok = True
    def chk(name, cond):
        nonlocal ok; ok &= bool(cond); print(f"   [{'ok' if cond else 'FAIL'}] {name}")
    print(f"s={s} K={K} p={p} m={m}  (H2: {(s-1)*(p-1)>2*K-2}, H3: {m==1 or s-2*m-1>m*(2*m-3)})")
    chk("moments mu(t^e), e<=K-2, p-integral", all(vp(x, p) >= 0 for x in mom[:K - 1]))
    chk("P p-integral", all(vp(P[i][j], p) >= 0 for i in range(1, K + 1) for j in range(1, K + 1)))
    delta = [None] + [c[j] * fmpq(j ** (s - 1)) / fmpq(p ** s) for j in range(1, K + 1)]
    h = [None] + [-fmpq(p ** s) * a[j] for j in range(1, K + 1)]
    chk("v(delta_j)=[(K-j)/p]+[(K+j)/p]-s <= 2m+1-s (p∤j), =2m-3 (p|j)", all((vp(delta[j], p) == (K - j) // p + (K + j) // p - s and vp(delta[j], p) <= 2 * m + 1 - s) if j % p else vp(delta[j], p) == 2 * m - 3 for j in range(1, K + 1)))
    def modp(q): return (int(q.p) * pow(int(q.q), -1, p)) % p
    chk("h_j p-integral; h_j = H_[j/p] (p∤j), H_d-1/(2d^s) (j=dp) mod p", all(vp(h[j], p) >= 0 and modp(h[j]) == modp(H[j // p] - (fmpq(1, 2 * (j // p) ** s) if j % p == 0 else 0)) for j in range(1, K + 1)))
    W = [j for j in range(1, K + 1) if j % p == 0] if m >= 2 else []
    R = [j for j in range(1, K + 1) if j not in W]
    theta = min(vp(P[j][i], p) - vp(delta[j], p) for j in R for i in range(1, K + 1))
    chk(f"theta={theta} >= s-2m-1={s-2*m-1}" + (" (m=1: >=1)" if m == 1 else ""), theta >= (s - 2 * m - 1 if m > 1 else 1))
    if W:
        Yp = fmpq_poly([0, 1]); n = len(W)
        vals = []
        for y in range(n + 1):
            vals.append(fmpq_mat(n, n, [P[i][j] + (delta[i] * (y - h[i]) if i == j else 0) for i in W for j in W]).det())
        dd = vals[:]
        for k in range(1, n + 1):
            for i in range(n, k - 1, -1): dd[i] = (dd[i] - dd[i - 1]) / fmpq(k)
        om = fmpq_poly([0])
        for k in range(n, -1, -1): om = om * fmpq_poly([-k, 1]) + fmpq_poly([dd[k]])
        vg = min(vp(x, p) for x in om.coeffs() if x != 0)
        chk(f"v_G(omega)={vg} <= m(2m-3)={m*(2*m-3)} ; needed < theta: {vg < theta}", vg <= m * (2 * m - 3))
    nb = sum(1 for j in R if modp(h[j]) != 0)
    print(f"   #(j in R: h_j != 0 mod p) = {nb};  n_p - m = {K-p+1-m} (n_p={K-p+1})   all checks: {ok}")
a = sys.argv[1:]
for i in range(0, len(a), 3): run(int(a[i]), int(a[i + 1]), int(a[i + 2]))
