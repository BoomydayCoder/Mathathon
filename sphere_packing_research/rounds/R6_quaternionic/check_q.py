#!/usr/bin/env python3
"""Numerical sanity checks for R6_quaternionic/proof.md  (runs in a few seconds).

 (1) Hurwitz units (24) acting on H^k = R^{4k} by RIGHT multiplication:
     D_k = dim of invariant symmetric forms = k(2k-1) for k=1,2,3;
     q(x) = |P(x (x) x)|^2/|x|^4 = 1/4 for random x (P = average over the 24 units);
     P(x (x) x) = L_{x x^*/4}  (left multiplication by the Hermitian matrix x x^*/4).
 (2) Icosians (120 units, two real places): D_k = 2k(2k-1) for k=1,2;
     q(x) = (u_1^2+u_2^2)/4.
 (3) Angular lemma:  E exp((b s/2) sum (u_v-1/s)^2) <= (1-b/a)^(-b s), u ~ Dirichlet(a,..,a).
 (4) Integral lemma (Lemma I): ratio of the double integral to (8/n^2) exp(n^2 S/8).
"""
import numpy as np
from itertools import product
from scipy import integrate, special

rng = np.random.default_rng(1)

def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return np.array([a0*b0 - a1*b1 - a2*b2 - a3*b3,
                     a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1,
                     a0*b3 + a1*b2 - a2*b1 + a3*b0])

def qconj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

E = np.eye(4)
def Rmat(u):   # matrix of x -> x u on R^4
    return np.column_stack([qmul(E[i], u) for i in range(4)])
def Lmat(u):   # matrix of x -> u x on R^4
    return np.column_stack([qmul(u, E[i]) for i in range(4)])

def hurwitz_units():
    U = []
    for i in range(4):
        for s in (1, -1):
            U.append(s*E[i])
    for signs in product((1, -1), repeat=4):
        U.append(0.5*np.array(signs, dtype=float))
    return U

def sym_basis(n):
    B = []
    for i in range(n):
        for j in range(i, n):
            M = np.zeros((n, n))
            if i == j:
                M[i, i] = 1.0
            else:
                M[i, j] = M[j, i] = 1/np.sqrt(2)
            B.append(M)
    return B

def invariant_dim(Rs, n):
    """rank of the projector B -> mean_g R_g B R_g^T on Sym_n (trace inner product)."""
    B = sym_basis(n)
    tr = 0.0
    for M in B:
        PM = sum(R @ M @ R.T for R in Rs)/len(Rs)
        tr += np.sum(PM*M)
    return tr

print("== (1) Hurwitz ==")
U = hurwitz_units()
assert len(U) == 24
# closure check
for a in U:
    for b in U:
        c = qmul(a, b)
        assert any(np.allclose(c, d) for d in U)
for k in (1, 2, 3):
    n = 4*k
    Rs = [np.kron(np.eye(k), Rmat(u)) for u in U]      # x=(x_1..x_k) -> (x_1 u,...,x_k u)
    for R in Rs:
        assert np.allclose(R @ R.T, np.eye(n))
    D = invariant_dim(Rs, n)
    qs = []
    errL = 0.0
    for _ in range(20):
        x = rng.standard_normal(n)
        X = np.outer(x, x)
        PX = sum(R @ X @ R.T for R in Rs)/24
        qs.append(np.sum(PX*PX)/np.dot(x, x)**2)
        # L_{x x^*/4}: block (i,j) = Lmat(x_i conj(x_j))/4
        xs = x.reshape(k, 4)
        LX = np.block([[Lmat(qmul(xs[i], qconj(xs[j])))/4 for j in range(k)] for i in range(k)])
        errL = max(errL, np.abs(LX - PX).max())
    print(f"k={k}: n={n}, D_k={D:.6f} (k(2k-1)={k*(2*k-1)}), q in [{min(qs):.12f},{max(qs):.12f}], "
          f"max|P(x(x)x)-L_(xx*/4)|={errL:.2e}")

print("== (2) icosians ==")
tau = (1+np.sqrt(5))/2; taup = (1-np.sqrt(5))/2
# generators as pairs (sigma_1(u), sigma_2(u)); sigma_2: sqrt5 -> -sqrt5
gens = [(np.array([.5, .5, .5, .5]), np.array([.5, .5, .5, .5])),
        (np.array([tau/2, (tau-1)/2, .5, 0.]), np.array([taup/2, (taup-1)/2, .5, 0.])),
        (np.array([0., 1., 0., 0.]), np.array([0., 1., 0., 0.]))]
def key(p):
    return tuple(np.round(np.concatenate(p), 8) + 0.0)
group = {key(g): g for g in gens}
frontier = list(gens)
while frontier:
    new = []
    for a in frontier:
        for g in gens:
            c = (qmul(a[0], g[0]), qmul(a[1], g[1]))
            kc = key(c)
            if kc not in group:
                group[kc] = c; new.append(c)
    frontier = new
I2 = list(group.values())
print("order of the generated group:", len(I2))
assert len(I2) == 120
for k in (1, 2):
    n = 8*k
    # coordinates: x = (x^(1), x^(2)), x^(v) in H^k
    Rs = [np.block([[np.kron(np.eye(k), Rmat(u1)), np.zeros((4*k, 4*k))],
                    [np.zeros((4*k, 4*k)), np.kron(np.eye(k), Rmat(u2))]]) for (u1, u2) in I2]
    D = invariant_dim(Rs, n)
    err = 0.0
    for _ in range(10):
        x = rng.standard_normal(n)
        X = np.outer(x, x)
        PX = sum(R @ X @ R.T for R in Rs)/120
        q = np.sum(PX*PX)/np.dot(x, x)**2
        u1 = np.dot(x[:4*k], x[:4*k])/np.dot(x, x)
        err = max(err, abs(q - (u1**2 + (1-u1)**2)/4))
    print(f"k={k}: n={n}, D_k={D:.6f} (2k(2k-1)={2*k*(2*k-1)}), max|q-(u1^2+u2^2)/4|={err:.2e}")

print("== (3) angular lemma (Dirichlet), exact LHS for s=2 by quadrature ==")
for (k, z) in [(50, 3.0), (200, 6.0), (1000, 8.0), (10000, 12.0)]:
    a = 2*k; s = 2
    b = int(np.ceil(2*z))
    logB = special.betaln(a, a)
    f = lambda u: np.exp((a-1)*np.log(u*(1-u)) - logB + z*s*2*(u-0.5)**2)
    lhs = integrate.quad(f, 0, 1, points=[0.5], limit=400)[0]
    rhs = (1-b/a)**(-b*s)
    print(f"k={k}, z={z}: E exp(z s sum eta^2) = {lhs:.6f} <= bound {rhs:.6f}")
    assert lhs <= rhs
# Monte Carlo for s=3
for (a, z) in [(100, 3.0), (1000, 5.0)]:
    s = 3; b = int(np.ceil(2*z))
    u = rng.dirichlet([a]*s, size=400000)
    lhs = np.mean(np.exp(z*s*np.sum((u-1/s)**2, axis=1)))
    print(f"s=3, a={a}, z={z}: MC {lhs:.5f} <= bound {(1-b/a)**(-b*s):.5f}")

print("== (4) integral lemma: ratio LHS / ((8/n^2) e^{n^2 S/8}) ==")
def logPhi(y):
    return np.log(0.5) if y <= 0 else min(np.log(0.5), -y*y/2 - np.log(np.sqrt(2*np.pi)*y))
def lemmaI(n, S, Y):
    """(n^2/8) e^{-n^2 S/8} * int_0^S (n sqrt(s)/2) int_0^Y Phi(y)(1-y sqrt s)^(-(n+2)/2) dy ds  (log-space)."""
    m = (n+2)/2
    E0 = n*n*S/8
    def inner(s):
        g = lambda y: np.exp(logPhi(y) - m*np.log1p(-y*np.sqrt(s)) - E0)
        a = n*np.sqrt(s)/2
        pts = sorted(set([1.0, max(1.0, a-8), a, a+8]))
        pts = [p for p in pts if 0 < p < Y]
        return n*np.sqrt(s)/2*integrate.quad(g, 0, min(Y, a+40), points=pts, limit=400)[0]
    val = integrate.quad(inner, 0, S, limit=200, points=[S*0.5, S*0.9, S*0.99])[0]
    return val*n**2/8
for n in (10**3, 10**5, 10**7, 10**9):
    S = 16*np.log(n)/n**2          # rescaled time S = q T  (n^2 S/8 = 2 log n)
    Y = 2*np.sqrt(n)
    print(f"n={n}: ratio = {lemmaI(n, S, Y):.4f}   (n^3 S^2 = {n**3*S**2:.3g})")

print("== constants ==")
e = np.e
for k in (10, 100, 1000, 10**5):
    n = 4*k
    print(f"Hurwitz k={k}: w q n^2/(4e) / n^2 = {24*0.25/(4*e):.5f} = 3/(2e) = {3/(2*e):.5f};"
          f"  w D_k/(2e)/n^2 = {24*k*(2*k-1)/(2*e)/n**2:.5f}")
print("icosians: w q0 n^2/(4e)/n^2 =", 120*(1/8)/(4*e), " 15/(4e) =", 15/(4*e))

print("== (2b) EXACT closure of <i, (1+i+j+k)/2, (tau + tau^{-1} i + j)/2> in H(Q(sqrt5)) ==")
from fractions import Fraction as Fr
def kadd(x, y): return (x[0]+y[0], x[1]+y[1])
def kneg(x): return (-x[0], -x[1])
def kmul(x, y): return (x[0]*y[0] + 5*x[1]*y[1], x[0]*y[1] + x[1]*y[0])     # (a+b sqrt5)(c+d sqrt5)
def ksum(*xs):
    r = (Fr(0), Fr(0))
    for x in xs: r = kadd(r, x)
    return r
def Qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    m = kmul; n_ = kneg
    return (ksum(m(a0, b0), n_(m(a1, b1)), n_(m(a2, b2)), n_(m(a3, b3))),
            ksum(m(a0, b1), m(a1, b0), m(a2, b3), n_(m(a3, b2))),
            ksum(m(a0, b2), n_(m(a1, b3)), m(a2, b0), m(a3, b1)),
            ksum(m(a0, b3), m(a1, b2), n_(m(a2, b1)), m(a3, b0)))
Z0 = (Fr(0), Fr(0)); O1 = (Fr(1), Fr(0)); H1 = (Fr(1, 2), Fr(0))
one = (O1, Z0, Z0, Z0)
gi = (Z0, O1, Z0, Z0)
gs = (H1, H1, H1, H1)
gt = ((Fr(1, 4), Fr(1, 4)), (Fr(-1, 4), Fr(1, 4)), H1, Z0)       # tau/2, (tau-1)/2 = tau^{-1}/2, 1/2, 0
gensX = [gi, gs, gt]
grp = {one}; front = [one]
while front:
    new = []
    for a in front:
        for g in gensX:
            c = Qmul(a, g)
            if c not in grp:
                grp.add(c); new.append(c)
    front = new
    assert len(grp) <= 1000
minus_one = (kneg(O1), Z0, Z0, Z0)
norms_ok = all(ksum(*[kmul(c, c) for c in x]) == O1 for x in grp)
has_ijk = all(e in grp for e in [(Z0, O1, Z0, Z0), (Z0, Z0, O1, Z0), (Z0, Z0, Z0, O1)])
irr = any(any(c[1] != 0 for c in x) for x in grp)
print("exact order:", len(grp), "| -1 in G:", minus_one in grp, "| all reduced norms 1:", norms_ok,
      "| i,j,k in G:", has_ijk, "| some element outside Q<1,i,j,k>:", irr)
assert len(grp) == 120 and minus_one in grp and norms_ok and has_ijk and irr
