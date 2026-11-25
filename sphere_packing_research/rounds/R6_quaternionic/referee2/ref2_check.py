# Referee 2 independent checks (algebra): invariant forms, freeness, orbit-constancy, q(x).
import numpy as np, itertools
rng = np.random.default_rng(7)
def qmul(a, b):
    a0,a1,a2,a3 = a; b0,b1,b2,b3 = b
    return np.array([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0])
def Rmat(u):  # matrix of x -> x*u on H = R^4
    return np.column_stack([qmul(e, u) for e in np.eye(4)])
def closure(gens, s):
    # gens: list of tuples of s quaternions (one per real place); componentwise product
    key = lambda g: tuple(np.round(np.concatenate(g), 9) + 0.0)
    elems = {key(g): g for g in gens}
    frontier = list(gens)
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                p = tuple(qmul(g[v], h[v]) for v in range(s))
                kk = key(p)
                if kk not in elems:
                    elems[kk] = p; new.append(p)
        frontier = new
    return list(elems.values())
i_ = np.array([0,1.,0,0]); om = np.array([.5,.5,.5,.5])
s5 = np.sqrt(5); tau = (1+s5)/2; taub = (1-s5)/2
g3  = np.array([tau/2, (1/tau)/2, .5, 0]); g3b = np.array([taub/2, (1/taub)/2, .5, 0])
hur = closure([(i_,), (om,)], 1)
ico = closure([(i_, i_), (om, om), (g3, g3b)], 2)
print("orders:", len(hur), len(ico))
def analyse(G, s, k, name):
    n = 4*s*k
    def R(g):  # right mult on (H^k)^s, coordinates ordered place, index, quaternion coordinate
        M = np.zeros((n, n))
        for v in range(s):
            Rv = Rmat(g[v])
            for a in range(k):
                o = (v*k + a)*4
                M[o:o+4, o:o+4] = Rv
        return M
    Rs = [R(g) for g in G]
    assert all(np.allclose(M.T @ M, np.eye(n)) for M in Rs), "not orthogonal"
    # freeness: R_u - I invertible for u != 1
    mins = [np.linalg.svd(M - np.eye(n), compute_uv=False).min() for M in Rs]
    nontriv = sorted(mins)[1:]
    print(name, "k=%d" % k, "min sing. value of R_u - 1 over u != 1:", round(min(nontriv), 6), "(#trivial:", sum(m < 1e-9 for m in mins), ")")
    # R-span of the scalars
    span = np.linalg.matrix_rank(np.array([np.concatenate(g) for g in G]))
    # invariant symmetric forms = image of projector P on Sym_n ; compute rank of P
    idx = [(a, b) for a in range(n) for b in range(a, n)]
    def vec(S):  # orthonormal coordinates of Sym_n w.r.t. Tr(AB)
        return np.array([S[a, b] * (1 if a == b else np.sqrt(2)) for a, b in idx])
    basis = []
    for a, b in idx:
        E = np.zeros((n, n)); E[a, b] = E[b, a] = 1 if a == b else 1/np.sqrt(2); basis.append(E)
    Pm = np.column_stack([vec(sum(M @ E @ M.T for M in Rs) / len(Rs)) for E in basis])
    dimQ = int(round(np.trace(Pm)))
    Dk = s*(k + 4*k*(k-1)//2)
    # q(x) and orbit constancy
    x = rng.standard_normal(n)
    Px = sum(M @ np.outer(x, x) @ M.T for M in Rs) / len(Rs)
    q = np.sum(Px*Px) / np.dot(x, x)**2
    uv = [np.dot(x[v*4*k:(v+1)*4*k], x[v*4*k:(v+1)*4*k]) / np.dot(x, x) for v in range(s)]
    qf = sum(u*u for u in uv) / 4
    y = Rs[5] @ x
    Py = sum(M @ np.outer(y, y) @ M.T for M in Rs) / len(Rs)
    orb = {tuple(np.round(M @ x, 8)) for M in Rs}
    print("   R-span of G:", span, " dim Q (trace of projector):", dimQ, " formula D_k:", Dk,
          " q:", round(q, 12), " formula:", round(qf, 12), " |P(xu(x)xu)-P(x(x)x)|:", float(np.abs(Px-Py).max()), " orbit size:", len(orb))
analyse(hur, 1, 1, "Hurwitz"); analyse(hur, 1, 2, "Hurwitz"); analyse(hur, 1, 3, "Hurwitz")
analyse(ico, 2, 1, "icosian"); analyse(ico, 2, 2, "icosian")
# constants
e = np.e
print("3/(2e) =", 3/(2*e), " 15/(4e) =", 15/(4*e), " 1/(2e) =", 1/(2*e))
for (w, delta, nm) in [(2, 1, "trivial"), (24, 4, "Hurwitz"), (120, 8, "icosian")]:
    print(nm, "w q0/(4e) =", w/(delta*4*e), " ratio to 1/(2e):", (w/(delta*4*e))/(1/(2*e)))
