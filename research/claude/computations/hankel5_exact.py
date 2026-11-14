#!/usr/bin/env python3
"""
hankel5_exact.py -- can ONE prime-independent pair of polynomial bases + entrywise valuations + the
weights lemma replace Fauzan's p-dependent bases?  Exact experiment (python-flint, gmpy2).

Family (Fauzan, s=5):  G = [mu_X(D_N^6 t^{i+j}/D_K)]_{i,j<h},  K=40n, N=3n, h=37n,  Delta_K = det G.
For a pair of Q-bases (f_i), (g_i') of Q[t]_{<h}:  C = [mu_X(D_N^6 f_i g_i'/D_K)],  det C = det F det G' Delta_K.
For every prime p we compute
   true_p  = Gauss valuation of Delta_K (exact),
   trop_p  = min over permutations of sum_i v_p(c_{i,sigma(i)})  -  v_p(det F) - v_p(det G')
           = the best lower bound for v_p(Delta_K) obtainable from entrywise valuations in this basis
             (LP duality: = best "weights lemma" bound),
   cls_p   = the same with v_p(c_ii') replaced by the "class model" lower bound
             min( min_{c>=1}(Z_c - P_c) - 4 , 2(Z_0 - P_0) + 1 )   [Fauzan (4.2),(4.3); Lean `tauX_direct`]
   loc_p   = class-model optimum over all p-adic bases (water filling = idealised Fauzan gamma_p^in),
   out_p   = Fauzan's gamma_p^out, formula (4.14), when his conditions (4.9) hold.
A basis element is a multiset of nodes (u,v): the polynomial prod (v^2 t + u^2)  (root t = -(u/v)^2).
No claims are made by this file; it only computes.
Usage: python3 hankel5_exact.py n [basis,basis,...]
"""
import sys, time, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import primes_upto, moments_F1, hankel, det_poly_affine, gauss_val, vp
import numpy as np
import gmpy2
from scipy.optimize import linear_sum_assignment
from flint import fmpq, fmpz, fmpq_mat, fmpz_poly

BIG = 10 ** 6


def poly_from_nodes(nodes):
    P = fmpz_poly([1])
    for (u, v) in nodes:
        P = P * fmpz_poly([u * u, v * v])
    return P


# ---------------------------------------------------------------- bases (lists of node multisets)
def newton(seq, h):
    return [list(seq[:i]) for i in range(h)]


def hermite(mult, h):
    """mult: dict node->(multiplicity); elements f_{j,i} = prod_{j'!=j}(t+j'^2)^{m_j'} (t+j^2)^i, i<m_j."""
    assert sum(mult.values()) == h
    out = []
    for j, m in mult.items():
        rest = [nd for jj, mm in mult.items() if jj != j for nd in [jj] * mm]
        for i in range(m):
            out.append(rest + [j] * i)
    return out


def bases(K, N, h):
    I = lambda j: (j, 1)
    low = [I(j) for j in range(N + 1, N + h + 1)]
    top = [I(j) for j in range(K, K - h, -1)]
    one = [I(j) for j in range(1, h + 1)]
    fq = [I(j) for j in range(0, h)]
    half = [(2 * j - 1, 2) for j in range(1, h + 1)]
    il_r, il_c = [], []
    for k in range(h):
        il_r += [I(N + 1 + k), I(K - k)]
        il_c += [I(K - k), I(N + 1 + k)]
    # low nodes with the zeros of W re-inserted late: N+1..K then nothing (same as low); variant: 2 passes
    low2 = [I(j) for j in range(N + 1, N + (h + 1) // 2 + 1) for _ in (0, 1)]
    lag = {I(j): 1 for j in range(N + 1, K + 1)}
    q, r = divmod(h, 3)
    h3 = {I(N + 1 + k): 3 for k in range(q)}
    if r:
        h3[I(N + 1 + q)] = r
    q2, r2 = divmod(h, 2)
    h2 = {I(N + 1 + k): 2 for k in range(q2)}
    if r2:
        h2[I(N + 1 + q2)] = r2
    B = {
        'mono': (newton([I(0)] * h, h), newton([I(0)] * h, h)),
        'tail': (newton(low, h), newton(top, h)),          # rows N+1.., columns K,K-1,..   (previous run)
        'symlow': (newton(low, h), newton(low, h)),
        'symtop': (newton(top, h), newton(top, h)),
        's3': (newton(one, h), newton(top, h)),            # the zeta(3) basis of hankel-zeta3.md
        'fq': (newton(fq, h), newton(fq, h)),              # Fauzan's q_i (3.11) up to row/column scalars
        'half': (newton(half, h), newton(half, h)),        # nodes 1/2,3/2,..
        'halftop': (newton(half, h), newton(top, h)),
        'interl': (newton(il_r, h), newton(il_c, h)),
        'low2': (newton(low2, h), newton(low2, h)),        # double nodes N+1,N+1,N+2,N+2,..
        'lagr': (hermite(lag, h), hermite(lag, h)),        # Lagrange basis on the poles
        'herm2': (hermite(h2, h), hermite(h2, h)),
        'herm3': (hermite(h3, h), hermite(h3, h)),
        'lagtail': (hermite(lag, h), newton(top, h)),
    }
    return B


# ---------------------------------------------------------------- class model
def ell(A, c, p):
    """#{1<=j<=A : j = +-c mod p}, 1<=c<=(p-1)/2."""
    return (A - c) // p + 1 + (A + c) // p


def class_counts(basis, p):
    m = (p - 1) // 2
    Z = np.zeros((len(basis), m + 1), dtype=np.int64)
    for i, nodes in enumerate(basis):
        for (u, v) in nodes:
            r = (u * pow(v, -1, p)) % p
            Z[i, min(r, p - r)] += 1
    return Z


def class_base(K, N, p):
    m = (p - 1) // 2
    base = np.array([0] + [6 * ell(N, c, p) - ell(K, c, p) for c in range(1, m + 1)], dtype=np.int64)
    base0 = 12 * (N // p) - 2 * (K // p) + 1
    return base, base0


def class_model_matrix(rows, cols, K, N, p):
    Zr, Zc = class_counts(rows, p), class_counts(cols, p)
    base, base0 = class_base(K, N, p)
    V = 2 * (Zr[:, 0][:, None] + Zc[:, 0][None, :]) + base0
    for c in range(1, (p - 1) // 2 + 1):
        V = np.minimum(V, Zr[:, c][:, None] + Zc[:, c][None, :] + base[c] - 4)
    return V


def gamma_loc(K, N, h, p):
    """class-model optimum over p-adic bases: 2*(sum of the h smallest half-weights)."""
    base, base0 = class_base(K, N, p)
    w = []
    for c in range(1, (p - 1) // 2 + 1):
        w += [(base[c] - 4) / 2 + i for i in range(h)]
    w += [base0 / 2 + 2 * i for i in range(h)]
    w.sort()
    return int(round(2 * sum(w[:h])))


def gamma_out(K, N, p):
    """Fauzan (4.14); None if (4.9) fails."""
    if not (p >= 7 and p <= K < 3 * p and p * p > 2 * K and 2 * N < p and 5 * N <= 2 * p - 2):
        return None
    v = K - p * (K // p)
    u = max(0, N + v - p + 1)
    tp = min(N, v) + u
    rp = max(0, K + 4 * N - 2 * p + 2)
    if K < 2 * p:
        return -7 * (K - p) + 6 * tp - 1 - min(rp, p - 1 - N + u)
    return -7 * (K - p) + 3 + 12 * N + 5 * tp - min(rp, p + u)


def trop(V):
    r, c = linear_sum_assignment(V)
    return int(round(V[r, c].sum()))


def vq(q, p):
    """valuation of fmpq via gmpy2."""
    num = gmpy2.mpz(int(q.p))
    if num == 0:
        return BIG
    a = gmpy2.remove(num, p)[1]
    if a:
        return int(a)
    return -int(gmpy2.remove(gmpy2.mpz(int(q.q)), p)[1])


if __name__ == "__main__":
    n = int(sys.argv[1])
    K, N, a = 40 * n, 3 * n, 6
    h = K - N
    t0 = time.time()
    al, be = moments_F1(5, K, N, a, h)
    Ha, Hb = hankel(al, h), hankel(be, h)
    cs = det_poly_affine(Ha, Hb)
    print(f"# s=5 K={K} N={N} h={h}  true det computed ({time.time()-t0:.0f}s)", flush=True)
    P = [p for p in primes_upto(K + 20)]
    true = {p: gauss_val(cs, p) for p in P}
    B = bases(K, N, h)
    names = sys.argv[2].split(',') if len(sys.argv) > 2 else list(B)
    # optional staircase Hermite profiles:  stair:J1:J2:J3  ->  m_j = #{i: J_i >= j} - 3*[j<=N]  (rows = columns)
    # optional asymmetric pair              pair:J1:J2:J3/J1':J2':J3'  (rows / columns)
    def stair(spec):
        st = [int(x) for x in spec.split(':')]
        mult = {}
        for j in range(1, max(st) + 1):
            mj = sum(1 for J in st if J >= j) - (3 if j <= N else 0)
            assert mj >= 0
            if mj:
                mult[(j, 1)] = mj
        return hermite(mult, h)
    for nm in names:
        if nm.startswith('stair:'):
            B[nm] = (stair(nm[6:]), stair(nm[6:]))
        if nm.startswith('pair:'):
            a_, b_ = nm[5:].split('/')
            B[nm] = (stair(a_), stair(b_))
    res, resm, viol = {}, {}, {}
    for name in names:
        rows, cols = B[name]
        Fp = [poly_from_nodes(nd) for nd in rows]
        Gp = [poly_from_nodes(nd) for nd in cols]

        def cm(polys):
            ent = []
            for f in polys:
                c = [int(x) for x in f.coeffs()]
                assert len(c) <= h
                ent += c + [0] * (h - len(c))
            return fmpq_mat(h, h, ent)
        F, G = cm(Fp), cm(Gp)
        dF, dG = F.det(), G.det()
        assert dF != 0 and dG != 0
        Gt = G.transpose()
        Ca, Cb = F * Ha * Gt, F * Hb * Gt
        ea = [(gmpy2.mpz(int(x.p)), gmpy2.mpz(int(x.q))) for x in Ca.entries()]
        eb = [(gmpy2.mpz(int(x.p)), gmpy2.mpz(int(x.q))) for x in Cb.entries()]
        res[name], resm[name], viol[name] = {}, {}, 0
        for p in P:
            V = np.empty(h * h)
            for k in range(h * h):
                v = BIG
                for (nu, de) in (ea[k], eb[k]):
                    if nu != 0:
                        x = int(gmpy2.remove(nu, p)[1])
                        if x == 0:
                            x = -int(gmpy2.remove(de, p)[1])
                        v = min(v, x)
                V[k] = v
            V = V.reshape(h, h)
            corr = vq(dF, p) + vq(dG, p)
            res[name][p] = trop(V) - corr
            if p >= 3:
                Vm = class_model_matrix(rows, cols, K, N, p)
                resm[name][p] = trop(Vm.astype(float)) - corr
                if p * p > 5 * K:
                    viol[name] += int((V < Vm - 1e-9).sum())
        print(f"# basis {name} done ({time.time()-t0:.0f}s)", flush=True)

    rng = lambda p: 'small' if p * p <= 5 * K else ('inner' if 3 * p <= K else ('outer' if p <= K else 'beyond'))
    print(f"\n== exact tropical bounds for v_p(Delta_K), K={K}  (true = exact valuation; out=(4.14); loc=class-model optimum)")
    print("   p  rng   | true |  out |  loc | " + " ".join(f"{nm[:7]:>7s}" for nm in names))
    for p in P:
        go = gamma_out(K, N, p)
        gl = gamma_loc(K, N, h, p) if p >= 3 else None
        f = lambda x: '   . ' if x is None else f"{x:5d}"
        print(f" {p:4d} {rng(p):6s}| {true[p]:5d}|{f(go)} |{f(gl)} | " + " ".join(f"{res[nm][p]:7d}" for nm in names))
    print(f"\n== class-model tropical bounds (same bases), K={K}")
    print("   p  rng   | true |  loc | " + " ".join(f"{nm[:7]:>7s}" for nm in names))
    for p in P:
        if p < 3:
            continue
        print(f" {p:4d} {rng(p):6s}| {true[p]:5d}|{gamma_loc(K, N, h, p):5d} | " + " ".join(f"{resm[nm][p]:7d}" for nm in names))
    print(f"\n== loss = sum_p (true_p - bound_p) log p / K^2 per range, K={K}")
    print("   range  | true cost |   out   |   loc   | " + " ".join(f"{nm[:7]:>7s}" for nm in names) + "   (exact tropical)")
    for R in ('small', 'inner', 'outer', 'beyond', 'all'):
        PR = [p for p in P if R == 'all' or rng(p) == R]
        tc = -sum(true[p] * math.log(p) for p in PR) / K ** 2
        lo = [p for p in PR if gamma_out(K, N, p) is not None]
        so = sum((true[p] - gamma_out(K, N, p)) * math.log(p) for p in lo) / K ** 2 if lo else float('nan')
        sl = sum((true[p] - gamma_loc(K, N, h, p)) * math.log(p) for p in PR if p >= 3) / K ** 2
        print(f"   {R:6s} | {tc:9.4f} | {so:7.4f} | {sl:7.4f} | " + " ".join(
            f"{sum((true[p]-res[nm][p])*math.log(p) for p in PR)/K**2:7.4f}" for nm in names))
    print("   (class-model tropical, p>=3)")
    for R in ('small', 'inner', 'outer', 'all'):
        PR = [p for p in P if p >= 3 and (R == 'all' or rng(p) == R)]
        print(f"   {R:6s} |           |         |         | " + " ".join(
            f"{sum((true[p]-resm[nm][p])*math.log(p) for p in PR)/K**2:7.4f}" for nm in names))
    print("\n== entries with exact valuation < class-model bound, primes with p^2>5K: " +
          ", ".join(f"{nm}:{viol[nm]}" for nm in names))
    print(f"# total time {time.time()-t0:.0f}s")
