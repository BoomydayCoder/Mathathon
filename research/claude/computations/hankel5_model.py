#!/usr/bin/env python3
"""
hankel5_model.py -- class-model study of GLOBAL (prime-independent) bases for Fauzan's zeta(5) determinant
at large K, inner range of primes.  Pure integer bookkeeping (numpy); nothing here is a proof.

Class model (Fauzan (4.2),(4.3); proved in Lean as `tauX_direct`): for p^2 > 5K, an entry
mu_X(D_N^6 f g / D_K) has Gauss valuation >= min( min_{c>=1}(ord_c f + ord_c g + 6 l_N(c) - l_K(c)) - 4,
2(ord_0 f + ord_0 g) + 12 m_N - 2 m_K + 1 ), where ord_c counts roots = +-c mod p.
   loc_p  = optimum of the model over all p-adic unimodular bases (water filling)  [= idealised gamma_p^in]
   H_p    = bound given by the weights lemma for the global HERMITE basis with multiplicities m_j
            f_{j,i} = prod_{j'!=j}(t+j'^2)^{m_j'} (t+j^2)^i (i<m_j), corrected by the generalised Vandermonde
            determinant prod_{j<j'} (j'^2-j^2)^{m_j m_j'} (squared, rows and columns).
   Newton pairs: tropical determinant (scipy) of the model matrix.
loss = sum_p (loc_p - bound_p) log p / K^2 over K/M < p <= K/3 (and p^2 > 5K).
Usage: python3 hankel5_model.py validate | scan K | newton K | search K nsteps iters seed | ext K
"""
import sys, math, time
import numpy as np

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def class_base(K, N, p):
    c = np.arange(0, (p - 1) // 2 + 1)
    lN = (N - c) // p + 1 + (N + c) // p
    lK = (K - c) // p + 1 + (K + c) // p
    b2 = 6 * lN - lK - 4                      # doubled half-weight 2*beta_c, c>=1
    b2[0] = 12 * (N // p) - 2 * (K // p) + 1  # zero class (step 4 per row in doubled units)
    return b2

def gamma_loc(K, N, h, p, b2=None):
    if b2 is None: b2 = class_base(K, N, p)
    nc = len(b2) - 1
    imax = (2 * h) // max(nc, 1) + 40
    i = np.arange(imax)
    w = (b2[1:, None] + 2 * i[None, :]).ravel()
    w = np.concatenate([w, b2[0] + 4 * i])
    w = np.partition(w, h - 1)[:h]
    assert w.max() < b2[1:].min() + 2 * (imax - 1)
    return int(w.sum())

def profile_from_blocks(blocks):
    """blocks: list of (A,B,height): m_j += height for A<j<=B."""
    J = max(B for A, B, hh in blocks)
    m = np.zeros(J + 1, dtype=np.int64)
    for A, B, hh in blocks:
        m[A + 1:B + 1] += hh
    return m

def hermite_bound(K, N, h, m, p, b2=None):
    """weights-lemma bound for v_p(Delta_K) in the global Hermite basis with multiplicities m[j], j>=1."""
    if b2 is None: b2 = class_base(K, N, p)
    nc = len(b2)
    js = np.nonzero(m)[0]
    mj = m[js]
    r = js % p
    cl = np.minimum(r, p - r)
    M = np.bincount(cl, weights=mj, minlength=nc).astype(np.int64)
    Q = np.bincount(cl, weights=mj * mj, minlength=nc).astype(np.int64)
    W2 = 2 * M + b2
    W2[0] = 4 * M[0] + b2[0]
    order = np.argsort(W2)
    w1, w2 = W2[order[0]], W2[order[1]]
    cap = np.where(cl == order[0], w2, w1)
    step = np.where(cl == 0, 4, 2)
    total = 0
    mm = int(mj.max())
    for d in range(1, mm + 1):           # rows with m_j - i = d
        sel = mj >= d
        own = W2[cl[sel]] - step[sel] * d
        total += int(np.minimum(own, cap[sel]).sum())
    vand = int(((M[1:] ** 2 - Q[1:]) // 2).sum() + (M[0] ** 2 - Q[0]))
    return total - 2 * vand

def inner_primes(K, M):
    return [p for p in primes_upto(K // 3) if p * M > K and p * p > 5 * K]

def loss_profile(K, N, h, m, P, locs=None, bases=None):
    tot = 0.0
    for k, p in enumerate(P):
        b2 = bases[k] if bases else None
        lo = locs[k] if locs else gamma_loc(K, N, h, p, b2)
        tot += (lo - hermite_bound(K, N, h, m, p, b2)) * math.log(p)
    return tot / K ** 2

def steps_profile(K, N, steps):
    """m = sum_i 1_[1,J_i] - 3*1_[1,N]."""
    m = profile_from_blocks([(0, J, 1) for J in steps] + [(0, N, -3)])
    assert m.min() >= 0, steps
    return m

# ------------------------------------------------------------ Newton pairs (model matrix + assignment)
def newton_trop(K, N, h, rows, cols, p):
    from scipy.optimize import linear_sum_assignment
    b2 = class_base(K, N, p)
    nc = len(b2)
    def cum(seq):
        Z = np.zeros((h, nc), dtype=np.int32)
        cur = np.zeros(nc, dtype=np.int32)
        for i in range(h):
            Z[i] = cur
            r = seq[i] % p
            cur[min(r, p - r)] += 1
        return Z
    Zr, Zc = cum(rows), cum(cols)
    V = 2 * (Zr[:, 0][:, None] + Zc[:, 0][None, :]) + int(b2[0])
    for c in range(1, nc):
        V = np.minimum(V, Zr[:, c][:, None] + Zc[:, c][None, :] + int(b2[c]))
    r, c = linear_sum_assignment(V)
    return int(V[r, c].sum())


# ------------------------------------------------------------ general designs (integer up/down steps, node 0, half-integers)
def design_nodes(K, N, h, up, down=(), half=(), m0=0):
    """integer multiplicities m_j = #{up>=j} - #{down>=j} - 3[j<=N] (j>=1), m_0 = m0 (node t=0),
       half-integer nodes (2k-1)/2 with multiplicity #{H in half: k<=H}.  Returns (m_int, m_half) or None."""
    J = max(list(up) + list(down) + [N])
    m = np.zeros(J + 1, dtype=np.int64)
    for x in up: m[1:x + 1] += 1
    for x in down: m[1:x + 1] -= 1
    m[1:N + 1] -= 3
    m[0] = m0
    mh = np.zeros(max(list(half) + [0]) + 1, dtype=np.int64)
    for x in half: mh[1:x + 1] += 1
    if m.min() < 0 or m.sum() + mh.sum() != h: return None
    return m, mh

def bound_design(K, N, h, nodes, p, b2):
    m, mh = nodes
    js = np.nonzero(m)[0]; mj = m[js]
    r = js % p; cl = np.minimum(r, p - r)
    if mh.sum():
        ks = np.nonzero(mh)[0]
        r2 = ((2 * ks - 1) * ((p + 1) // 2)) % p
        cl = np.concatenate([cl, np.minimum(r2, p - r2)]); mj = np.concatenate([mj, mh[ks]])
    nc = len(b2)
    M = np.bincount(cl, weights=mj, minlength=nc).astype(np.int64)
    Q = np.bincount(cl, weights=mj * mj, minlength=nc).astype(np.int64)
    W2 = 2 * M + b2; W2[0] = 4 * M[0] + b2[0]
    o = np.argpartition(W2, 1)[:2]; o = o[np.argsort(W2[o])]
    cap = np.where(cl == o[0], W2[o[1]], W2[o[0]])
    step = np.where(cl == 0, 4, 2)
    total = 0
    for d in range(1, int(mj.max()) + 1):
        sel = mj >= d
        total += int(np.minimum(W2[cl[sel]] - step[sel] * d, cap[sel]).sum())
    return total - 2 * int(((M[1:] ** 2 - Q[1:]) // 2).sum() + (M[0] ** 2 - Q[0]))

class Ctx:
    def __init__(self, K, M, hi=3.0):
        self.K = K; self.N = 3 * K // 40; self.h = K - self.N; self.M = M
        self.P = [p for p in primes_upto(int(K / hi)) if p * M > K and p * p > 5 * K]
        self.b = [class_base(K, self.N, p) for p in self.P]
        self.loc = [gamma_loc(K, self.N, self.h, p, b) for p, b in zip(self.P, self.b)]
    def loss(self, nodes):
        return sum((lo - bound_design(self.K, self.N, self.h, nodes, p, b)) * math.log(p)
                   for p, b, lo in zip(self.P, self.b, self.loc)) / self.K ** 2

def scale_design(K, th_up, th_down=(), th_half=(), m0=0):
    N = 3 * K // 40; h = K - N
    up = [int(round(t * K)) for t in th_up]; down = [int(round(t * K)) for t in th_down]
    half = [int(round(t * K)) for t in th_half]
    up[-1] += h + 3 * N - (sum(up) - sum(down) + sum(half) + m0)
    return design_nodes(K, N, h, up, down, half, m0), (up, down, half, m0)

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == 'validate':
        K, N, h = 120, 9, 111
        tgt = {'herm3': {29: -303, 31: -255, 37: -346}, 'herm2': {29: -322, 31: -318, 37: -327},
               'lagr': {29: -374, 31: -345, 37: -438}}
        prof = {'herm3': [(N, N + 37, 3)], 'herm2': [(N, N + 55, 2), (N + 55, N + 56, 1)], 'lagr': [(N, K, 1)]}
        for nm in prof:
            m = profile_from_blocks(prof[nm]); assert m.sum() == h
            print(nm, {p: (hermite_bound(K, N, h, m, p), tgt[nm][p]) for p in (29, 31, 37)},
                  ' (fast formula, tropical of model matrix from hankel5_exact)')
        print('loc', {p: gamma_loc(K, N, h, p) for p in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)},
              ' expected 385,228,51,-6,-113,-218,-255,-300,-348,-372,-385')
        rows = list(range(N + 1, K + 1)); cols = list(range(K, N, -1))
        print('tail  model trop', {p: newton_trop(K, N, h, rows, cols, p) for p in (29, 31, 37)}, ' expected -484,-492,-498')
        print('symlow model trop', {p: newton_trop(K, N, h, rows, rows, p) for p in (29, 31, 37)}, ' expected -372,-444,-384')
    elif mode == 'scan':
        K = int(sys.argv[2]); N = 3 * K // 40; h = K - N
        for M in (10, 20, 40):
            P = inner_primes(K, M)
            if not P: continue
            bases = [class_base(K, N, p) for p in P]
            locs = [gamma_loc(K, N, h, p, b) for p, b in zip(P, bases)]
            cost = -sum(l * math.log(p) for l, p in zip(locs, P)) / K ** 2
            print(f"K={K} M={M}: {len(P)} primes in ({max(K/M,(5*K)**.5):.0f},{K//3}];  -sum loc_p log p/K^2 = {cost:.4f}")
            cands = {
                'lagr  m=1 on (N,K]': [(N, K, 1)],
                'herm2 m=2 on (N,N+h/2]': [(N, N + h // 2, 2)] + ([(N + h // 2, N + h // 2 + 1, h % 2)] if h % 2 else []),
                'herm3 m=3 on (N,N+h/3]': [(N, N + h // 3, 3)] + ([(N + h // 3, N + h // 3 + 1, h % 3)] if h % 3 else []),
                'herm4 m=4 on (N,N+h/4]': [(N, N + h // 4, 4)] + ([(N + h // 4, N + h // 4 + 1, h % 4)] if h % 4 else []),
                'm=1 on [1,h] (zeros of W re-used)': [(0, h, 1)],
                'm=3 on (0,h/3]': [(0, h // 3, 3)] + ([(h // 3, h // 3 + 1, h % 3)] if h % 3 else []),
            }
            for nm, bl in cands.items():
                m = profile_from_blocks(bl); assert m.sum() == h, (nm, m.sum())
                print(f"    {nm:36s} loss {loss_profile(K, N, h, m, P, locs, bases):.4f}")
    elif mode == 'newton':
        K = int(sys.argv[2]); N = 3 * K // 40; h = K - N
        low = list(range(N + 1, K + 1)); top = list(range(K, N, -1))
        low2 = [j for j in range(N + 1, N + h // 2 + 2) for _ in (0, 1)][:h]
        one = list(range(1, h + 1))
        pairs = {'tail': (low, top), 'symlow': (low, low), 'low2': (low2, low2), 's3': (one, top)}
        for M in (10, 20):
            P = inner_primes(K, M)
            locs = [gamma_loc(K, N, h, p) for p in P]
            for nm, (r, c) in pairs.items():
                t0 = time.time()
                L = sum((lo - newton_trop(K, N, h, r, c, p)) * math.log(p) for lo, p in zip(locs, P)) / K ** 2
                print(f"K={K} M={M} Newton pair {nm:7s}: model loss {L:.4f}  ({time.time()-t0:.0f}s)", flush=True)
    elif mode == 'search':
        K = int(sys.argv[2]); nst = int(sys.argv[3]); iters = int(sys.argv[4]); seed = int(sys.argv[5])
        M = int(sys.argv[6]) if len(sys.argv) > 6 else 20
        N = 3 * K // 40; h = K - N
        P = inner_primes(K, M)
        bases = [class_base(K, N, p) for p in P]
        locs = [gamma_loc(K, N, h, p, b) for p, b in zip(P, bases)]
        rng = np.random.default_rng(seed)
        S = K + 2 * N                       # sum of the step positions
        def rand_steps():
            while True:
                x = np.sort(rng.dirichlet(np.ones(nst)))  # positions as fractions of S
                st = [int(round(v * S)) for v in x]
                st[-1] += S - sum(st)
                st.sort()
                if st[-3] >= N and st[0] >= 0 and st[-1] <= 3 * K: return st
        def f(st):
            return loss_profile(K, N, h, steps_profile(K, N, st), P, locs, bases)
        best, bv = None, 1e9
        for it in range(iters):
            if best is None or rng.random() < 0.3:
                st = rand_steps()
            else:
                st = list(best)
                i, j = rng.choice(nst, 2, replace=False)
                d = int(rng.integers(1, max(2, K // 20)))
                st[i] += d; st[j] -= d; st.sort()
                if st[0] < 0 or st[-3] < N: continue
            v = f(st)
            if v < bv: best, bv = st, v
        print(f"K={K} M={M} steps={nst} iters={iters} seed={seed}: best loss {bv:.4f} at steps {best} (N={N}, K={K}, sum={sum(best)})")
    elif mode == 'ext':
        # by-product: is the water-filling (inner) bound better than (4.14) for K/3 < p < K/2 ?
        K = int(sys.argv[2]); N = 3 * K // 40; h = K - N
        def gamma_out(p):
            v = K - p * (K // p); u = max(0, N + v - p + 1); tp = min(N, v) + u; rp = max(0, K + 4 * N - 2 * p + 2)
            if K < 2 * p: return -7 * (K - p) + 6 * tp - 1 - min(rp, p - 1 - N + u)
            return -7 * (K - p) + 3 + 12 * N + 5 * tp - min(rp, p + u)
        g = 0.0; cnt = 0
        for p in primes_upto(K):
            if 3 * p > K and p * p > 5 * K:
                lo, go = gamma_loc(K, N, h, p), gamma_out(p)
                if lo > go:
                    g += (lo - go) * math.log(p); cnt += 1; last = K / p
                    first = K / p if cnt == 1 else first
        print(f"K={K}: sum over p>K/3 of (loc_p - out_p)_+ log p / K^2 = {g/K**2:.4f}  ({cnt} primes, K/p in [{last:.3f},{first:.3f}])")
    elif mode == 'search2':
        # hill climbing: search2 K M n_up n_down n_half use_m0 iters seed
        K, M, nu, nd, nh, um0, iters, seed = [int(x) for x in sys.argv[2:10]]
        C = Ctx(K, M); N, h = C.N, C.h
        rng = np.random.default_rng(seed)
        sg = [1] * nu + [-1] * nd + [1] * nh + [1]
        tot = h + 3 * N
        def build(x):
            if min(x) < 0: return None
            return design_nodes(K, N, h, x[:nu], x[nu:nu + nd], x[nu + nd:nu + nd + nh], x[-1])
        def start():
            while True:
                x = [int(v) for v in rng.integers(N, K, nu)] + [int(v) for v in rng.integers(0, K // 2, nd)] + \
                    [int(v) for v in rng.integers(0, K // 2, nh)] + [int(rng.integers(0, 6)) if um0 else 0]
                x[0] += tot - sum(a * b for a, b in zip(sg, x))
                nd_ = build(x)
                if nd_ is not None: return x, C.loss(nd_)
        best, bv = None, 1e9
        cur, cv = start()
        scale = K // 8; stall = 0
        for it in range(iters):
            x = list(cur)
            na = len(x) - (0 if um0 else 1)
            a, b = rng.choice(na, 2, replace=False)
            d = int(rng.integers(1, max(2, scale))) * (1 if rng.random() < .5 else -1)
            if a == len(x) - 1 or b == len(x) - 1: d = 1 if d > 0 else -1
            x[a] += d; x[b] -= d * sg[a] * sg[b]
            nd_ = build(x)
            if nd_ is None: continue
            v = C.loss(nd_)
            if v < cv:
                cur, cv, stall = x, v, 0
            else:
                stall += 1
                if stall > 40:
                    scale = max(2, scale // 2); stall = 0
                    if scale == 2 and rng.random() < 0.15:
                        if cv < bv: best, bv = cur, cv
                        cur, cv = start(); scale = K // 8
            if cv < bv: best, bv = list(cur), cv
        print(f"K={K} M={M} up={nu} down={nd} half={nh} m0={um0} iters={iters} seed={seed}: best loss {bv:.4f}"
              f"  up={[round(v/K,4) for v in best[:nu]]} down={[round(v/K,4) for v in best[nu:nu+nd]]}"
              f" half={[round(v/K,4) for v in best[nu+nd:nu+nd+nh]]} m0={best[-1]}  raw={best}")
    elif mode == 'eval':
        # eval K "th_up" "th_down" "th_half" m0  [hi]   -> losses for several M ; hi: upper end p<=K/hi
        K = int(sys.argv[2])
        f = lambda t: [float(x) for x in t.split(',') if x]
        thu, thd, thh = f(sys.argv[3]), f(sys.argv[4]), f(sys.argv[5]); m0 = int(sys.argv[6])
        hi = float(sys.argv[7]) if len(sys.argv) > 7 else 3.0
        nodes, raw = scale_design(K, thu, thd, thh, m0)
        assert nodes is not None, raw
        out = []
        Ms = [int(x) for x in sys.argv[8].split(',')] if len(sys.argv) > 8 else (6, 10, 20, 40, 80, 200, 1000)
        for M in Ms:
            if K / M < (5 * K) ** .5 and M > 6 and K / (M / 2) < (5 * K) ** .5: break
            C = Ctx(K, M, hi)
            out.append(f"M={M}:{C.loss(nodes):.4f}({len(C.P)}pr)")
        print(f"K={K} design up={thu} down={thd} half={thh} m0={m0} raw={raw} p<=K/{hi}: loss vs loc  " + "  ".join(out), flush=True)
    elif mode == 'fz':
        # consistency: sum_{K/20<p<=K/3} -(loc_p + v_p(S_K)) log p / K^2  vs Fauzan (5.18) = 0.042788...
        K = int(sys.argv[2]); N = 3 * K // 40; h = K - N
        tot = 0.0
        i = np.arange(1, h)
        for p in primes_upto(K // 3):
            if 20 * p > K and p * p > 5 * K:
                vS = 2 * h * (K // p) - 12 * h * (N // p) - 2 * int((2 * i // p).sum())
                tot -= (gamma_loc(K, N, h, p) + vS) * math.log(p)
        print(f"K={K}: sum_(K/20<p<=K/3) -(loc_p+v_p(S_K)) log p/K^2 = {tot/K**2:.5f}   [Fauzan (5.18): 0.04279]")
    elif mode == 'gammacheck':
        # loc_p / p  versus Fauzan's limiting function Gamma(x), formula (5.4)
        K = int(sys.argv[2]); N = 3 * K // 40; h = K - N
        al, Hc = 3 / 40, 23 / 20
        P = primes_upto(K // 3)
        for xt in (3.3, 4.7, 6.1, 8.9, 13.3, 19.7, 31.0):
            p = min(P, key=lambda q: abs(K / q - xt)); x = K / p
            z = (np.arange(1, (p - 1) // 2 + 1)) / p
            ell = lambda y: np.floor(y - z) + np.floor(y + z) + 1
            T = math.floor(2 * Hc * x); s_ = Hc * x - T / 2; q = math.floor(2 * x); npl = (2 * x - q) / 2
            b = 3 * ell(al * x)
            G = ((T - b) * (T + b - ell(x) - 5)).sum() / p + s_ * (2 * T - q - 5) + max(0.0, s_ - npl)
            print(f"   p={p} x={x:.3f}: loc_p/p = {gamma_loc(K, N, h, p)/p:.4f}   Gamma(x) = {G:.4f}")
