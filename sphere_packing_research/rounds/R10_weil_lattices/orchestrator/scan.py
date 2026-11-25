"""Round 13 (orchestrator): structured-sublattice scan for the SL_2(p) Weil lattices, p = 1 mod 4, (p+1)/2 prime.
For an odd-order subgroup C_m of the split torus the fixed lattice N_H(Lambda) has rank 2(p-1)/m and can be generated
from norm-map images of RANDOM orbit vectors (no need for the full orbit or the full lattice).  We assume D = r^2 p
(true for p = 5, 13, 37 where the full lattice was computed) and CHECK it through the calibration log R_U = 0 (Lemma B).
Then: exhaustive search for violators (Gamma < theta_p) over the whole t-circle, and identification of the unimodular
member L = (fixed lattice)/sqrt(m) at t = +-log(eps): roots and their s-values."""
import numpy as np, math, sys, time, json, collections
import proto as P, run37 as R
V = P.V

def fund_unit(p):
    """fundamental unit of Q(sqrt p), p = 1 mod 4 prime: smallest (a + b sqrt p)/2 > 1 with a^2 - p b^2 = -4 or +4."""
    b = 1
    while True:
        for sgn in (-4, 4):
            a2 = p * b * b + sgn
            a = math.isqrt(a2)
            if a * a == a2:
                return (a + b * math.sqrt(p)) / 2, a, b, sgn
        b += 1

def random_orbit_vectors(S, T, v0, nsamp, rng, steps=3, burn=25, wl=7):
    """independent short random words applied to v0 (keeps floating error at the 1e-13 level)."""
    h = S.shape[1]; p = 2 * h + 1
    Tph = np.array([[(a * x * x) % p for x in range(1, h + 1)] for a in range(1, S.shape[0] + 1)])   # exponents of T
    def Tpow(k):
        return np.exp(2j * np.pi * ((k * Tph) % p) / p)
    out = []
    for i in range(nsamp):
        u = np.broadcast_to(v0.astype(complex), (S.shape[0], h)).copy()
        for _ in range(wl):
            u = u * Tpow(int(rng.integers(0, p)))
            u = np.einsum('aij,aj->ai', S, u)
        out.append(u * Tpow(int(rng.integers(0, p))))
    return np.array(out)

def scan(p, ranks=(8, 24, 40), nsamp=300, tstep=0.1, seed=0, verbose=True):
    # nsamp raised from 160 to 300 after the round-13 verification: 160 samples generated only an index-2
    # sublattice of the fixed lattice at p = 541 (see rounds/R13_scan_verify).
    t0 = time.time(); rng = np.random.default_rng(seed)
    h = (p - 1) // 2; n = 2 * (p - 1); r = (p - 1) // 4
    S, T = P.build(p, -1)
    rel = P.check_relations(S, T, p)
    assert rel['S^2=-I'] < 1e-8 and rel['(ST)^3=-I'] < 1e-8, rel
    v0 = P.v0_vec(p)
    eps, ua, ub, usgn = fund_unit(p); L = math.log(eps)
    D = r * r * p; theta = (p - 1) ** 2 / (4 * p)
    O = random_orbit_vectors(S, T, v0, nsamp, rng)
    del S
    g0 = next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in set(R.P_factors(p - 1))))
    res = {'p': p, 'r': r, 'eps': eps, 'unit': [ua, ub, usgn], 'theta': theta, 'subs': []}
    print('p=%d r=%d (p+1)/2=%d  eps=(%d+%d sqrt p)/2, log eps=%.3f, theta=%.4f, build %.1fs' % (p, r, (p + 1) // 2, ua, ub, L, theta, time.time() - t0), flush=True)
    odd = [m for m in range(3, p, 2) if (p - 1) % m == 0 and (2 * (p - 1)) // m in ranks]
    for m in sorted(odd, reverse=True):
        rk_expect = 2 * (p - 1) // m
        a = pow(g0, (p - 1) // m, p); M = R.torus_matrix(p, a)
        N = np.eye(h); Mk = np.eye(h)
        for _ in range(m - 1):
            Mk = Mk @ M; N = N + Mk
        ON = np.einsum('ij,naj->nai', N, O)
        Yp, Ym = P.realify(ON, p)
        keep = np.where((np.abs(Yp).sum(1) + np.abs(Ym).sum(1)) > 1e-7)[0]
        try:
            BasH, rH = P.lattice_basis_exact(np.concatenate([Yp[keep], Ym[keep]], axis=1))
        except Exception as e:
            print('  C_%d: lattice_basis failed: %s' % (m, e)); continue
        dp = Yp.shape[1]
        GpH = BasH[:, :dp] @ BasH[:, :dp].T; GmH = BasH[:, dp:] @ BasH[:, dp:].T
        Bh = GpH + GmH
        integ = float(np.abs(Bh - np.round(Bh)).max())
        s2, ld = np.linalg.slogdet(np.round(Bh))
        logR = rH * math.log(n / rH) + (rH / 2) * math.log(D) - ld      # det(Lambda)^{rH/n} = D^{rH/2}
        info = {'m': m, 'rank': rH, 'integrality': integ, 'logR_assuming_D=r^2p': logR}
        print('  C_%d: rank %d (expected %d), integrality %.1e, log R_U (assuming D=r^2 p) = %.6f' % (m, rH, rk_expect, integ, logR), flush=True)
        if rH != rk_expect or abs(logR) > 1e-6 or integ > 1e-6:
            print('     WARNING: rank / calibration / integrality check failed -- the generated lattice is probably NOT the full fixed lattice (too few samples or precision loss); results below are unreliable', flush=True)
        # exhaustive violator search over the t-circle
        bound = 2 * math.sqrt(theta * D) * math.cosh(tstep / 2) * (1 + 1e-9)
        best = (1e99, None); cnt = 0; fail = 0; t1 = time.time()
        for t in np.arange(-L, L + 1e-9, tstep):
            Gt = math.exp(t) * GpH + math.exp(-t) * GmH
            try:
                Xs, vals, flag = V.short_vectors(Gt, bound, maxc=300000, maxnodes=300_000_000)
            except Exception as e:
                fail += 1; continue
            if flag == 2: fail += 1     # enumeration truncated by the cap: count as incomplete
            if len(Xs):
                mp = np.einsum('ij,jk,ik->i', Xs, GpH, Xs); mm = np.einsum('ij,jk,ik->i', Xs, GmH, Xs)
                gam = mp * mm / D; cnt += len(Xs)
                if gam.min() < best[0]: best = (float(gam.min()), float(t))
        ratio = best[0] / theta if best[1] is not None else float('nan')
        info.update({'min_gamma_over_theta': ratio, 'at_t': best[1], 'count': cnt, 'enum_fail': fail})
        print('     %s t-circle search: %d vectors below bound, min Gamma/theta = %.6f at t=%s, enum failures %d  [%.1fs]' % ('exhaustive' if fail == 0 else 'INCOMPLETE (not exhaustive)', cnt, ratio, best[1], fail, time.time() - t1), flush=True)
        # unimodular member: L = fixed lattice / sqrt(m) at t = +-log eps
        for tstar in (-L, L):
            Gu = (math.exp(tstar) * GpH + math.exp(-tstar) * GmH) / math.sqrt(D) / m
            dist = float(np.abs(Gu - np.round(Gu)).max())
            if dist < 1e-6:
                Gi = np.round(Gu); sg, ldu = np.linalg.slogdet(Gi)
                try:
                    Xs, vals, flag = V.short_vectors(Gi, 4.01, maxc=400000, maxnodes=400_000_000)
                except Exception as e:
                    print('     unimodular member: enumeration failed'); break
                mp = np.einsum('ij,jk,ik->i', Xs, GpH, Xs) * math.exp(tstar) / math.sqrt(D) / m
                mm = np.einsum('ij,jk,ik->i', Xs, GmH, Xs) * math.exp(-tstar) / math.sqrt(D) / m
                nrm = np.round(mp + mm).astype(int); s = np.round((mp - mm) * math.sqrt(p) / 2 * 1).astype(int)
                c = collections.Counter(zip(nrm.tolist(), np.abs(s).tolist()))
                desc = ', '.join('norm %d |s|=%d: %d pairs' % (k[0], k[1], v) for k, v in sorted(c.items()))
                nroots = int(2 * (nrm == 2).sum())
                info.update({'unimodular_logdet': float(ldu), 'roots': nroots, 'short': desc})
                print('     L = fixed/sqrt(%d) at t*=%+.3f: integral (%.0e), logdet %.1e, even %s, #roots = %d ; vectors of norm<=4: %s' % (m, tstar, dist, ldu, bool((np.diag(Gi) % 2 == 0).all()), nroots, desc[:400]), flush=True)
                break
        else:
            print('     no integral unimodular member found at t=+-log eps (dist %.2e)' % dist)
        res['subs'].append(info)
    print('  total %.1fs' % (time.time() - t0), flush=True)
    return res

if __name__ == '__main__':
    out = []
    for a in sys.argv[1:]:
        p = int(a)
        try:
            out.append(scan(p))
        except Exception as e:
            print('p=%d FAILED: %r' % (p, e), flush=True)
        json.dump(out, open('scan_results.json', 'w'), indent=1)
