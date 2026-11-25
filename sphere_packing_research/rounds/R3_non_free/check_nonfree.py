#!/usr/bin/env python3
"""R3_non_free: numerical sanity checks for notes.md.  Light (<10 s)."""
import itertools, math
import numpy as np
from mpmath import mp, mpf, loggamma, log, pi, exp, sqrt, gamma

mp.dps = 30

def logvol(m, R):  # log volume of m-ball of radius R
    return (mpf(m) / 2) * log(pi) + m * log(R) - loggamma(mpf(m) / 2 + 1)

print("== (1) asymptotic form of the lower-rank term T_U = Vol_{fn}(R)/c_U^k ==")
# W: dim_Q = D = r*d, c_W=1 wlog, U of Q-dim s*d, tau2 given; n = D*k
for (D, sdim, k, tau2, V) in [(4, 2, 50, 2.0, 1e4), (6, 2, 40, 3.0, 1e3), (8, 4, 30, 1.7, 1e2), (4, 2, 200, 2.0, 1e6)]:
    n = D * k; f = mpf(sdim) / D; fn = sdim * k
    # R from V = Vol_n(R)/c_W^k, c_W = 1
    logR = (log(V) + loggamma(mpf(n) / 2 + 1) - (mpf(n) / 2) * log(pi)) / n
    R = exp(logR)
    cU = mpf(tau2) ** (mpf(sdim) / 2)  # c_U^{2/sdim} = tau2 * c_W^{2/D}
    exact = exp(logvol(fn, R) - k * log(cU))
    approx = mpf(V) ** f * (pi * n) ** (f / 2) / sqrt(pi * f * n) * ((1 / f) / tau2) ** (mpf(fn) / 2)
    print(f"  D={D} dimU={sdim} k={k} n={n} tau2={tau2} V={V:g}: exact T_U={float(exact):.6g} approx={float(approx):.6g} ratio={float(exact/approx):.5f}")

print("== (2) examples: tau_U^2 = c_U^{2/dimU}/c_W^{2/dimW} vs 1/f = dimW/dimU ==")
def covol(B):
    B = np.array(B, dtype=float)
    return math.sqrt(abs(np.linalg.det(B @ B.T)))

def report(name, cU, dU, cW, dW):
    tau2 = cU ** (2 / dU) / cW ** (2 / dW)
    need = dW / dU
    print(f"  {name}: c_U={cU:.5g} dimU={dU} c_W={cW:.5g} dimW={dW}  tau^2={tau2:.5f}  need>={need:.3f}  ratio f*tau^2={tau2/need:.4f}  {'PASS' if tau2 >= need - 1e-9 else 'FAIL'}")

# Z^2, D_8, K=Q r=2
report("Z^2 line (1,0)", 1.0, 1, 1.0, 2)
# A_2, D_12
A2 = [[1, 0], [-0.5, math.sqrt(3) / 2]]
A2 = (np.array(A2) * math.sqrt(2)).tolist()
report("A_2 line (min vector)", math.sqrt(2), 1, covol(A2), 2)
# D_4 over Q (2_+^{1+4}, r=4)
D4 = [[1, 1, 0, 0], [1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]]
report("D_4/Q line (min vector)", math.sqrt(2), 1, covol(D4), 4)
D4U = [[1, 1, 0, 0], [1, -1, 0, 0]]  # D_4 cap fixed plane of diag(1,1,-1,-1)
report("D_4/Q fixed plane of diag(1,1,-1,-1)", covol(D4U), 2, covol(D4), 4)
print("     index [D4 : (D4 cap U)+(D4 cap U')] =", covol(D4U) * covol(D4U) / covol(D4), " (max possible 2^a = 4)")
# E_8 over Q r=8 : min norm 2, covol 1
report("E_8/Q line (min vector)", math.sqrt(2), 1, 1.0, 8)
# D_4 as Hermitian Z[sqrt(-2)]-lattice of rank 2: K-line = real 2-dim lattice w*O_K, covol = sqrt(|d_K|)/2 * |w|^2 ... check directly
# realize: Hurwitz quaternions scaled; s = (i+j) acts by left mult, s^2=-2.  w=1+i (norm 2 in D4 = {x in Z^4, sum even} model)
def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2, a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0)
s = (0, 1, 1, 0)
best = None
for w in itertools.product(range(-2, 3), repeat=4):
    if sum(w) % 2 or not any(w):
        continue
    sw = qmul(s, w)
    assert sum(sw) % 2 == 0  # D_4 stable under s
    c = covol([w, sw])  # w*Z[sqrt-2] : is it saturated? index-check below
    if best is None or c < best[0] - 1e-12:
        best = (c, w)
report(f"D_4 over Z[sqrt-2] (r=2,d=2), K-line through w={best[1]}", best[0], 2, covol(D4), 4)

print("== (3) entropy inequality  m*h(a/m) <= 2*min(a,b)*log q ,  m=a+b(q-1) ==")
def h(x):
    return 0.0 if x <= 0 or x >= 1 else -x * math.log(x) - (1 - x) * math.log(1 - x)
for q in [2, 3, 5, 7, 11]:
    ok = []
    Bq = 2000
    for a in range(1, 12 * Bq):
        m = a + Bq * (q - 1)
        if m * h(a / m) <= 2 * min(a, Bq) * math.log(q) + 1e-9:
            ok.append(a / Bq)
    if ok:
        lo, hi = min(ok), max(ok)
        rho = lambda t: (t - 1) / (t + q - 1)
        print(f"  q={q}: allowed a/b in [{lo:.4f}, {hi:.4f}]  (e/q={math.e/q:.4f}); normalized character (a-b)/m in [{rho(lo):.4f}, {rho(hi):.4f}], free value -1/(q-1)={-1/(q-1):.4f}")
    else:
        print(f"  q={q}: no admissible a/b")

print("== (4) Heisenberg p^{1+2} on O_K^p, K=Q(zeta_p): tau_line^2=1 (need p). Upper bound for ANY O_K-lattice of rank r=p ==")
def blichfeldt(m):  # gamma_m <= (2/pi) Gamma(2+m/2)^{2/m}
    return float((2 / pi) * exp(2 * loggamma(mpf(2) + mpf(m) / 2) / m))
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
    d = p - 1; r = p; m = r * d
    rd = p ** ((p - 2) / (p - 1))
    ub_bl = blichfeldt(m) * rd / d
    ub_kl = 1.744 * m / (2 * math.pi * math.e) * rd / d
    mh = m / (2 * math.pi * math.e) * rd / d
    print(f"  p={p}: dim={m} rd_K={rd:.3f} need tau^2>={r}; Blichfeldt UB {ub_bl:.2f}; KL-asympt UB {ub_kl:.2f}; MH-level (random W_Z) {mh:.2f}")

print("== (5) SL_2(p) Weil lattice = Craig A_{p-1}^{(m)}, m=(p+1)/4: min norm and f*tau^2 ==")
def craig_min(p, mm, bound):
    # A_{p-1}^{(m)} = (1-x)^m Z[x]/(Phi_p) ... realize inside Z^p: vectors c in Z^p with sum_i i^j c_i = 0 mod p for j=0..m-1 and sum c_i = 0 (j=0 exact)
    # search vectors with small support and entries +-1 (norm = support size)
    best = None
    for wt in range(2, bound + 1, 2):
        half = wt // 2
        # entries: half +1's and half -1's (sum zero), first +1 at position 0 (cyclic shift invariance)
        for pos in itertools.combinations(range(1, p), half - 1):
            P = (0,) + pos
            rest = [i for i in range(p) if i not in P]
            for neg in itertools.combinations(rest, half):
                good = True
                for j in range(1, mm):
                    if (sum(pow(i, j, p) for i in P) - sum(pow(i, j, p) for i in neg)) % p:
                        good = False; break
                if good:
                    return wt
    return None
for p, bound in [(7, 4), (11, 6), (19, 10)]:
    mm = (p + 1) // 4
    mn = craig_min(p, mm, bound)
    r = (p - 1) // 2
    # Z-lattice: dim p-1, det = p^{2m-1}; K=Q(sqrt(-p)), tau_line^2 = gamma*sqrt(p)/2 (valid if lattice is p-modular i.e. 2m-1=(p-1)/2)
    gam = mn / p ** ((2 * mm - 1) / (p - 1)) if mn else float('nan')
    tau2 = gam * math.sqrt(p) / 2
    print(f"  p={p}: Craig m={mm}, det=p^{2*mm-1}, found min norm (+-1 vectors, wt<= {bound}) = {mn} (theory >= {2*mm}); gamma={gam:.4f}; tau_line^2<={tau2:.4f} need r={r}; f*tau^2={tau2/r:.4f}")

print("== (6) unramified base change: K=Q(sqrt(-15)), E=Q(sqrt(-3),sqrt(5)) (Hilbert class field), W_Z=O_E, r=2 ==")
# integral basis of O_E: 1, w3=(1+sqrt(-3))/2, w5=(1+sqrt5)/2, w3*w5 ; embeddings (sqrt-3 -> +-i sqrt3, sqrt5 -> +-sqrt5)
def emb(e3, e5):
    w3 = (1 + e3 * 1j * math.sqrt(3)) / 2; w5 = (1 + e5 * math.sqrt(5)) / 2
    return [1, w3, w5, w3 * w5]
E = np.array([emb(a, b) for a in (1, -1) for b in (1, -1)])  # 4 embeddings x 4 basis
G = (E.conj().T @ E).real  # Gram of Tr_{E/Q}(x conj y)
cW = math.sqrt(abs(np.linalg.det(G)))
# K = Q(sqrt(-15)), O_K basis 1, (1+sqrt(-15))/2 ; sqrt(-15)=sqrt(-3)*sqrt5 = (2w3-1)(2w5-1) = 4w3w5-2w3-2w5+1 -> (1+sqrt-15)/2 = 1 - w3 - w5 + 2 w3 w5
BK = np.array([[1, 0, 0, 0], [1, -1, -1, 2]], dtype=float)
cU = math.sqrt(abs(np.linalg.det(BK @ G @ BK.T)))
print(f"  covol(O_E)={cW:.6f} (sqrt|d_E|=15), covol(O_K in E-trace form)={cU:.6f} (=2*sqrt(15)={2*math.sqrt(15):.6f})")
report("O_E over O_K, U=K (subfield)", cU, 2, cW, 4)
# brute force all K-lines through small w: c_U for U = K w  (lattice O_E cap K w found by saturation: use w*O_K and divide by index via gcd trick -> lower bound on tau via |w|^2/2)
mn = min(float(np.array(v) @ G @ np.array(v)) for v in itertools.product(range(-3, 4), repeat=4) if any(v))
print(f"  lambda_1(O_E)^2 = {mn:.4f}; for d=2: tau_line^2 = min h(w) = lambda_1^2/2 = {mn/2:.4f} (need r=2): borderline exactly")

print("== (7) g = |G|/dim_Q A = sum rho^2 : cyclic check and fpf examples ==")
from sympy import totient, mobius, divisors
for m in [6, 30, 210, 2310, 30030]:
    ssum = sum(int(totient(dd)) * (int(mobius(dd)) / int(totient(dd))) ** 2 for dd in divisors(m))
    print(f"  C_{m} on Q(zeta_{m}): sum rho^2={ssum:.6f}  m/phi(m)={m/int(totient(m)):.6f}  e^gamma loglog m={math.exp(0.5772156649)*math.log(math.log(m)):.4f}")
# D_8 on Q^2: character 2,-2,0,0,0 ; classes sizes 1,1,2,2,2
print("  D_8 on Q^2: sum rho^2 =", (1 + 1 + 0) , " |G|/dim A = 8/4 =", 8 / 4)
# D_12 on Q^2 (A_2): chi: e:2, r,r^5:1, r^2,r^4:-1, r^3:-2, reflections 0
print("  D_12 on Q^2: sum rho^2 =", (4 + 2 * 1 + 2 * 1 + 4) / 4, " |G|/dim A = 12/4 =", 12 / 4)
# SL_2(5) on icosians W=Q^8: chi_W = 2*(chi_2+chi_2') ; chi_2 values by element order: 1:2(1), 2:-2(1), 3:-1(20), 6:1(20), 4:0(30), 5:(-1+-sqrt5)/2 (12+12), 10:(1+-sqrt5)/2 (12+12)
phi = (1 + math.sqrt(5)) / 2
vals = [(1, 2, 2), (1, -2, -2), (20, -1, -1), (20, 1, 1), (30, 0, 0), (12, phi - 1, -phi), (12, -phi, phi - 1), (12, phi, 1 - phi), (12, 1 - phi, phi)]
tot = sum(c * ((2 * (x + y)) / 8) ** 2 for c, x, y in vals)
print(f"  SL_2(5) on icosians (Q^8): sum rho^2 = {tot:.6f}  |G|/dim_Q A = 120/8 = {120/8}")
