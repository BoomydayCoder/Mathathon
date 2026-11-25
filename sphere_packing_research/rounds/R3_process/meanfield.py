"""Mean-field model of Klartag's evolving-ellipsoid process (heuristic, Poisson lattice model).
Units: lattice covolume = Vol(B^n); V = relative volume = det(A)^(-1/2); D = n(n+1)/2.
Checks:
 (1) closed form K(u) = V0*(2 e^u Phi(sqrt(2u)) - 1) for the expected number of contact POINTS
     (sharp reflection principle, clock variable u = n^2*tau/8) against numerical quadrature;
 (2) exact identity E|P_E(w (x) w)|^2 = (2 dim E + |P_E Id|^2)/(n(n+2)) (w uniform on sphere), Monte Carlo;
 (3) solve model K (full variance, Klartag/Rezvyakova accounting) and model R (variance reduced by N_t/D).
"""
import numpy as np
from scipy import integrate, optimize, stats

Phi = stats.norm.cdf

def K_closed(u, V0=1.0):
    return V0 * (2 * np.exp(u) * Phi(np.sqrt(2 * u)) - 1)

def K_quad(n, tau, a0=1.0):
    # exact (non-linearised) density of lattice points in depth s = a0 - 1/|x|^2: (n/2)(a0-s)^(-n/2-1)
    f = lambda s: (n / 2) * (a0 - s) ** (-n / 2 - 1) * 2 * stats.norm.sf(s / np.sqrt(tau))
    val, _ = integrate.quad(f, 0, min(a0 * 0.5, 12 * np.sqrt(tau)), limit=400)
    return val

print("(1) contact-count formula: n, u, closed form, quadrature with exact density (a0=1)")
for n in (50, 200, 1000):
    for u in (0.5, 2.0, 2 * np.log(n)):
        tau = 8 * u / n ** 2
        print(f"   n={n:5d} u={u:7.3f} closed={K_closed(u):12.4f} quad={K_quad(n, tau):12.4f}")

print("(2) projection identity, Monte Carlo")
rng = np.random.default_rng(1)
n = 8
iu = np.triu_indices(n)
def vec(S):  # isometry Sym_n -> R^D for the trace inner product
    M = S.copy() * np.sqrt(2); M[np.diag_indices(n)] /= np.sqrt(2)
    return M[iu]
for k in (1, 5, 20, 30):
    X = rng.standard_normal((k, n)); X /= np.linalg.norm(X, axis=1)[:, None]
    B = np.array([vec(np.outer(x, x)) for x in X]).T
    Q, _ = np.linalg.qr(B)
    W = rng.standard_normal((40000, n)); W /= np.linalg.norm(W, axis=1)[:, None]
    vals = [np.sum((Q.T @ vec(np.outer(w, w))) ** 2) for w in W[:20000]]
    pid = np.sum((Q.T @ vec(np.eye(n))) ** 2)
    print(f"   n={n} k={k:2d}: MC={np.mean(vals):.5f}  formula={(2 * k + pid) / (n * (n + 2)):.5f}  |P_E Id|^2={pid:.3f} (n={n})")

def cpairs(u, V0):
    return V0 * (np.exp(u) * Phi(np.sqrt(2 * u)) - 0.5)

print("(3) mean-field predictions (V0 = (1-1/n)^n as in Klartag; T_K = 16 log n / n^2)")
print("    n      D   | model R: V(T_K)/n^2  Vfinal/n^2  Vfinal/D | model K: V(T_K)/n^2  max_t V/n^2 (Rezvyakova 1/(2e)=0.1839 needs V0->0)")
for n in (10, 30, 100, 1000, 10 ** 4, 10 ** 6):
    D = n * (n + 1) / 2
    V0 = (1 - 1 / n) ** n
    TK = 16 * np.log(n) / n ** 2
    # model R: u = n^2 tau/8 is the clock of the contact count; log V = log V0 + (D/4) tau = log V0 + (2D/n^2) u
    # d tau/dt = 1 - c/D  =>  du/dt = (n^2/8)(1 - c(u)/D)
    rhsR = lambda t, y: [(n ** 2 / 8) * max(0.0, 1 - cpairs(y[0], V0) / D)]
    solR = integrate.solve_ivp(rhsR, [0, TK], [0.0], rtol=1e-10, atol=1e-12)
    uTK = solR.y[0, -1]
    VR_TK = V0 * np.exp(2 * D / n ** 2 * uTK)
    uf = optimize.brentq(lambda u: cpairs(u, V0) - D, 0, 60)
    VR_fin = V0 * np.exp(2 * D / n ** 2 * uf)
    # model K: c(t) = cpairs(n^2 t/8); log V = log V0 + (1/4) int (D - c) dt
    logVK = lambda T: np.log(V0) + 0.25 * integrate.quad(lambda t: D - cpairs(n ** 2 * t / 8, V0), 0, T)[0]
    VK_TK = np.exp(logVK(TK))
    Tstar = 8 / n ** 2 * uf   # c = D
    VK_max = np.exp(logVK(Tstar))
    print(f"{n:8d} {D:12.0f} | {VR_TK / n ** 2:10.4f} {VR_fin / n ** 2:10.4f} {VR_fin / D:9.4f} | {VK_TK / n ** 2:10.4f} {VK_max / n ** 2:10.4f}")

print("(4) logistic closed form for model R (large n): V(t) = D/(1+(D/V0-1) exp(-D t/4)); at T_K -> n^2/3 when V0=1")
for n in (100, 10 ** 4):
    D = n * (n + 1) / 2
    for V0 in (1.0, np.exp(-1.0)):
        TK = 16 * np.log(n) / n ** 2
        print(f"   n={n} V0={V0:.3f}: V(T_K)/n^2 = {D / (1 + (D / V0 - 1) * np.exp(-D * TK / 4)) / n ** 2:.4f}")

print("(5) Poisson-model volume law of a vertex: Gamma(D, scale 2): mean 2D; extreme over exp(S) vertices (Cramer form lam-1-log lam = S/D)")
for n in (100, 1000, 10 ** 6):
    D = n * (n + 1) / 2
    for kappa in (0.5, 1.0, 2.0):
        lam = optimize.brentq(lambda l: l - 1 - np.log(l) - kappa * np.log(n), 1.0 + 1e-9, 1e4)
        print(f"   n={n} #vertices=exp({kappa} D log n): lam={lam:.3f}, V=2D*lam={2 * D * lam / (n ** 2 * np.log(n)):.3f} n^2 log n")
# sticky-exploration barrier in the Poisson model: P(V_T >= V0 + 2 lam D) <= P(Poisson(lam D) <= D)
for D in (55, 5050):
    for lam in (1.5, 2.0, 3.0):
        print(f"   D={D} lam={lam}: P(Poi(lam D) <= D) = {stats.poisson.cdf(D, lam * D):.3e}  vs exp(-D(lam-1-log lam)) = {np.exp(-D * (lam - 1 - np.log(lam))):.3e}")
# neighbourhood bound of Lemma N: log P_{V_n/2} >= log( n (4n+2)^(-n^2) kappa_{n^2} / prod zeta )
from scipy.special import gammaln, zeta
for n in (2, 8, 24, 100):
    lb = np.log(n) - n ** 2 * np.log(4 * n + 2) + (n ** 2 / 2) * np.log(np.pi) - gammaln(n ** 2 / 2 + 1) - sum(np.log(zeta(k)) for k in range(2, n + 1))
    print(f"   Lemma N: n={n}: log P_(V_n/2) >= {lb:.1f} ; claimed simplification -2n^2 log n - n^2 = {-2 * n ** 2 * np.log(n) - n ** 2:.1f}")
