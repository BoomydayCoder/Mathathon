#!/usr/bin/env python3
"""
zeta5opt_equilibrium.py -- exact (Cauchy-transform) equilibrium measure for the Fauzan-type field, a >= 4.

Variable t = y^2/K^2 (as in Fauzan section 6).  Field and problem:
    V(t) = 2 pi sqrt(t) + J(1,t) - a J(alpha,t),   J(c,t) = int_0^c log(t+u^2) du = c log(t+c^2) - 2c + 2 sqrt(t) arctan(c/sqrt(t)),
    minimise  int V dmu - I(mu)  over positive measures mu on [0,inf) of mass lam;   I(mu) = iint log|t-s| dmu dmu.
Charge representation:  V'(t)/2 = G(t) = int_0^inf om(u) du/(t+u^2),
    om(u) = (3-a)/2 on (0,alpha),  3/2 on (alpha,1),  1 on (1,inf).
For a > 3 the support is ONE interval [A,B] with 0 < A < B (soft edges), R(z) = sqrt((z-A)(z-B)), r_u = sqrt((A+u^2)(B+u^2)):
    Cauchy transform   F(z) = int dmu(s)/(z-s) = G(z) + R(z) Hh(z),   Hh(z) = int_0^inf om(u) du / (r_u (z+u^2)),
    density            mu'(t) = -(1/pi) sqrt((t-A)(B-t)) Hh(t)      (A<t<B),
    endpoint equations (C1)  int_0^inf om(u)/r_u du = 0,
                       (C2)  -((B-A)^2/4) int_0^inf om(u) du / (r_u (r_u + u^2 + (A+B)/2)) = lam.
Given (C1): Hh(t) = int (om(u)/r_u) [1/(t+u^2) - 1/(t+alpha^2)] du < 0 for all t >= 0 (integrand <= 0 pointwise), hence
    mu' > 0 on (A,B);  d/dt (2U^mu - V) = 2 R Hh  is < 0 for t > B and > 0 for 0 <= t < A   (Euler-Lagrange inequality, analytic).
Constants:  ell = V - 2U^mu on [A,B] = (Chebyshev mean of V over [A,B]) - 2 lam log((B-A)/4),
            E* = (int V dmu + lam ell)/2,   C_* = -2 lam + 2 a alpha lam (1 - log alpha) + 3 lam^2 - 2 lam^2 log(2 lam),
            U = -E* + C_*   (log F_K(zeta(5)) <= U K^2 + O(K log K)).
Two implementations: float (numpy, Gauss-Legendre / Gauss-Chebyshev) for scans, and mpmath (any precision) for certificates
(non-rigorous quadrature; "numerical to N digits").
Usage: python3 zeta5opt_equilibrium.py            (Fauzan's parameters + the optimised ones, 30 digits)
"""
import sys, math
import numpy as np
import mpmath as mp


# ------------------------------------------------------------------ float version
_gl = {}


def _glnodes(n):
    if n not in _gl:
        _gl[n] = np.polynomial.legendre.leggauss(n)
    return _gl[n]


def _charges(A, B, a, alpha, n=300):
    """nodes u_i and weights W_i with  int om(u) f(u)/r_u du ~ sum W_i f(u_i)."""
    x, w = _glnodes(n)
    sA = math.sqrt(A)
    s1, s2 = math.asinh(alpha / sA), math.asinh(1 / sA)
    us, Ws = [], []
    for lo, hi, om in ((0.0, s1, (3 - a) / 2), (s1, s2, 1.5)):
        s = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
        u = sA * np.sinh(s)
        us.append(u); Ws.append(om * 0.5 * (hi - lo) * w / np.sqrt(B + u * u))
    v = 0.5 * x + 0.5
    us.append(1 / v); Ws.append(1.0 * 0.5 * w / np.sqrt((1 + A * v * v) * (1 + B * v * v)))
    return np.concatenate(us), np.concatenate(Ws)


def _C1(A, B, a, alpha):
    return _charges(A, B, a, alpha)[1].sum()


def _C2(A, B, a, alpha):
    u, W = _charges(A, B, a, alpha)
    r = np.sqrt((A + u * u) * (B + u * u))
    return -((B - A) ** 2 / 4) * (W / (r + u * u + (A + B) / 2)).sum()


def Vt(t, a, alpha):
    st = np.sqrt(t)
    J = lambda c: c * np.log(t + c * c) - 2 * c + 2 * st * np.arctan(c / st)
    return 2 * np.pi * st + J(1.0) - a * J(alpha)


def Cstar(a, alpha, lam):
    return -2 * lam + 2 * a * alpha * lam * (1 - math.log(alpha)) + 3 * lam ** 2 - 2 * lam ** 2 * math.log(2 * lam)


def solve_float(a, alpha, lam, n=400):
    from scipy.optimize import brentq
    def A_of_B(B):
        f = lambda lA: _C1(math.exp(lA), B, a, alpha)
        lo, hi = math.log(B) - 40, math.log(B) - 1e-9
        return math.exp(brentq(f, lo, hi, xtol=1e-15, rtol=1e-15))
    def mass(B):
        return _C2(A_of_B(B), B, a, alpha) - lam
    # V'(B) > 0 is needed for (C1) to be solvable: find the minimum t0 of V first
    g = lambda t: math.pi + math.atan(1 / math.sqrt(t)) - a * math.atan(alpha / math.sqrt(t))
    t0 = brentq(g, 1e-12, 10.0)
    Blo, Bhi = t0 * (1 + 1e-6), 50.0
    B = brentq(mass, Blo, Bhi, xtol=1e-15, rtol=1e-15)
    A = A_of_B(B)
    u, W = _charges(A, B, a, alpha)
    a0, b0 = math.sqrt(A), math.sqrt(B)
    mx, hx = (a0 + b0) / 2, (b0 - a0) / 2
    th = np.arange(1, n + 1) * np.pi / (n + 1)
    x = mx + hx * np.cos(th)
    t = x * x
    psi = -(W[None, :] / (t[:, None] + u[None, :] ** 2)).sum(1)
    wq = hx ** 2 / (n + 1) * np.sin(th) ** 2 * 2 * x * np.sqrt((x + a0) * (x + b0)) * psi   # weights of mu at t = x^2
    massq = wq.sum()
    IV = (wq * Vt(t, a, alpha)).sum()
    nt = 40000
    tc = (A + B) / 2 + (B - A) / 2 * np.cos((2 * np.arange(1, nt + 1) - 1) * np.pi / (2 * nt))
    ell = Vt(tc, a, alpha).mean() - 2 * lam * math.log((B - A) / 4)
    E = (IV + lam * ell) / 2
    return dict(A=A, B=B, ell=ell, IV=IV, E=E, U=-E + Cstar(a, alpha, lam), mass=massq, psimin=psi.min())


# ------------------------------------------------------------------ mpmath version
def solve_mp(a, alpha, lam, dps=30, n=None, verbose=True):
    mp.mp.dps = dps + 10
    a = mp.mpf(a); alpha = mp.mpf(alpha); lam = mp.mpf(lam)
    om1, om2, om3 = (3 - a) / 2, mp.mpf(3) / 2, mp.mpf(1)

    def charge_int(A, B, f_s, f_v):
        """int om(u) f/r_u du with f given in the two parametrisations (u = sqrt(A) sinh s; u = 1/v)."""
        sA = mp.sqrt(A)
        s1, s2 = mp.asinh(alpha / sA), mp.asinh(1 / sA)
        g = lambda s: f_s(sA * mp.sinh(s)) / mp.sqrt(B + A * mp.sinh(s) ** 2)
        gv = lambda v: f_v(v) / mp.sqrt((1 + A * v * v) * (1 + B * v * v))
        return om1 * mp.quad(g, [0, s1]) + om2 * mp.quad(g, [s1, s2]) + om3 * mp.quad(gv, [0, 1])

    C1 = lambda A, B: charge_int(A, B, lambda u: 1, lambda v: 1)

    def C2(A, B):
        m_ = (A + B) / 2
        fs = lambda u: 1 / (mp.sqrt((A + u * u) * (B + u * u)) + u * u + m_)
        fv = lambda v: v * v / (mp.sqrt((1 + A * v * v) * (1 + B * v * v)) + 1 + m_ * v * v)
        return -((B - A) ** 2 / 4) * charge_int(A, B, fs, fv)

    f0 = solve_float(float(a), float(alpha), float(lam))
    sol = mp.findroot(lambda lA, lB: (C1(mp.e ** lA, mp.e ** lB), C2(mp.e ** lA, mp.e ** lB) - lam),
                      (mp.log(f0['A']), mp.log(f0['B'])), tol=mp.mpf(10) ** (-dps - 5))
    A, B = mp.e ** sol[0], mp.e ** sol[1]
    resid = (C1(A, B), C2(A, B) - lam)
    Hh = lambda t: charge_int(A, B, lambda u: 1 / (t + u * u), lambda v: v * v / (1 + t * v * v))

    def V(t):
        st = mp.sqrt(t)
        J = lambda c: c * mp.log(t + c * c) - 2 * c + (2 * st * mp.atan(c / st) if st != 0 else 0)
        return 2 * mp.pi * st + J(mp.mpf(1)) - a * J(alpha)

    m_, hw = (A + B) / 2, (B - A) / 2
    a0, b0 = mp.sqrt(A), mp.sqrt(B)
    mx, hx = (a0 + b0) / 2, (b0 - a0) / 2
    # mu-quadrature in x = sqrt(t) (Chebyshev, second kind); nearest singularity of the integrand: x = 0
    rho = (mx + mp.sqrt(a0 * b0)) / hx
    if n is None:
        n = int(1.2 * (dps + 8) * mp.log(10) / (2 * mp.log(rho))) + 30
    th = [k * mp.pi / (n + 1) for k in range(1, n + 1)]
    xs = [mx + hx * mp.cos(x) for x in th]
    t = [x * x for x in xs]
    psi = [-Hh(x) for x in t]
    wq = [hx ** 2 / (n + 1) * mp.sin(q) ** 2 * 2 * x * mp.sqrt((x + a0) * (x + b0)) * p for q, x, p in zip(th, xs, psi)]
    Vn = [V(x) for x in t]
    mass = mp.fsum(wq)
    IV = mp.fsum(w * v for w, v in zip(wq, Vn))
    # Chebyshev mean of V over [A,B] in the t-variable (V is closed form; singularity t = 0 is close, so many nodes)
    rho_t = (m_ + mp.sqrt(A * B)) / hw
    nt = int(1.2 * (dps + 8) * mp.log(10) / (2 * mp.log(rho_t))) + 50
    thc = [(2 * j - 1) * mp.pi / (2 * nt) for j in range(1, nt + 1)]
    Vc = [V(m_ + hw * mp.cos(x)) for x in thc]
    ell = mp.fsum(Vc) / nt - 2 * lam * mp.log(hw / 2)
    E = (IV + lam * ell) / 2
    Cs = -2 * lam + 2 * a * alpha * lam * (1 - mp.log(alpha)) + 3 * lam ** 2 - 2 * lam ** 2 * mp.log(2 * lam)
    # --- Euler-Lagrange on the support via Chebyshev coefficients (t-variable):  V - 2U = ell  <=>  v_k + (4/k) c_k = 0
    kmax = 12
    el_err = mp.mpf(0)
    for k in range(1, kmax + 1):
        vk = 2 * mp.fsum(v * mp.cos(k * x) for v, x in zip(Vc, thc)) / nt
        ck = mp.fsum(w * mp.cos(k * mp.acos((s_ - m_) / hw)) for w, s_ in zip(wq, t))
        el_err = max(el_err, abs(vk + 4 * ck / k))
    # --- off the support: phi(t) = 2U(t) - V(t) + ell < 0 (direct quadrature of the potential)
    U = lambda x: mp.fsum(w * mp.log(abs(x - s_)) for w, s_ in zip(wq, t))
    off = []
    for x in [mp.mpf(0), A / 4, A / 2, B + (B - A) / 4, 2 * B, 5 * B, 20 * B, 200 * B]:
        off.append((x, 2 * U(x) - V(x) + ell))
    # --- independent check of ell: phi(2B) must equal int_B^{2B} 2 R Hh  (uses neither ell nor U)
    phi2B = mp.quad(lambda s_: 2 * mp.sqrt((s_ - A) * (s_ - B)) * Hh(s_), [B, 2 * B])
    E2 = phi2B - (2 * U(2 * B) - V(2 * B) + ell)
    # --- density maximum (for the smoothing error) and positivity
    dens = [mp.sqrt((x - A) * (B - x)) * p / mp.pi for x, p in zip(t, psi)]
    out = dict(A=A, B=B, resid=resid, mass=mass, ell=ell, IV=IV, E=E, E2=E2, Cs=Cs, U=-E + Cs, el_err=el_err, off=off,
               psimin=min(psi), densmax=max(dens), n=n, rho=rho, HhA=Hh(A), HhB=Hh(B), Hh0=Hh(mp.mpf(0)))
    if verbose:
        d = dps
        print(f"  parameters a={a}, alpha={mp.nstr(alpha, 12)}, lam={mp.nstr(lam, 12)}   (working precision {dps + 10} digits, {n} nodes in x=sqrt(t) (rho={mp.nstr(rho, 6)}), {nt} nodes in t)")
        print(f"  support [A,B] (t-variable): A = {mp.nstr(A, d)}\n                              B = {mp.nstr(B, d)}   (sqrt: [{mp.nstr(mp.sqrt(A), 12)}, {mp.nstr(mp.sqrt(B), 12)}])")
        print(f"  residuals of (C1),(C2): {mp.nstr(resid[0], 3)}, {mp.nstr(resid[1], 3)};   quadrature mass - lam = {mp.nstr(mass - lam, 3)}")
        print(f"  -Hh > 0 at all nodes: min psi = {mp.nstr(min(psi), 8)};  Hh(0) = {mp.nstr(out['Hh0'], 8)}, Hh(A) = {mp.nstr(out['HhA'], 8)}, Hh(B) = {mp.nstr(out['HhB'], 8)}; max density = {mp.nstr(max(dens), 8)}")
        print(f"  E-L equality on [A,B]: max_k<= {kmax} |v_k + 4 c_k/k| = {mp.nstr(el_err, 3)}")
        print("  E-L strict inequality off the support, phi = 2U - V + ell: " + ", ".join(f"phi({mp.nstr(x, 4)})={mp.nstr(v, 6)}" for x, v in off))
        print(f"  ell = {mp.nstr(ell, d)}\n  int V dmu = {mp.nstr(IV, d)}")
        print(f"  E* = (int V dmu + lam ell)/2 = {mp.nstr(E, d)}")
        print(f"  independent check of ell: int_B^2B 2 R Hh dt - phi(2B) = {mp.nstr(E2, 3)}")
        print(f"  C_* = {mp.nstr(Cs, d)}\n  U = -E* + C_* = {mp.nstr(-E + Cs, d)}")
    return out


if __name__ == "__main__":
    cases = [(6, mp.mpf(3) / 40, mp.mpf(37) / 40, "Fauzan (a=6, alpha=3/40, lam=37/40); his proved U = -2733991/2000000 = -1.3669955")]
    for arg in sys.argv[1:]:
        aa, num, den = arg.split(',')
        cases.append((int(aa), mp.mpf(int(num)) / int(den), 1 - mp.mpf(int(num)) / int(den), f"a={aa}, alpha={num}/{den}, lam=1-alpha"))
    for a, al, la, name in cases:
        print("=== " + name)
        f = solve_float(a, float(al), float(la))
        print(f"  float solver: A={f['A']:.12f} B={f['B']:.12f} E*={f['E']:.12f} U={f['U']:.12f} mass={f['mass']:.12f}")
        solve_mp(a, al, la, dps=30)
