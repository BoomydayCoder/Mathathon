import math
from math import log, exp, sqrt, pi, e
def zeta(s, N=20000):
    return sum(n**-s for n in range(1, N)) + N**(1-s)/(s-1)
def kron(a, p):  # Legendre symbol, p = 1 mod 4 prime -> chi_p(n) = (n/p)
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1
def Lchi(s, p, N=20000):
    return sum(kron(n, p)*n**-s for n in range(1, N))
def zetaF(s, p): return zeta(s)*Lchi(s, p)
print("== A1: Solomon zeta of O^r: prod_{i=0}^{2r-1} zeta_F(2k-i), converges iff k>=r+1 ==")
for p in (13, 37, 61, 157):
    r = (p-1)//4
    for k in (r+1, r+2, 2*r):
        prod = 1.0
        for i in range(2*r):
            s = 2*k - i
            prod *= zetaF(s, p) if s < 40 else 1.0
        print(f"p={p} r={r} k={k}: zeta_(O^r)(k) = {prod:.4f}")
print("uniform bound prod_{j>=3} zeta(j)^2 =", (math.prod(zeta(j) for j in range(3, 60)))**2)
print("== A3: breathing mode. typical q/q0 = 1+(r-1/2)/k ; trace part = r/k = n/(2D) ==")
def Vfinal_log_p(c, V0exp):   # exponent of p in V_final for k=c*r, V0 ~ p^V0exp, wD ~ p^5, n^2 ~ p^4
    y = 1.0/c
    return V0exp*y/(1+y) + 5.0/(1+y)
for c in (1, 1.5, 2, 3, 5, 10, 100):
    a = Vfinal_log_p(c, 3); b = Vfinal_log_p(c, 2)
    print(f"k={c}r: V~p^{a:.3f} (V0~|G|) -> n^{a/2:.3f} ; V~p^{b:.3f} (V0~n) -> n^{b/2:.3f}   [Klartag n^2 = p^4]")
print("== A2: borderline classes, f=1/3 and lines, at V=wD/2 (contact reach) ==")
for p in (13, 37, 61, 157, 1009):
    r = (p-1)//4
    for c in (1, 3):
        k = max(r+1, int(c*r)); n = 8*r*k; w = p*(p*p-1); D = 4*k*k-2*k
        Vr = w*D/2
        def T(f, V): return f**-0.5 * pi**((f-1)/2) * V**f * n**(-(1-f)/2)
        T3 = T(1/3, Vr); Tl = T(1/r, Vr)
        con_C3 = T3*3/(2*(p-1)); con_S = T3/24
        V0max_S = n*(6/1.18)**3
        print(f"p={p} r={r} k={k} n={n} D={D} |G|={w}: T_(f=1/3)={T3:.1f} -> constraints C3-type {con_C3:.2f}, SL2(3)-type {con_S:.2f} (vs D);"
              f" T_line={Tl:.3f}; lines needed to fill D: {D*w/Tl:.2e}; V0max from SL2(3)-type={V0max_S:.0f} vs |G|/2={w/2:.0f};"
              f" sum_U T_U(V)/V at V=|G|/4 (SL2(3)-type, p^3/24 subspaces): {(w/24)*T(1/3, w/4)/(w/4):.2f}")
print("== A6: isotypic components. R_i = (|disc_trd Gamma_i| / (e_i chi_i(1))^{a_i})^m * 2^{(free 2-exp) - e} ==")
print("SL2(3) psi_0 component: R = 2^-e, 0<=e<=4r/3 (borderline iff e=0)")
for p in (61, 421, 541, 661):
    r = (p-1)//4; m = 2*r/15
    print(f"2I icosian component p={p}: r={r} m={m:.2f}: R(e=0) = (625/256)^m = {(625/256)**m:.3g}; need 2-adic gluing e >= {1.2877*m:.2f} (max {8*m:.1f})")
