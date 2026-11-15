"""Numerical checks (report, Parts II-III): J_n = -sum_{t>=1} R_n'(t) for n = 0..5 with 40-digit
arithmetic, the ratios J_{n-1}/J_n, the bound J_n <= 2 zeta(3) 33^{-n}, and J_n = 2 a_n zeta(3) - 2 b_n
for n = 1, 2 (J_1 = 10 zeta(3) - 12, J_2 = 146 zeta(3) - 351/2)."""
import sympy as sp, mpmath as mp
mp.mp.dps = 40
t = sp.symbols('t')
def r(N): return sp.prod([j - t for j in range(1, N+1)]) / sp.prod([t + j for j in range(0, N+1)])
J = []
for N in range(6):
    f = sp.lambdify(t, sp.diff(r(N)**2, t), 'mpmath')
    J.append(-mp.nsum(f, [1, mp.inf]))
for N, v in enumerate(J):
    print(f"J_{N} = {mp.nstr(v, 15)}   bound 2 zeta(3) 33^-n holds: {v <= 2*mp.zeta(3)*mp.mpf(33)**(-N)}")
print("ratios J_{n-1}/J_n:", [mp.nstr(J[k-1]/J[k], 6) for k in range(1, 6)], " 17+12 sqrt2 =", mp.nstr(17+12*mp.sqrt(2), 6))
print("J_1 - (10 zeta(3) - 12) =", mp.nstr(J[1] - (10*mp.zeta(3) - 12), 5))
print("J_2 - (146 zeta(3) - 351/2) =", mp.nstr(J[2] - (146*mp.zeta(3) - mp.mpf(351)/2), 5))
