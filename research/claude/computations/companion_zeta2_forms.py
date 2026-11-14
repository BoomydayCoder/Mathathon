"""By-product check: the companion sums G_n = sum_{t>=1} R_n(t) (linear forms in 1, zeta(2))
satisfy Apery's zeta(3) recurrence with right-hand side -8(2n+1):
   (n+1)^3 G_{n+1} - (2n+1)(17n^2+17n+5) G_n + n^3 G_{n-1} = -8(2n+1),  n>=1.
Derivation: sum Zudilin's identity (6) over t>=1; S_n(t) -> -8(2n+1) as t->oo, S_n(1)=0.
Exact check via partial fractions: G_n = (u_n/2) zeta(2) - w_n,
   w_n = sum_k A2k H_k^(2) + sum_k A1k H_k   (the zeta(2)-part must satisfy L=0, the rational part L[w] = +8(2n+1))."""
from fractions import Fraction as Fr
from math import comb
def H(k, j=1): return sum(Fr(1, l**j) for l in range(1, k+1))
def parts(n):
    A2 = [comb(n+k, n)**2 * comb(n, k)**2 for k in range(n+1)]
    A1 = [2*A2[k]*(2*H(k) - H(n+k) - H(n-k)) for k in range(n+1)]
    return A2, A1
def uw(n):
    A2, A1 = parts(n)
    return sum(A2), sum(A2[k]*H(k, 2) + A1[k]*H(k) for k in range(n+1))
ok = True
vals = [uw(n) for n in range(0, 32)]
for n in range(1, 31):
    b = (2*n+1)*(17*n*n+17*n+5)
    Lu = (n+1)**3*vals[n+1][0] - b*vals[n][0] + n**3*vals[n-1][0]
    Lw = (n+1)**3*vals[n+1][1] - b*vals[n][1] + n**3*vals[n-1][1]
    if Lu != 0 or Lw != 8*(2*n+1): ok = False
print("G_0 = zeta(2) - %s, G_1 = %s zeta(2) - %s, G_2 = %s zeta(2) - %s" % (vals[0][1], vals[1][0], vals[1][1], vals[2][0], vals[2][1]))
print("CHECK L[u/2]=0 and L[w]=8(2n+1) exactly, n=1..30 (so L[G] = -8(2n+1)):", "PASS" if ok else "FAIL")
import mpmath as mp
mp.mp.dps = 40
for n in (1, 5, 10, 20):
    print("n=%d  G_n = %s" % (n, mp.nstr(vals[n][0]*mp.zeta(2) - mp.mpf(vals[n][1].numerator)/vals[n][1].denominator, 12)))
