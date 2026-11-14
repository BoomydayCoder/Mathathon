"""Numerical evidence (NOT a proof) for the bound D_n < 3^n used in Zudilin's proof,
and for the weaker requirement log D_n < (4/3) log(1+sqrt 2) n = 1.1752 n actually needed.
D_n = lcm(1..n); log D_n = psi(n) (Chebyshev function)."""
from math import log, sqrt
N = 10**5
spf = list(range(N + 1))
for i in range(2, int(N**0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, N + 1, i):
            if spf[j] == j:
                spf[j] = i
psi = 0.0; worst = (0, 0.0); ok3 = True
c_need = 4 / 3 * log(1 + sqrt(2))
for n in range(2, N + 1):
    p = spf[n]; m = n
    while m % p == 0:
        m //= p
    if m == 1:            # n is a prime power p^k: lcm gains a factor p
        psi += log(p)
    r = psi / n
    if r > worst[1]:
        worst = (n, r)
    if not psi < n * log(3):
        ok3 = False
print("max over 2<=n<=%d of log(D_n)/n = %.6f at n=%d" % (N, worst[1], worst[0]))
print("log 3 = %.6f ; needed constant (4/3)log(1+sqrt2) = %.6f" % (log(3), c_need))
print("CHECK D_n < 3^n for 1<=n<=%d:" % N, "PASS" if ok3 else "FAIL")
print("CHECK log D_n < 1.1752 n     :", "PASS" if worst[1] < c_need else "FAIL")
print("psi(N)/N at N=%d: %.6f (PNT predicts -> 1)" % (N, psi / N))
