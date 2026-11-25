import math, itertools, cmath
from sympy import totient, divisors
# Parseval: sum_k tr(tau^k)^2 = M*C, tr(tau^k | Q(zeta_n)) = Ramanujan sum c_n(k) (computed from roots of unity directly)
def ram(n, k): return sum(cmath.exp(2j*math.pi*a*k/n) for a in range(1, n+1) if math.gcd(a, n) == 1).real
for N in ([1,3,15], [1,5,7,35,105], [3,25,45], [4,12,20,60], [1,3,6,9]):
    M = math.lcm(*N); C = sum(int(totient(n)) for n in N)
    s = sum(sum(ram(n, k) for n in N)**2 for k in range(M))
    print(f"N={N}: M/C={M/C:.10f}  sum_g rho(g)^2={s/C**2:.10f}")
# count of N within divisors(1155) with lcm=1155 and C<480
D = [d for d in divisors(1155)]; phi = {d: int(totient(d)) for d in D}
cnt = tot = 0
for mask in range(1, 1 << len(D)):
    N = [D[i] for i in range(len(D)) if mask >> i & 1]
    if math.lcm(*N) != 1155: continue
    tot += 1
    if sum(phi[n] for n in N) < 480: cnt += 1
print("M=1155: index sets with lcm=M:", tot, "; with C<phi(M)=480:", cnt)
