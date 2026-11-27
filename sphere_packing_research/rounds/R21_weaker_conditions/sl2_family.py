# R21: admissible rational characters of SL_2(q), q odd prime power, from the generic character table (exact integer arithmetic).
# Rational characters:  DS(t) = sum of the discrete series X_theta, theta of order t | q+1, t > 2   (degree phi(t)(q-1)/2)
#                       HD    = eta_1 + eta_2   (degree q-1)        PS(t), t | q-1, t > 2  (degree phi(t)(q+1)/2)
#                       HP    = xi_1 + xi_2     (degree q+1)        St (degree q)
# For each: test every element of prime order (central involution, unipotent, split, non-split):  regular (trace 0) or free.
from sympy import isprime, factorint, totient, divisors, mobius
from math import gcd
def ram(t, m):                      # Ramanujan sum c_t(m)
    d = t // gcd(t, m)
    return mobius(d) * totient(t) // totient(d)
def status(n, tr, l):               # tr = rational trace on an element of prime order l (same on all its nontrivial powers)
    if tr == 0 and n % l == 0: return 'reg'
    if tr * (l - 1) == -n: return 'free'
    return None
rows = []
for q in range(3, 2000, 2):
    f = factorint(q)
    if len(f) != 1: continue
    p = list(f)[0]; e = f[p]
    eps = 1 if (q - 1) // 2 % 2 == 0 else -1
    size = q * (q * q - 1)
    cands = []
    for t in divisors(q + 1):
        if t > 2:
            c = int(totient(t)) // 2; n = c * (q - 1)
            thz = 1 if ((q + 1) // 2) % t == 0 else -1          # theta(-1)
            cands.append(('DS(%d)' % t, n, c, n * thz, lambda m, t=t: -int(ram(t, m)), lambda l: 0, -c, thz))
    for t in divisors(q - 1):
        if t > 2:
            c = int(totient(t)) // 2; n = c * (q + 1)
            alz = 1 if ((q - 1) // 2) % t == 0 else -1
            cands.append(('PS(%d)' % t, n, c, n * alz, lambda m: 0, lambda l, t=t: int(ram(t, l)), c, alz))
    sq = int(round(q ** 0.5)); issq = sq * sq == q
    if not issq:
        cands.append(('HD', q - 1, 2, -eps * (q - 1), lambda m: -2 * (-1) ** m, lambda l: 0, -1, -eps))
        cands.append(('HP', q + 1, 2, eps * (q + 1), lambda m: 0, lambda l: 2 * (-1) ** l, 1, eps))
    cands.append(('St', q, 1, q, lambda m: -1, lambda l: 1, 0, 1))
    for (name, n, c, trz, fb, fa, tru, zsign) in cands:
        if zsign != -1: continue                                  # faithful only
        if e > 1 and name in ('HD', 'HP'): continue               # prime powers: handled by hand in notes (not admissible)
        ok = True; free = []; reg = []
        st = status(n, trz, 2)
        if st is None: ok = False
        else: (free if st == 'free' else reg).append(2)
        st = status(n, tru, p)                                    # unipotent (q = p: all powers conjugate up to Galois)
        if e > 1: st = None if tru != 0 else st
        if st is None: ok = False
        else: (free if st == 'free' else reg).append(p)
        for l in factorint(q + 1):
            if l == 2: continue
            st = status(n, fb((q + 1) // l), l)
            if st is None: ok = False
            else: (free if st == 'free' else reg).append(l)
        for l in factorint(q - 1):
            if l == 2: continue
            st = status(n, fa((q - 1) // l), l)
            if st is None: ok = False
            else: (free if st == 'free' else reg).append(l)
        if ok:
            deg = n // c
            rows.append((q, name, n, c, size / (2 * c * deg * deg), sorted(set(free)), sorted(set(reg))))
print('q  char  dim(rational)  #conj  gain  free primes  regular primes')
for r in rows:
    if r[4] >= 2.5 or r[0] < 40: print(r[0], r[1], r[2], r[3], round(r[4], 3), r[5], r[6])
big = [r for r in rows if r[4] >= 2.5]
print('pairs with gain >= 2.5:', len(big), '; all are HD with q prime, q = 1 mod 4, (q+1)/2 prime:',
      all(r[1] == 'HD' and isprime(r[0]) and r[0] % 4 == 1 and isprime((r[0] + 1) // 2) for r in big if r[0] > 7))
print('primes q < 2000, q = 1 mod 4, (q+1)/2 prime:', [q for q in range(5, 2000, 4) if isprime(q) and isprime((q + 1) // 2)])
print('other admissible (gain < 2.5):', sorted(set((r[0], r[1]) for r in rows if r[4] < 2.5))[:40])
