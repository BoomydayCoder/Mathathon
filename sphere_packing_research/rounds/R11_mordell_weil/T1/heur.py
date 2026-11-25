# Heuristic mean-value counts for violating Delta-lines and Delta-planes (generic-lattice model).
from mpmath import mp, mpf, pi, gamma, log, sqrt, e
mp.dps = 40
def funit(p):
    # fundamental unit of Q(sqrt p), p = 1 mod 4: smallest (a+b sqrt p)/2 >1 with a^2 - p b^2 = -4 (norm -1)
    b = 1
    while True:
        for s in (-4, 4):
            a2 = p*b*b + s
            a = int(mp.sqrt(a2))
            for aa in (a-1, a, a+1):
                if aa > 0 and aa*aa == a2:
                    return (aa + b*sqrt(p))/2, s
        b += 1
def omega(h): return pi**(mpf(h)/2)/gamma(mpf(h)/2+1)
def A1(s, p, le):   # number of vectors mod eps^Z with N(q_F) < N0 in a rank-s O-lattice of covol 1: A1 * N0^(2s)
    h = 4*s
    return h*omega(h)**2*(mpf(4)/p)**(2*s)*le
print("p r  log(eps)  N1_exact  N1_stirling  |G|  vec-orbits=N1/|G|   planes P2 (u=2)  plane-orbits=P2/(|G|/2)")
for p in (5, 13, 17, 29, 37, 41, 61, 73, 157, 193):
    r = (p-1)//4
    eps, s = funit(p)
    le = log(eps)
    G = p*(p*p-1)
    N1 = A1(r, p, le)*mpf(r)**(4*r)          # threshold N(q_F) < r^2, i.e. N0^(2r) = r^(4r)
    N1s = le/pi*(pi*e/sqrt(p))**(p-1)
    out = "%4d %3d %8.4f %12.4e %12.4e %9d %12.4e" % (p, r, le, N1, N1s, G, N1/G)
    if r >= 3:
        u = 2
        A2 = (r-1)*A1(r, p, le)*A1(r-1, p, le)/(2*u*A1(2, p, le))
        P2 = A2*(mpf(r)/2)**(8*r)            # threshold N(Mdet H) < (r/2)^4, count ~ A2*Y^(2r)
        out += "  %12.4e %12.4e" % (P2, P2/(mpf(G)/2))
    print(out, " norm(eps)=", -1 if s == -4 else 1)
