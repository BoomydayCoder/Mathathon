import itertools, math
from fractions import Fraction
from scipy.stats import beta
# (1) p=13 enumeration: a_j = w.g^j w in {-2..2}, j=1..6, sum=-2, circulant PSD  -> possible s
p=13; QR={pow(x,2,p) for x in range(1,p)}
leg=lambda j: 1 if j%p in QR else -1
cos=[[math.cos(2*math.pi*j*k/p) for j in range(1,7)] for k in range(1,7)]
vals={}
for a in itertools.product(range(-2,3),repeat=6):
    if sum(a)!=-2: continue
    if all(4+2*sum(a[j]*cos[k][j] for j in range(6))>-1e-9 for k in range(6)):
        s=sum(leg(j+1)*a[j] for j in range(6)); vals[s]=vals.get(s,0)+1
print("p=13: admissible s values and number of sequences:",sorted(vals.items()))
# (2) heuristic number of minimal vectors violating N>=r^2 (Beta model), and variance from the 4-design identity E s^2 = mu^2/4
for p,mu,K,name in [(13,4,196560,"Leech"),(17,4,146880,"extremal 32"),(29,6,15590400,"extremal 56"),(37,8,6218175600,"extremal 72"),(37,6,None,"min 6 in dim 72"),(41,8,1250172000,"extremal 80"),(61,12,None,"extremal 120 (hyp.)"),(73,14,None,"extremal 144 (hyp.)")]:
    r2=Fraction(p-1,4)**2
    # N = (p mu^2/4 - s^2)/4 >= r2  <=> s^2 <= p mu^2/4 - 4 r2 ; s = mu/2 mod 2
    bound=Fraction(p*mu*mu,4)-4*r2
    smax=max([s for s in range(0,10**4) if s*s<=bound and (s-mu//2)%2==0],default=None)
    allowed=[s for s in range(-10**3,10**3) if s*s<=bound and (s-mu//2)%2==0]
    d=(p-1)/2
    if smax is None: print(p,mu,name,"NO admissible s at all: every minimal vector violates"); continue
    # s = (sqrt p mu/2)(2x-1), x~Beta(d,d); violation |s| > smax  (next admissible value smax+2; use threshold smax+1)
    y=(smax+1)/(math.sqrt(p)*mu/2)
    tail=2*beta.sf((1+y)/2,d,d) if y<1 else 0.0
    print("p=%d mu=%d %s: need N>=%s, i.e. |s|<=%d (allowed s: %s); design variance E s^2 = %.2f ; Beta tail P(|s|>%d)=%.3e ; expected # violating minimal vectors = %s"%(
        p,mu,name,r2,smax,allowed if len(allowed)<12 else "...",mu*mu/4,smax,tail, "%.3e"%(K*tail) if K else "n/a"))
# (3) Lemma A' numbers: R_U lower bounds for fixed spaces of elements of odd prime order q | (p+1)/2
for p in [17,29,41,53,89,97,101,109,113]:
    P=(p+1)//2; out=[]
    for q in range(3,P+1,2):
        if P%q==0 and all(q%t for t in range(3,int(q**.5)+1,2)):
            N=P//q
            if N>1:
                a=4*(N-1); m=2*(p-1); out.append((q,a,round((m/(a*q))**a,3)))
    print("p=%d (p+1)/2=%d: (q, dim fixed space, lower bound for R_U) ="%(p,P),out)
# (4) Eichler mass (1/2) zeta_F(-1) for the unramified totally definite quaternion algebra over Q(sqrt p)
def sigma(n): return sum(d for d in range(1,n+1) if n%d==0)
for p in [5,13,17,29,37,41,53,61,73]:
    z=Fraction(sum(sigma((p-b*b)//4) for b in range(-int(p**.5),int(p**.5)+1) if b%2 and b*b<p),60)
    print("p=%d zeta_F(-1)=%s mass=%s"%(p,z,z/2))
# (5) constants: C = |G| D /(2e n^2) with D=2k(2k-1), n=2(p-1)k -> p(p+1)/(2e(p-1))
for p in [5,13,37,61,73]:
    print("p=%d |G|=%d  C_infty = p(p+1)/((p-1)2e) = %.3f  ratio to 1/(2e): %.2f"%(p,p*(p*p-1),p*(p+1)/((p-1)*2*math.e),p*(p+1)/(p-1)))
