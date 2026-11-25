import numpy as np, math
from scipy import integrate, stats, special
rng = np.random.default_rng(1)
print("== (a) Siegel cone volume n=2: expect zeta(2)/2 =", math.pi**2/12)
# F = reduced bases b1=(r,0) rotated by theta in [0,pi), b2=(x,y), |x|<=r/2, x^2+y^2>=r^2, 0<r*y<=1
f = lambda s: 0.25/math.sqrt(1-s*s)
val = math.pi*integrate.quad(f,-0.5,0.5)[0]
print("   analytic reduction:", val)
# brute-force Monte Carlo of the same 4-dim cone (independent of my algebra)
N=4_000_000
r = rng.uniform(0,1.3,N); x = rng.uniform(-0.65,0.65,N); y = rng.uniform(0,6,N)
ok = (np.abs(x)<=r/2)&(x*x+y*y>=r*r)&(r*y<=1)
vol = math.pi*np.mean(ok*r)*1.3*1.3*6      # measure r dr dtheta dx dy, theta in [0,pi)
print("   Monte Carlo (y truncated at 6):", vol)
print("== (b) Lemma N exact bound")
def lemmaN(n, eps=None):
    eps = 1/(2*n) if eps is None else eps
    N2=n*n
    logk = (N2/2)*math.log(math.pi)-math.lgamma(N2/2+1)
    lz = sum(math.log(float(special.zeta(k))) for k in range(2,n+1))
    return math.log(n)+N2*math.log(eps/(2+2*eps))+logk-lz
for n in (2,3,8,24,100):
    print("   n=%d exact %.1f  simplified %.1f  (3/4 n^2 log n = %.1f)"%(n,lemmaN(n),-2*n*n*math.log(n)-n*n, -0.75*n*n*math.log(n)))
print("== (c) spherical identity")
n=7; D=n*(n+1)//2
def sym_basis_proj(vecs):
    # orthonormal basis (Frobenius) of span of symmetric matrices
    M=np.array([v.reshape(-1) for v in vecs]); q,_=np.linalg.qr(M.T); return q
for c in (1,5,12,20,27):
    ws = rng.normal(size=(c,n)); ws/=np.linalg.norm(ws,axis=1)[:,None]
    Q = sym_basis_proj([np.outer(w,w) for w in ws])
    Id=np.eye(n).reshape(-1); PId2=np.sum((Q.T@Id)**2)
    w = rng.normal(size=(200000,n)); w/=np.linalg.norm(w,axis=1)[:,None]
    T = (w[:,:,None]*w[:,None,:]).reshape(len(w),-1)
    mc = np.mean(np.sum((T@Q)**2,axis=1))
    print("   c=%2d  MC %.5f  formula %.5f   |P_E Id|^2=%.3f vs cn/(n-1+c)=%.3f ; 1-E = %.4f vs N/D = %.4f"%(c,mc,(2*c+PId2)/(n*(n+2)),PId2,c*n/(n-1+c),1-mc,1-c/D))
print("== (d) contact count closed form vs quadrature (a0=1, V0=1)")
for n,u in ((1000,0.5),(1000,2.0),(10**5,2.0)):
    tau=8*u/n**2
    g=lambda s: (n/2)*(1-s)**(-n/2-1)*2*stats.norm.sf(s/math.sqrt(tau))
    q=integrate.quad(g,0,min(0.9,40*math.sqrt(tau)),limit=400)[0]/2
    print("   n=%d u=%.1f  pairs quad %.4f  closed %.4f"%(n,u,q,math.exp(u)*stats.norm.cdf(math.sqrt(2*u))-0.5))
print("== (e) model K / R constants")
V0=1/math.e
print("   model K at T_K: %.4f n^2 ; max D/e -> %.4f n^2"%(V0*math.exp(-2*V0), 1/(2*math.e)))
print("   model R at T_K: V0=1: %.4f, V0=1/e: %.4f"%(0.5/(1+0.5),0.5/(1+0.5/V0)))
print("== (f) Gamma/Poisson tail vs Chernoff")
for Dd,lam in ((55,2.0),(210,1.5),(210,3.0)):
    p=stats.gamma.sf(lam*Dd,Dd+1); ch=math.exp(-Dd*(lam-1-math.log(lam)))
    print("   D=%d lam=%.1f  P(Gamma(D+1)>=lam D)=%.3e  P(Poi<=D)=%.3e  Chernoff %.3e"%(Dd,lam,p,stats.poisson.cdf(Dd,lam*Dd),ch))
print("== (g) Prop U: det bound (1+1/n)^n<=e and threshold constants")
for n in (10,100,1000):
    D=n*(n+1)/2
    for rho in (1/n, math.sqrt(math.log(n)/n), 0.5):
        thr = 2*math.sqrt(math.e)*D*math.log(1+4*rho*n/(1-rho))
        print("   n=%d rho=%.3g: Poisson first-moment threshold of Prop U  V = %.3f n^2 log n"%(n,rho,thr/(n*n*math.log(n))))
