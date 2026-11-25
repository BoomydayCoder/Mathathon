import numpy as np, itertools, math
from math import gcd
rng=np.random.default_rng(1)
def phi(m): return sum(1 for a in range(1,m+1) if gcd(a,m)==1)
# (i) cyclic: sum_j c_m(j)^2/phi(m)^2 = m/phi(m)
for m in [3,4,5,6,7,8,9,12,15,30,210]:
    c=[sum(math.cos(2*math.pi*a*j/m) for a in range(1,m+1) if gcd(a,m)==1) for j in range(m)]
    print("cyclic m=%d: sum rho^2=%.6f  m/phi=%.6f"%(m,sum(x*x for x in c)/phi(m)**2,m/phi(m)))
# quaternion groups: elements as (a,b,c,d); generate closure
def qmul(p,q):
    a,b,c,d=p; e,f,g,h=q
    return (a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g, a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e)
def closure(gens):
    S=[(1.0,0,0,0)]
    def key(q): return tuple(round(x,6)+0.0 for x in q)
    seen={key(S[0])}; frontier=list(S)
    while frontier:
        new=[]
        for x in frontier:
            for g in gens:
                y=qmul(x,g); k=key(y)
                if k not in seen: seen.add(k); new.append(y)
        S+=new; frontier=new
    return S
i=(0,1.0,0,0); j=(0,0,1.0,0)
om=(-0.5,0.5,0.5,0.5)
tau=(1+5**.5)/2
ico=(0.5*tau, 0.5*(tau-1), 0.5, 0.0)   # (tau + sigma*i + j)/2 with sigma=tau-1=1/tau, an icosian unit of order 10
oct_=(2**-.5,2**-.5,0,0)
def order(q):
    x=(1.0,0,0,0); n=0
    while True:
        x=qmul(x,q); n+=1
        if abs(x[0]-1)<1e-6 and abs(x[1])+abs(x[2])+abs(x[3])<1e-6: return n
def mu2_over_phi2(o):
    # mu(o)^2/phi(o)^2
    n=o; p=2; sqfree=True
    while p*p<=n:
        if n%p==0:
            n//=p
            if n%p==0: sqfree=False
        p+=1
    return (1.0/phi(o)**2) if sqfree else 0.0
for name,gens,dimQ in [("Q8",[i,j],4),("2T",[i,j,om],4),("2O",[i,j,om,oct_],8),("2I",[i,om,ico],8)]:
    G=closure(gens)
    s=sum(mu2_over_phi2(order(g)) for g in G)
    print(name,"|G|=",len(G)," sum mu^2/phi^2(o(g))=%.6f"%s," |G|/dimQ(A)=%.6f"%(len(G)/dimQ)," gain over Klartag=|G|/(2dim)=%.4f"%(len(G)/(2*dimQ)))
# (ii)-(iii) cyclotomic: real rep of mu_m on (C^k)^s, s=phi(m)/2
def cyc_rep(m,k):
    units=[a for a in range(1,m) if gcd(a,m)==1 and a<m-a] if m>2 else [1]
    s=len(units); n=2*k*s
    def mat(jj):
        M=np.zeros((n,n))
        for si,a in enumerate(units):
            th=2*math.pi*a*jj/m; R=np.array([[math.cos(th),-math.sin(th)],[math.sin(th),math.cos(th)]])
            for kk in range(k):
                o=2*(si*k+kk); M[o:o+2,o:o+2]=R
        return M
    return [mat(jj) for jj in range(m)],s,n
for m,k in [(5,2),(7,2),(8,3),(12,3),(9,2),(15,1),(7,4)]:
    G,s,n=cyc_rep(m,k)
    # D = (1/|G|) sum (tr(g)^2+tr(g^2))/2
    D=sum((np.trace(g)**2+np.trace(g@g))/2 for g in G)/len(G)
    # MC of q
    X=rng.standard_normal((4000,n)); X/=np.linalg.norm(X,axis=1)[:,None]
    q_direct=np.mean([np.einsum('ij,ij->i',X@g.T,X)**2 for g in G],axis=0)
    blocks=(X.reshape(4000,s,2*k)**2).sum(axis=2)
    q_formula=0.5*(blocks**2).sum(axis=1)
    print("m=%d k=%d n=%d: D=%.3f (k^2 phi/2=%.1f); max|q_direct-q_formula|=%.2e; mean q=%.5f, (k+1)/(n+2)=%.5f, (2D+n)/(n(n+2))=%.5f; min q seen=%.5f >= 1/phi=%.5f"%(m,k,n,D,k*k*phi(m)/2,np.abs(q_direct-q_formula).max(),q_formula.mean(),(k+1)/(n+2),(2*D+n)/(n*(n+2)),q_formula.min(),1/phi(m)))
