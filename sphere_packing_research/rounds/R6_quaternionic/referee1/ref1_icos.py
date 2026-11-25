import numpy as np, math, itertools
rng = np.random.default_rng(3)
def qmul(p,q):
    a,b,c,d = p; e,f,g,h = q
    return np.array([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g, a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])
def Rq(q):
    a,b,c,d = q
    return np.array([[a,-b,-c,-d],[b,a,d,-c],[c,-d,a,b],[d,c,-b,a]],float)
# check Rq is right multiplication
x = rng.normal(size=4); u = rng.normal(size=4)
assert np.allclose(Rq(u)@x, qmul(x,u))
def gens(tau):
    return [np.array([0,1.,0,0]), np.array([.5,.5,.5,.5]), np.array([tau/2, (1/tau)/2, .5, 0])]
t1 = (1+math.sqrt(5))/2; t2 = (1-math.sqrt(5))/2
G = [(np.array([1.,0,0,0]), np.array([1.,0,0,0]))]
g = list(zip(gens(t1), gens(t2)))
def key(p): return tuple(np.round(np.concatenate(p),8)+0.0)
seen = {key(G[0])}
front = list(G)
while front:
    new = []
    for (a1,a2) in front:
        for (b1,b2) in g:
            c = (qmul(a1,b1), qmul(a2,b2)); kk = key(c)
            if kk not in seen:
                seen.add(kk); G.append(c); new.append(c)
    front = new
print("|2I| =", len(G), " all unit norms:", all(abs(np.dot(a,a)-1)<1e-9 and abs(np.dot(b,b)-1)<1e-9 for a,b in G))
for k in (1,2):
    n = 8*k
    Rs = []
    for (u1,u2) in G:
        R = np.zeros((n,n)); R[:4*k,:4*k] = np.kron(np.eye(k),Rq(u1)); R[4*k:,4*k:] = np.kron(np.eye(k),Rq(u2)); Rs.append(R)
    def P(B): return sum(R@B@R.T for R in Rs)/len(Rs)
    for _ in range(3):
        x = rng.normal(size=n)*np.concatenate([np.ones(4*k)*rng.uniform(.2,2), np.ones(4*k)])
        PB = P(np.outer(x,x)); q = np.sum(PB*PB)/np.dot(x,x)**2
        u1 = np.dot(x[:4*k],x[:4*k])/np.dot(x,x); u2 = 1-u1
        print(f" k={k}: q={q:.12f}  (u1^2+u2^2)/4={(u1*u1+u2*u2)/4:.12f}  >= 1/8: {q>=0.125-1e-12}")
    basis = []
    for i in range(n):
        for j in range(i,n):
            E = np.zeros((n,n)); E[i,j]=E[j,i]=1; basis.append(P(E).ravel())
    print(f" k={k}: dim Q = {np.linalg.matrix_rank(np.array(basis),tol=1e-8)} expected {2*k*(2*k-1)}")
# free action: no nontrivial u fixes a nonzero vector  <=> R_u - I invertible
k=1
mins = []
for (u1,u2) in G[1:]:
    R = np.zeros((8,8)); R[:4,:4]=Rq(u1); R[4:,4:]=Rq(u2)
    mins.append(np.linalg.svd(R-np.eye(8),compute_uv=False).min())
print(" min singular value of R_u - I over u != 1:", min(mins))
# angular lemma check s=2
from scipy import special, integrate
def lhs(a,b,s=2):
    f = lambda u: math.exp((b*s/2)*((u-.5)**2+(1-u-.5)**2) + (a-1)*math.log(u)+(a-1)*math.log(1-u) - special.betaln(a,a))
    return integrate.quad(f,0,1,limit=400, points=[.25,.5,.75])[0]
for (a,b) in ((20,4),(100,30),(400,60),(2000,56)):
    print(f" Dirichlet a={a} b={b}: E exp = {lhs(a,b):.6f} <= bound {(1-b/a)**(-2*b):.6f}")
