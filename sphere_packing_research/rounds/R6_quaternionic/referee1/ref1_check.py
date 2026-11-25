# Referee 1 independent checks (analysis + constants) for R6 proof.md
import numpy as np, itertools, math
from scipy import integrate

rng = np.random.default_rng(7)

def Rq(q):
    a,b,c,d = q
    return np.array([[a,-b,-c,-d],[b,a,d,-c],[c,-d,a,b],[d,c,-b,a]],float)

# --- (1) Hurwitz units, right-multiplication matrices, P(x⊗x), q, dim Q
units = []
for s in (1,-1):
    for i in range(4):
        e = [0,0,0,0]; e[i] = s; units.append(tuple(e))
for sg in itertools.product((0.5,-0.5), repeat=4):
    units.append(sg)
assert len(units) == 24
for k in (1,2,3):
    n = 4*k
    Rs = [np.kron(np.eye(k), Rq(u)) for u in units]
    for R in Rs: assert np.allclose(R@R.T, np.eye(n))
    def P(B): return sum(R@B@R.T for R in Rs)/24
    qs = []
    for _ in range(5):
        x = rng.normal(size=n); B = np.outer(x,x); PB = P(B)
        qs.append(np.sum(PB*PB)/np.dot(x,x)**2)
    # dimension of Q: rank of P on a basis of Sym_n
    basis = []
    for i in range(n):
        for j in range(i,n):
            E = np.zeros((n,n)); E[i,j]=E[j,i]=1; basis.append(P(E).ravel())
    rk = np.linalg.matrix_rank(np.array(basis), tol=1e-9)
    print(f"Hurwitz k={k}: q values {np.round(qs,12)}  dimQ={rk} expected {k*(2*k-1)}")

# --- (2) Ito drift: sum_{i,j} |pi(u_i ⊗_s u_j)|^2 = Tr pi  (pi = P here), k=2
k=2; n=8
Rs = [np.kron(np.eye(k), Rq(u)) for u in units]
def P(B): return sum(R@B@R.T for R in Rs)/24
U,_ = np.linalg.qr(rng.normal(size=(n,n)))
tot = 0
for i in range(n):
    for j in range(n):
        S = (np.outer(U[:,i],U[:,j])+np.outer(U[:,j],U[:,i]))/2
        PS = P(S); tot += np.sum(PS*PS)
print("sum |P(u_i (x)_s u_j)|^2 =", tot, " expected D_k =", k*(2*k-1))

# --- (3) shell integral (Prop 3.5) for constant q, direct numerical evaluation
def Phi_log(y):
    # log of min(1/2, exp(-y^2/2)/(sqrt(2pi) y))
    if y <= 0: return math.log(0.5)
    return min(math.log(0.5), -y*y/2 - math.log(math.sqrt(2*math.pi)*y))

def shell(n, q, w, C0=2.0):
    rho = math.log(math.log(n)); alpha = (1-rho/n)**-2
    z = math.log(w*q*n*n/4) + rho; T = 8*z/(q*n*n)
    m = (n+2)/2
    Y = C0*math.sqrt(n/q)
    def inner(t):
        s = math.sqrt(q*t)
        f = lambda y: math.exp(Phi_log(y) - m*math.log1p(-y*s/alpha))
        peak = n*s/2
        pts = sorted(set([1.0, max(1.0,peak-8), peak, peak+8]))
        val,_ = integrate.quad(f, 0, min(Y, peak+60), points=[p for p in pts if p < min(Y,peak+60)], limit=400)
        return n*s/2 * alpha**(-m) * val
    val,_ = integrate.quad(inner, 0, T, limit=400, points=[T*0.5,T*0.8,T*0.9,T*0.97])
    bound = math.exp(-rho)*8/(q*n*n)*math.exp(z)
    return val, bound, val/bound, z, T

for n in (10**3, 10**4, 10**5, 10**6):
    for (q,w,name) in ((0.25,24,'Hurwitz'),(1.0,2,'trivial')):
        val,bound,ratio,z,T = shell(n,q,w)
        print(f"n={n:>8} {name:8s} shell integral={val:.5g} claimed bound={bound:.5g} ratio={ratio:.4f} (z={z:.2f}, T n^2={T*n*n:.1f})")

# --- (4) final optimisation: minimise 2rho - theta z + (8 e^{-rho}/(w q0 n^2)) e^z
for (w,delta,name) in ((24,4,'Hurwitz'),(120,8,'icosian'),(2,1,'trivial')):
    q0 = 1/delta
    n = 10**6; rho = math.log(math.log(n))
    zs = np.linspace(5, 60, 2000001)
    f = 2*rho - 2*zs + 8*math.exp(-rho)/(w*q0*n*n)*np.exp(zs)
    fmin = f.min()
    V = math.exp(-fmin/2)   # det^{-1/2}
    print(f"{name}: numeric min -> V/n^2 = {V/n**2:.6f}; formula w/(4 e delta) = {w/(4*math.e*delta):.6f}")
print("3/(2e) =", 3/(2*math.e), " 15/(4e) =", 15/(4*math.e), " 1/(2e) =", 1/(2*math.e))
