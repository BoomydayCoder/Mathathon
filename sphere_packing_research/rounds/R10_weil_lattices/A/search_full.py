import sys, json
from weil import *
p = int(sys.argv[1]); off = int(sys.argv[2]); nproc = int(sys.argv[3]); beta1 = int(sys.argv[4]); beta2 = int(sys.argv[5]); tlim = float(sys.argv[6])
B = np.load(os.path.join(HERE, f'B_p{p}_L0.npy')); E = np.load(os.path.join(HERE, f'E_p{p}_L0.npy'))
n = len(B); sp = math.sqrt(p)
detB = bareiss_det(B); D = round(math.exp(math.log(detB) / (p - 1)))
assert D ** (p - 1) == detB, 'D not an exact integer root'
theta = (p - 1) ** 2 / (4 * p); THR = (p - 1) ** 2 * D   # violation iff N = p b^2 - e^2 < THR (exact integers)
Gp = (B + E / sp) / 2; Gm = (B - E / sp) / 2
Bo = B.astype(object); Eo = E.astype(object)
def Nex(x):
    xo = np.array([int(z) for z in x], dtype=object)
    b = int(xo.dot(Bo).dot(xo)); e = int(xo.dot(Eo).dot(xo))
    return p * b * b - e * e, b, e
a_, b_, eps = fund_unit(p); le = math.log(eps)
ts = np.arange(-le, le + 0.1, 0.1)
ts = ts[off::nproc]
best = None; out = os.path.join(HERE, f'search_p{p}_proc{off}.json'); rec = []
t0 = time.time()
rng = np.random.default_rng(100 + off)
for t in ts:
    if time.time() - t0 > tlim: break
    A = math.exp(t) * Gp + math.exp(-t) * Gm
    U = None
    for beta in (beta1, beta2):
        if time.time() - t0 > tlim: break
        U, r, nt = c_bkz(A, beta, tours=3, maxnodes=3_000_000, U0=U)
        V = U.copy()
        # candidates: basis rows and pairwise sums/differences of the 40 shortest
        q = np.einsum('ij,jk,ik->i', V.astype(float), A, V.astype(float))
        order = np.argsort(q)[:40]
        cands = [V[i] for i in range(n)]
        for i in range(len(order)):
            for j in range(i + 1, len(order)):
                cands.append(V[order[i]] + V[order[j]]); cands.append(V[order[i]] - V[order[j]])
        Cn = np.array(cands)
        bb = np.einsum('ij,jk,ik->i', Cn.astype(float), B.astype(float), Cn.astype(float))
        ee = np.einsum('ij,jk,ik->i', Cn.astype(float), E.astype(float), Cn.astype(float))
        gam = (p * bb * bb - ee * ee) / (4 * p * D)
        i = int(np.argmin(gam))
        N, b, e = Nex(Cn[i])
        g = N / (4 * p * D)
        log(f'p={p} t={t:+.3f} beta={beta} tours={nt}: |b1|^2_Qt={q.min():.3f} (2sqrt(theta D)={2*math.sqrt(theta*D):.3f}), min Gamma over {len(Cn)} candidates = {g:.6f} (theta {theta:.6f}) N={N} THR={THR} violator={N < THR}')
        rec.append(dict(t=float(t), beta=beta, b1=float(q.min()), min_gamma=g, N=N, violator=bool(N < THR)))
        if best is None or N < best['N_exact']:
            best = dict(p=p, x=[int(z) for z in Cn[i]], B_norm=b, E_norm=e, N_exact=N, THR_exact=THR, D=D, Gamma=g, theta=theta,
                        m_plus=(b + e / sp) / 2, m_minus=(b - e / sp) / 2, t=float(t), beta=beta, violator=bool(N < THR))
        json.dump(dict(best=best, rec=rec), open(out, 'w'), indent=1)
log('search finished', 'best Gamma', None if best is None else best['Gamma'])
