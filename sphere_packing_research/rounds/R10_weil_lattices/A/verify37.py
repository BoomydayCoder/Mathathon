import sys, json
from weil import *
p = int(sys.argv[1])
c = setup(p)
OV, par, seen, key = orbit(c.V0, [c.S, c.T])
assert np.abs(OV - c.OV).max() < 1e-12
c.BasC = np.load(os.path.join(HERE, f'BasC_p{p}_L0.npy'))
c.rank = c.n
finalize(c)
assert (c.B == np.load(os.path.join(HERE, f'B_p{p}_L0.npy'))).all() and (c.E == np.load(os.path.join(HERE, f'E_p{p}_L0.npy'))).all()
Dint = round(c.D); assert Dint ** (p - 1) == c.detB
log(f'EXACT: det B = D^(p-1) with D = {Dint}')
bests = [json.load(open(os.path.join(HERE, f))) ['best'] for f in sorted(os.listdir(HERE)) if f.startswith(f'search_p{p}_proc')]
best = min(bests, key=lambda b: b['N_exact'])
x = np.array(best['x'], dtype=np.int64)
xo = np.array([int(z) for z in x], dtype=object)
b = int(xo.dot(c.B.astype(object)).dot(xo)); e = int(xo.dot(c.E.astype(object)).dot(xo))
N = p * b * b - e * e; THR = (p - 1) ** 2 * Dint
log(f'EXACT: b = x B x^T = {b}, e = x E x^T = {e}, N = p b^2 - e^2 = {N}, threshold (p-1)^2 D = {THR}, violator: {N < THR}; Gamma = N/(4pD) = {N/(4*p*Dint):.10f}, theta = {c.theta:.10f}, ratio {N/THR:.10f}')
# independent: from the embedding
w = np.tensordot(x.astype(float), c.BasC, axes=1)
rw = remb(w, c); mp = float(rw[:c.hp] @ rw[:c.hp]); mmn = float(rw[c.hp:] @ rw[c.hp:])
log(f'embedding check: m_plus(w) = {mp:.8f}, m_minus(w) = {mmn:.8f}, product/D = {mp*mmn/Dint:.10f}; from (b,e): {(b+e/math.sqrt(p))/2:.8f}, {(b-e/math.sqrt(p))/2:.8f}')
# direct construction of O.w via the commutant: delta_u(g v0) = g u
F, sat = fixed_lattice(c, c.Hgen, c.Hord)
FB = F.astype(object).dot(c.B.astype(object)).dot(F.astype(object).T)
dH = bareiss_det(FB)
x0 = np.rint(coords(c, c.V0[None])).astype(np.int64)[0]
b0 = int(x0 @ c.B @ x0); e0 = int(x0 @ c.E @ x0)
log(f'v0: B = {b0}, E = {e0} (m_plus = m_minus = {b0/2})')
log(f'EXACT order test: det_B(Lambda^H) * 4^4 * 2^8 == p^4 * b0^8 : {dH * 256 * 256 == p**4 * b0**8}   (i.e. det of Tr nrd form = (p/4)^4; valid since E(v0,v0) = {e0})')
r = (p - 1) // 4
log(f'EXACT calibration: det_B(Lambda^H) == r^8 D^4 : {dH == r**8 * Dint**4}')
U8 = np.tensordot(F.astype(float), c.BasC, axes=1)          # 8 x h x h, u_k = delta_k v0
X = remb(OV[:800], c)
Q, R, piv = scipy.linalg.qr(X.T, mode='economic', pivoting=True)
sel = piv[:c.n]
Mx = X[sel]
log(f'cond of chosen orbit basis: {np.linalg.cond(Mx):.2f}')
a = np.linalg.lstsq(Mx.T, rw, rcond=None)[0]
log(f'w in terms of orbit vectors: residual {np.abs(a @ Mx - rw).max():.2e}')
gens = [c.S, c.T]
def apply_g(i, Uarr):
    path = []
    while par[i][0] >= 0:
        path.append(par[i][1]); i = par[i][0]
    for gi in reversed(path):
        Uarr = act(gens[gi], Uarr)
    return Uarr
chk = max(np.abs(apply_g(int(i), c.V0[None])[0] - OV[i]).max() for i in sel)
log(f'path reconstruction check (g_i v0 = orbit vector): {chk:.2e}')
dW = np.zeros_like(U8)
for coef, i in zip(a, sel):
    dW += coef * apply_g(int(i), U8)
Cw = coords(c, dW)
log(f'O.w inside Lambda: coordinates integral, max frac = {np.abs(Cw-np.rint(Cw)).max():.2e}')
Y = np.rint(Cw).astype(np.int64)
GY = Y.astype(object).dot(c.B.astype(object)).dot(Y.astype(object).T)
dOw = bareiss_det(GY)
pred = (N * N) // 1  # ((p/4) m+ m-)^4 = (N/16)^4 / ... check below
log(f'EXACT: det_B of the rank-8 sublattice O.w = {dOw}; predicted ((p/4) m+ m-)^4 = (N/16)^4 = {N**4}/{16**4} -> equal: {dOw * 16**4 == N**4}')
log(f'EXACT: r^8 D^4 = {r**8 * Dint**4};  R_U >= r^8 D^4 / det_B(O.w) = {r**8 * Dint**4 / dOw:.8f}  (>1 means unstable)')
# is w (and x) in the span? rank of [Y; x]
rk = np.linalg.matrix_rank(np.vstack([Y, x]).astype(float))
log(f'rank of O.w = {np.linalg.matrix_rank(Y.astype(float))}, with w added: {rk}')
# G-equivariance spot check of delta: delta_k(S w) = S delta_k(w)
Sw = act(c.S, w); rs = remb(Sw, c); a2 = np.linalg.lstsq(Mx.T, rs, rcond=None)[0]
dSW = np.zeros_like(U8)
for coef, i in zip(a2, sel):
    dSW += coef * apply_g(int(i), U8)
log(f'equivariance check |delta(S w) - S delta(w)| = {np.abs(dSW - act(c.S, dW)).max():.2e}')
cert = dict(p=p, basis_files=[f'B_p{p}_L0.npy', f'E_p{p}_L0.npy', f'BasC_p{p}_L0.npy'], x=[int(z) for z in x], B_norm=b, E_norm=e,
            m_plus=(b + e / math.sqrt(p)) / 2, m_minus=(b - e / math.sqrt(p)) / 2, D=Dint, detB=str(c.detB), Gamma=N / (4 * p * Dint), theta=c.theta,
            N_exact=N, THR_exact=THR, violator=bool(N < THR), R_U_lower_bound=r**8 * Dint**4 / dOw, det_B_Ow=str(dOw), Ow_rows=Y.tolist(),
            formula='Gamma = (p b^2 - e^2)/(4 p D); violator iff p b^2 - e^2 < (p-1)^2 D; Ow_rows = Z-basis (integer coords) of the rank-8 sublattice O.w of Lambda')
json.dump(cert, open(os.path.join(HERE, f'violator_p{p}.json'), 'w'), indent=1)
log('certificate written')
