import sys, json
from weil import *
p = int(sys.argv[1])
rng = np.random.default_rng(1)
c = setup(p)
build_lattice(c, rng)
np.save(os.path.join(HERE, f'B_p{p}_L0.npy'), c.B); np.save(os.path.join(HERE, f'E_p{p}_L0.npy'), c.E)
np.save(os.path.join(HERE, f'BasC_p{p}_L0.npy'), c.BasC)
res = dict(p=p, rank=c.rank, detB=str(c.detB), D=c.D, theta=c.theta, detfac=c.detfac)
# order check
F, sat = fixed_lattice(c, c.Hgen, c.Hord)
log(f'fixed lattice of H (order {c.Hord}): rank {len(F)}, saturated {sat}')
x0 = np.rint(coords(c, c.V0[None])).astype(np.int64)
g0, mp0, mm0, b0, e0 = gamma_int(c, x0)
log(f'v0: m_plus={mp0[0]:.6f} m_minus={mm0[0]:.6f} B={b0[0]} Gamma(v0)={g0[0]:.6f}')
FG = (F @ c.Gp @ F.T) / mp0[0] + (F @ c.Gm @ F.T) / mm0[0]
dO = np.linalg.det(FG)
log(f'order test: det(Tr nrd form on Lambda^H) = {dO:.8f}; maximal iff = (p/4)^4 = {(p/4)**4:.8f}; ratio = {dO/(p/4)**4:.8f}')
res['order_det_ratio'] = dO / (p / 4) ** 4
dH = sub_det(c, F)
r = (p - 1) // 4
if len(F) == 8:
    target = r ** 8 * c.D ** 4
    log(f'calibration H-line: det_B(Lambda^H) = {dH}; r^8 D^4 = {target:.6f}; R_U = {target/dH:.8f}')
    res['calib_RU_H'] = target / dH
out = {}
if len(F) <= 24:
    out['H'] = grid_enum(c, F, label=f'p={p} H-fixed (order {c.Hord})')
if p % 8 == 5 and p > 5:
    m = c.Hord
    for q in [d for d in range(2, m) if m % d == 0]:
        g = c.I.copy()
        for _ in range(m // q): g = mm(g, c.Hgen)
        Fq, satq = fixed_lattice(c, g, q)
        dq = sub_det(c, Fq); rk = len(Fq)
        RU = (c.n / rk) ** rk * math.exp(math.log(c.detB) * rk / c.n) / dq
        log(f'C_{q} fixed lattice: rank {rk}, saturated {satq}, det_B = {dq}, R_U = {RU:.8f}')
        res[f'calib_RU_C{q}'] = RU; res[f'rank_C{q}'] = rk
        if rk <= 40:
            out[f'C{q}'] = grid_enum(c, Fq, label=f'p={p} C_{q}-fixed', tl=300)
            np.save(os.path.join(HERE, f'F_p{p}_C{q}.npy'), Fq)
np.save(os.path.join(HERE, f'F_p{p}_H.npy'), F)
if c.n <= 24:
    I = np.eye(c.n, dtype=np.int64)
    out['full'] = grid_enum(c, I, label=f'p={p} FULL lattice', tl=600)
    # minimum of B itself
    Y, fl, nd = enum_short(c.B.astype(float), b0[0] - 0.5, maxnodes=3_000_000_000)
    log(f'vectors with B-norm < B(v0,v0)={b0[0]}: {len(Y)} (flag {fl}, nodes {nd})')
    res['n_shorter_than_v0'] = len(Y); res['enum_flag_min'] = fl
for k, o in out.items():
    res[k] = dict(min_gamma=o['min_gamma'], complete=o['complete'], n_enum=o['n_enum'], n_viol=o['n_viol'], rank=o['rank'])
    if o['best'][1] is not None:
        res[k]['cert'] = cert(c, o['best'])
json.dump(res, open(os.path.join(HERE, f'result_p{p}_L0.json'), 'w'), indent=1)
log('done')
