import sys, json
from weil import *
p = int(sys.argv[1]); Ffile = sys.argv[2] if len(sys.argv) > 2 else None
c = Ctx(); c.p = p
c.B = np.load(os.path.join(HERE, f'B_p{p}_L0.npy')); c.E = np.load(os.path.join(HERE, f'E_p{p}_L0.npy'))
c.n = len(c.B); sp = math.sqrt(p)
c.Gp = (c.B + c.E / sp) / 2; c.Gm = (c.B - c.E / sp) / 2
detB = bareiss_det(c.B); D = round(math.exp(math.log(detB) / (p - 1))); assert D ** (p - 1) == detB
c.D = float(D); c.theta = (p - 1) ** 2 / (4 * p)
F = np.eye(c.n, dtype=np.int64) if Ffile is None else np.load(os.path.join(HERE, Ffile))
a, b, eps = fund_unit(p); le = math.log(eps); step = 0.1
ts = np.arange(-le, le + step, step)
bound0 = 2 * math.sqrt(c.theta * c.D) * math.cosh(step / 2) * (1 + 1e-9)
GpF = F @ c.Gp @ F.T; GmF = F @ c.Gm @ F.T
res = []; THR = (p - 1) ** 2 * D; Bo = c.B.astype(object); Eo = c.E.astype(object)
minN = None; tot = 0; viol = 0; allfl = 0; nodes = 0
for t in ts:
    A = math.exp(t) * GpF + math.exp(-t) * GmF
    Y, fl, nd = enum_short(A, bound0, maxc=300000, maxnodes=int(sys.argv[3]) if len(sys.argv) > 3 else 2_000_000_000)
    allfl |= fl; nodes += nd; tot += len(Y)
    mN = None
    for y in Y:
        xo = np.array([int(z) for z in (y @ F)], dtype=object)
        bb = int(xo.dot(Bo).dot(xo)); ee = int(xo.dot(Eo).dot(xo)); N = p * bb * bb - ee * ee
        assert bb > 0
        if mN is None or N < mN: mN = N
        if N < THR: viol += 1
    if mN is not None and (minN is None or mN < minN): minN = mN
    log(f't={t:+.3f}: {len(Y)} vectors, flag {fl}, nodes {nd}, min N = {mN} (THR {THR})')
log(f'SUMMARY p={p} F={Ffile} rank {len(F)}: t-values {len(ts)}, total vectors {tot}, nodes {nodes}, flags {allfl}, exact min N = {minN}, THR = (p-1)^2 D = {THR}, min Gamma/theta = {None if minN is None else minN/THR}, violators {viol}')
