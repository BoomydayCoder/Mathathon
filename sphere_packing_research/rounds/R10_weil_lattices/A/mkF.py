import sys
from weil import *
p = int(sys.argv[1]); q = int(sys.argv[2])
c = setup(p)
c.BasC = np.load(os.path.join(HERE, f'BasC_p{p}_L0.npy')); c.rank = c.n
finalize(c)
g = c.I.copy()
for _ in range(c.Hord // q): g = mm(g, c.Hgen)
Fq, sat = fixed_lattice(c, g, q)
dq = sub_det(c, Fq); rk = len(Fq); Dint = round(c.D)
log(f'C_{q} fixed lattice p={p}: rank {rk}, saturated {sat}; EXACT calibration det_B == (n/rk)^rk * D^(rk/2): {dq == (c.n // rk) ** rk * Dint ** (rk // 2)}')
np.save(os.path.join(HERE, f'F_p{p}_C{q}.npy'), Fq)
