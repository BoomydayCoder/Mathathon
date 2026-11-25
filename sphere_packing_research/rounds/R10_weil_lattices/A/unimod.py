import sys
from weil import *
p = int(sys.argv[1]); mult = float(sys.argv[2])
B = np.load(os.path.join(HERE, f'B_p{p}_L0.npy')); E = np.load(os.path.join(HERE, f'E_p{p}_L0.npy'))
sp = math.sqrt(p); Gp = (B + E / sp) / 2; Gm = (B - E / sp) / 2
detB = bareiss_det(B); D = round(math.exp(math.log(detB) / (p - 1)))
a, b, eps = fund_unit(p); t = math.log(eps)
A = (math.exp(t) * Gp + math.exp(-t) * Gm) / math.sqrt(D)      # determinant-1 member at t = log eps
Y, fl, nd = enum_short(A, mult * (1 + 1e-9), maxc=200000, maxnodes=6_000_000_000)
q = np.einsum('ij,jk,ik->i', Y.astype(float), A, Y.astype(float)) if len(Y) else np.zeros(0)
log(f'p={p}: det-1 member at t=log(eps)={t:.4f}: det={np.linalg.det(A):.6f}; vectors (pairs) with norm <= {mult}: {len(Y)} (flag {fl}, nodes {nd}); min norm {q.min() if len(q) else None}; count at min: {int((q < q.min()+1e-6).sum()) if len(q) else 0}')
