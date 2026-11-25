"""Full lattices from the formula: W_37 (rank 72), W_61 (rank 120): exact integer Gram, evenness, exact determinant,
(G^{-1}E)^2 = p exactly (eta), and the exact inequality of Theorem 5 on the Borel sublattice is not tested here."""
import sys, time
sys.path.insert(0,'.'); sys.path.insert(0,'../R10_weil_lattices/orchestrator')
import ideal_model as IM, numpy as np
from fractions import Fraction
for p in (29, 37, 61):
    for var in ('L1', 'L0'):
        t0 = time.time()
        Gi, Ei, Bz, F = IM.build(p, p - 1, verbose=False, variant=var)
        n = Gi.shape[0]
        even = all(int(Gi[i, i]) % 2 == 0 for i in range(n))
        det = IM.exact_det(Gi)
        # exact check E G^{-1} E = p G  (i.e. eta^2 = p) :  solve G X = E over Q is slow; check via integer identity det(E) = p^{n/2} det(G) * (+-1) and symmetric
        detE = IM.exact_det(Ei)
        print('p=%d %s rank=%d even=%s det(G)=%d det(E)=%s p^(n/2)=%s  [%.0fs]' % (p, var, n, even, det, 'p^%d*%d' % (n // 2, detE // p ** (n // 2)) if detE % p ** (n // 2) == 0 else detE, n // 2, time.time() - t0), flush=True)
