#!/usr/bin/env python3
# joint isotypic components for chains of subgroups H1 < H2 (commuting projectors): C_6, C_10, C_4, SL_2(3) inside SL_2(5); C_6, C_4 inside SL_2(3)
import sys, os, json
from iso import *
p = int(sys.argv[1])
results = []
Wl = setup(p); rep = Rep(Wl); rep.test_hom(); assert check_modular(Wl)
n = Wl.n; mI = (p - 1, 0, 0, p - 1)
def find2T(els):
    o6 = [M for M in els if morder(M, p) == 6]
    for a in o6:
        for b in o6:
            ab = mmul(a, b, p)
            if (ab[0] + ab[3]) % p == 0 and len(closure((a, b), p)) == 24: return (a, b)
def joint(A, B, name):
    for la, (EA, NA, ra) in A.items():
        for lb, (EB, NB, rb) in B.items():
            Et = EA @ EB
            assert np.array_equal(Et, EB @ EA), 'projectors do not commute'
            N = NA * NB
            g0 = int(np.gcd.reduce(Et.ravel())); g0 = gcd(g0, N); Et = Et // g0; N = N // g0
            if int(np.trace(Et)) == 0: continue
            comp_stats(Wl, Et, N, '%s: [%s] & [%s]' % (name, la, lb), results)
def relab(d): return {str(k): v for k, v in d.items()}
if p % 10 in (1, 9):
    gens = find_gens(p, '2I')[0]; els = closure(gens, p)
    o2I = group_components(rep, gens, 'SL_2(5)', CH_2I, results, expect_order=120)
    g2T = find2T(els)
    o2T = group_components(rep, g2T, 'SL_2(3)<SL_2(5)', CH_2T, results, expect_order=24, relcheck=relcheck_poly(3))
    c6 = relab(cyclic_components(rep, g2T[0], 'C_6<SL_2(3)<SL_2(5)', results, sums=False))
    c10 = relab(cyclic_components(rep, gens[1], 'C_10<SL_2(5)', results, sums=False))
    x4 = [M for M in els if morder(M, p) == 4][0]
    c4 = relab(cyclic_components(rep, x4, 'C_4<SL_2(5)', results, sums=False))
    joint(o2I, o2T, 'SL_2(5)>SL_2(3)'); joint(o2I, c6, 'SL_2(5)>C_6'); joint(o2I, c10, 'SL_2(5)>C_10'); joint(o2T, c6, 'SL_2(3)>C_6')
    # triple chain
    for la, (EA, NA, ra) in o2I.items():
        for lb, (EB, NB, rb) in o2T.items():
            for lc, (EC, NC, rc) in c6.items():
                Et = EA @ EB @ EC; N = NA * NB * NC
                g0 = gcd(int(np.gcd.reduce(Et.ravel())), N); Et = Et // g0; N = N // g0
                if int(np.trace(Et)) == 0: continue
                comp_stats(Wl, Et, N, 'SL_2(5)>SL_2(3)>C_6: [%s] & [%s] & [%s]' % (la, lb, lc), results)
else:
    gens = find_gens(p, '2T')[0]; els = closure(gens, p)
    o2T = group_components(rep, gens, 'SL_2(3)', CH_2T, results, expect_order=24, relcheck=relcheck_poly(3))
    c6 = relab(cyclic_components(rep, gens[0], 'C_6<SL_2(3)', results, sums=False))
    joint(o2T, c6, 'SL_2(3)>C_6')
json.dump(results, open(os.path.join(HERE, 'extra2_p%d.json' % p), 'w'), indent=1, default=str)
log('DONE extra2', p)
