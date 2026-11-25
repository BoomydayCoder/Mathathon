#!/usr/bin/env python3
# extra tests: (a) p=29 cross-validation against builder B's exact R_U>1 values (non-split C_3, C_5 fixed lattices);
# (b) second conjugacy class of SL_2(5), SL_2(3) (conjugation by diag(g,1), outer automorphism); (c) all sums of components
import sys, os, json
from iso import *
p = int(sys.argv[1]); mode = sys.argv[2]
results = []; searches = []
Wl = setup(p); rep = Rep(Wl); rep.test_hom(); assert check_modular(Wl)
n = Wl.n; g0 = Wl.g
def conj(M):
    a, b, c, d = M; gi = pow(g0, p - 2, p)
    return (a, b * g0 % p, c * gi % p, d)     # diag(g,1) M diag(g,1)^-1
def nonsplit():
    rng = np.random.default_rng(5)
    while True:
        a, b, c = [int(v) for v in rng.integers(0, p, 3)]
        if a == 0 or b == 0: continue
        d = (1 + b * c) * pow(a, p - 2, p) % p; M = (a, b, c, d)
        if morder(M, p) == p + 1: return M
if mode == 'ns':
    ns = cyclic_components(rep, nonsplit(), 'non-split torus C_%d' % (p + 1), results, sums=True)
    split = cyclic_components(rep, (g0, 0, 0, pow(g0, p - 2, p)), 'split torus C_%d' % (p - 1), results, sums=True)
    poly = {}
    gl = find_gens(p, '2T'); poly['2T'] = group_components(rep, gl[0], 'SL_2(3)', CH_2T, results, expect_order=24, relcheck=relcheck_poly(3))
    if p % 10 in (1, 9):
        gl = find_gens(p, '2I'); poly['2I'] = group_components(rep, gl[0], 'SL_2(5)', CH_2I, results, expect_order=120, relcheck=relcheck_poly(5))
if mode == 'class2':
    for kind, nm, CH, od, k in (('2T', 'SL_2(3)', CH_2T, 24, 3), ('2I', 'SL_2(5)', CH_2I, 120, 5)):
        if kind == '2I' and p % 10 not in (1, 9): continue
        gl = find_gens(p, kind)
        for i, gens in enumerate(gl[:3]):
            g2 = (conj(gens[0]), conj(gens[1]))
            out = group_components(rep, g2, nm + ' [outer-twisted class, gens %d]' % i, CH, results, expect_order=od, relcheck=relcheck_poly(k))
            # all pairwise sums of components
            labs = list(out)
            if len(labs) == 3 and i == 0:
                for a_, b_ in itertools.combinations(labs, 2):
                    (E1, N1, r1), (E2, N2, r2) = out[a_], out[b_]
                    L = N1 * N2 // gcd(N1, N2)
                    comp_stats(Wl, E1 * (L // N1) + E2 * (L // N2), L, nm + ' twisted: sum [%s]+[%s]' % (a_, b_), results)
            if i == 0:
                for lab, (Et, N, r) in out.items():
                    if r and r['dimU'] <= 40:
                        C, exact = sat_basis(Wl, Et, N, r['dimU'], int(r['glue_index']))
                        sr = Wl.tsearch(C, label='p%d %s' % (p, r['label']), tlim=100, maxnodes=2 * 10**9); sr['det_certified_exact'] = exact; searches.append(sr)
json.dump(dict(results=results, searches=searches), open(os.path.join(HERE, 'extra_p%d_%s.json' % (p, mode)), 'w'), indent=1, default=str)
log('DONE extra', p, mode)
