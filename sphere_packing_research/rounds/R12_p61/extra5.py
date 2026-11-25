#!/usr/bin/env python3
# complete the remaining t-points (66..74) of the exhaustive Gamma search on Fix(C_3) = split eigen-orders {4,20} (rank 40), p=61
import sys, os, json
from iso import *
p = 61; results = []
Wl = setup(p); rep = Rep(Wl); g0 = Wl.g
split = cyclic_components(rep, (g0, 0, 0, pow(g0, p - 2, p)), 'split torus C_60', results, sums=True)
Et, N, r = split[(4, 20)]
C, exact = sat_basis(Wl, Et, N, r['dimU'], int(r['glue_index']))
sr = Wl.tsearch(C, label='p61 Fix(C_3) rank 40, t-points 64..74', tlim=170, maxnodes=4 * 10**9, trange=(64, 74))
json.dump(sr, open(os.path.join(HERE, 'search_p61_fixC3_tail.json'), 'w'), indent=1, default=str)
log('DONE extra5')
