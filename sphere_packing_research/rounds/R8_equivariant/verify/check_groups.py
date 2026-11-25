"""(1) Structural checks: group order, dim Q (character formula + Reynolds rank) vs theory and vs eq.py's D,
lattice = O^k (covolume of O, index), invariance of eq.py's Q basis under my group."""
import numpy as np, math, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import eq as their

phi = {5: 4, 7: 6, 8: 4, 12: 4, 13: 12}


def theory_D(name, k):
    if name in ('eis', 'gauss'):
        return k * k
    if name == 'hur':
        return k * (2 * k - 1)
    if name == 'ico':
        return 2 * k * (2 * k - 1)
    return k * k * phi[int(name[4:])] // 2


cases = [('eis', k) for k in (1, 2, 3, 4, 5, 6, 8)] + [('gauss', k) for k in range(1, 7)] + [('hur', k) for k in range(1, 7)] + \
        [('zeta5', 2), ('zeta5', 3), ('zeta7', 2), ('zeta8', 2), ('zeta8', 3), ('zeta8', 4), ('zeta12', 2), ('zeta12', 3),
         ('zeta12', 4), ('zeta13', 2), ('ico', 1), ('ico', 2), ('ico', 3)]
out = []
cache = {}
for name, k in cases:
    if name not in cache:
        els, w, n0, extra = vlib.base_group(name)
        cache[name] = (els, w, n0, extra, vlib.is_group(els))
        if name == 'ico':
            idx = vlib.lattice_index(extra['unit_coords'])
            print('icosian: |det Bm| =', extra['covol'], '(maximal order: 25/16 = 1.5625); index of Z-span of my 120 units in their lattice =', idx)
            Bm = extra['Bm']
            # second half conjugate of first? test: each basis column c: first4 = p+q*sqrt5, second4 = p-q*sqrt5 with p,q in (1/4)Z
            s5 = math.sqrt(5)
            P = (Bm[:4] + Bm[4:]) / 2
            Qm = (Bm[:4] - Bm[4:]) / (2 * s5)
            print('   basis columns of form (x, conj x) with x in (1/4)Z[sqrt5]^4:', np.abs(4 * P - np.rint(4 * P)).max() < 1e-9 and np.abs(4 * Qm - np.rint(4 * Qm)).max() < 1e-9)
        if name == 'hur':
            print('hurwitz: covol(O) =', extra['covol'], '(expected 1/2)')
    els, w, n0, extra, isg = cache[name]
    n = n0 * k
    Dchar = vlib.dimQ_character(els, k)
    QB = vlib.Qbasis_reynolds(els, k) if n <= 24 else None
    Drey = len(QB)
    E = their.Eq(name, k)
    # is their Q basis invariant under ALL my group elements, and same span as mine?
    gs = [np.kron(np.eye(k), g.astype(float)) for g in els]
    inv = max(np.abs(g.T @ B @ g - B).max() for g in gs for B in E.QB)
    # span equality: project their basis on mine
    Mm = QB.reshape(Drey, -1)
    Tm = E.QBf
    res = np.abs(Tm - (Tm @ Mm.T) @ Mm).max()
    row = dict(name=name, k=k, n=n, w_mine=len(els), w_theirs=E.w, group_closed=bool(isg), D_theory=theory_D(name, k), D_char=str(Dchar),
               D_reynolds=Drey, D_theirs=E.D, their_group_order=len(E.group), theirQ_invariance_err=float(inv), span_residual=float(res))
    out.append(row)
    ok = (len(els) == E.w == len(E.group) and Dchar == theory_D(name, k) == Drey == E.D and inv < 1e-9 and res < 1e-8 and isg)
    print(('OK  ' if ok else 'BAD ') + json.dumps(row), flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'check_groups.json'), 'w'), indent=1)
