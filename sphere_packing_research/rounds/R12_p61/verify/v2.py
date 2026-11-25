import sys
src = open('v.py').read().split('# ---- split torus fixed lattices')[0]
sys.argv = [sys.argv[0], sys.argv[1], '--nodet']
exec(src)
P('O_F-stability: (I+M)/2 integral:', bool(((M + I) % 2 == 0).all()))
chi3 = {1: 2, 2: -2, 4: 0, 3: -1, 6: 1}
if p == 37:
    ents = [e for e in json.load(open(R12 + '/search_p37.json')) if 'Hurwitz' in e['label']]
    w = np.array(ents[0]['w'], dtype=np.int64)
    for gi, (s, t) in enumerate(iso.find_gens(p, '2T')):
        els = gen_group(s, t, 24)
        Et = sum(chi3[order2(x)] * Px for x, Px in els.items())
        assert (Et @ Et == 12 * Et).all()
        P('tester SL_2(3) gens #%d %s: dim %d, saved Hurwitz violator in component: %s' % (gi, (s, t), np.trace(Et) // 12, not (w @ (12 * I - Et)).any()))
        if gi == 0:
            res, C = certify(Et, 12, 'tester SL_2(3) #0 Hurwitz')
    b = int(w.astype(object) @ Bo @ w.astype(object)); e = int(w.astype(object) @ Eo @ w.astype(object))
    P('vector b,e =', b, e, 'Gamma/theta =', Fraction(b * b * p - e * e, 4 * p * D) / Fraction((p - 1) ** 2, 4 * p))
if p == 61:
    chi_ico = {1: 4, 2: -4, 4: 0, 3: -2, 6: 2, 5: -1, 10: 1}; chi4 = {1: 4, 2: -4, 4: 0, 3: 1, 6: -1, 5: -1, 10: 1}
    chi6 = {1: 6, 2: -6, 4: 0, 3: 0, 6: 0, 5: 1, 10: -1}
    s, t = poly_group(5, 1); els = gen_group(s, t, 120)
    Ei = sum(chi_ico[order2(x)] * Px for x, Px in els.items()); E4 = sum(chi4[order2(x)] * Px for x, Px in els.items()); E6 = sum(chi6[order2(x)] * Px for x, Px in els.items())
    # e_ico = Ei/60, e_4 = E4/30, e_6 = E6/20 ; sum must be identity
    P('e_ico + e_4 + e_6 == 1 exactly:', bool((Ei + 2 * E4 + 3 * E6 == 60 * I).all()), ' dims', np.trace(Ei) // 60, np.trace(E4) // 30, np.trace(E6) // 20)
    res, C = certify(Ei + 2 * E4, 60, 'SL_2(5) ico+deg4 (complement of deg-6)')
    k = 72; g6 = 2 ** 8 * 5 ** 40
    P('  => deg-6 comp: same glue; claimed g=2^8 5^40; R(72) = %.4g ; R(48) with that g = %.4g' % (math.exp(72 * math.log(120 / 72) - math.log(g6)), math.exp(48 * math.log(120 / 48) - math.log(g6))))
    # SL_2(3) inside this SL_2(5): Hurwitz comp + joint
    o6 = [x for x in els if order2(x) == 6]; o4 = [x for x in els if order2(x) == 4]
    found = None
    for u in o6[:3]:
        for v in o4:
            sub = {ID}; fr = [ID]
            while fr and len(sub) <= 24:
                new = []
                for x in fr:
                    for gg in (u, v):
                        y = mm(x, gg)
                        if y not in sub: sub.add(y); new.append(y)
                fr = new
            if len(sub) == 24: found = sub; break
        if found: break
    Eh = sum(chi3[order2(x)] * els[x] for x in found)
    res, C = certify(Eh, 12, 'SL_2(3) inside SL_2(5): Hurwitz comp')
P('done')
