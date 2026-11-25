import json, numpy as np, math, sys
R = json.load(open('recheck_main.json')) + json.load(open('recheck_scan.json'))
print('records', len(R), 'with LP', sum(r['has_lp'] for r in R), 'records with error keys', sum(1 for r in R if r['errors']))
mx = lambda xs: max(xs) if xs else None
print('Q: max |det-1| = %.2e, max asym = %.2e, min eig = %.3g' % (mx([abs(r['Q_det'] - 1) for r in R]), mx([r['Q_sym'] for r in R]), min(r['Q_mineig'] for r in R)))
print('lam1sq: max rel diff vs stored = %.2e ; V0 max rel diff = %.2e' % (mx([abs(r['lam1sq'] / r['lam1sq_stored'] - 1) for r in R]), mx([abs(r['V0'] / r['V0_stored'] - 1) for r in R])))
for tag in ('LP', 'K', 'Kd', 'extLP', 'extK'):
    E = [(r, r[tag]) for r in R if tag in r]
    if not E: continue
    print(f'--- {tag}: {len(E)} forms')
    print('   max |min-1| = %.2e' % mx([abs(e['min'] - 1) for r, e in E]))
    print('   pairs mismatches:', sum(e['pairs'] != e['pairs_stored'] for r, e in E), ' V max rel diff = %.2e' % mx([abs(e['V'] / e['V_stored'] - 1) for r, e in E]))
    print('   max V/Vbest = %.9f' % mx([e['V_over_Vbest'] for r, e in E]), ' #forms with V/Vbest > 1-1e-6:', sum(e['V_over_Vbest'] > 1 - 1e-6 for r, e in E))
    print('   non-perfect (rank<N):', sum(e['rank'] != r['n'] * (r['n'] + 1) // 2 for r, e in E), ' max asym %.1e' % mx([e['sym'] for r, e in E]))
    if tag in ('LP', 'K', 'Kd'):
        print('   v max rel diff = %.2e, B0sq max abs diff = %.2e, V via det(Q) rel diff %.2e' % (mx([abs(e['v'] / e['v_stored'] - 1) for r, e in E]), mx([abs(e['B0sq'] - e['B0sq_stored']) for r, e in E]), mx([abs(e['V_viaQ'] / e['V'] - 1) for r, e in E])))
    if tag == 'LP':
        print('   DUAL CERTIFICATE: max rel residual = %.2e ; max |sum lam - Tr(Q^-1 A)|/obj = %.2e ; max |obj - stored lp_obj|/obj = %.2e' % (
            mx([e['dual_resid'] for r, e in E]), mx([abs(e['dual_sum'] - e['primal_obj']) / e['primal_obj'] for r, e in E]), mx([abs(e['primal_obj'] - e['lp_obj_stored']) / e['primal_obj'] for r, e in E])))
        bad = [(r['file'], r['idx'], e['dual_resid']) for r, e in E if e['dual_resid'] > 1e-8]
        print('   records with dual residual > 1e-8:', bad[:10], len(bad))
        print('   eig_ratio max rel diff = %.2e ; rank_stored != N: %d' % (mx([abs(e['eig_ratio'] / e['eig_ratio_stored'] - 1) for r, e in E]), sum(e['rank_stored'] != r['n'] * (r['n'] + 1) // 2 for r, e in E)))
        print('   LP rounds: min %d max %d ; fallback>0: %d ; enlarge>0: %d ; ncons min %d max %d' % (min(e['rounds'] for r, e in E), max(e['rounds'] for r, e in E), sum(e['fallback'] > 0 for r, e in E), sum(e['enlarge'] > 0 for r, e in E), min(e['ncons'] for r, e in E), max(e['ncons'] for r, e in E)))
        print('   dual npos: #(npos == N) = %d ; #(npos < N) = %d ; #(npos > N) = %d ; min positive multiplier overall %.2e' % (
            sum(e['dual_npos'] == r['n'] * (r['n'] + 1) // 2 for r, e in E), sum(e['dual_npos'] < r['n'] * (r['n'] + 1) // 2 for r, e in E), sum(e['dual_npos'] > r['n'] * (r['n'] + 1) // 2 for r, e in E), min(e['dual_minpos'] for r, e in E)))
    if tag in ('extLP', 'extK'):
        eu = [(r['file'], r['idx'], e['eutaxy'], e['steps']) for r, e in E]
        print('   eutaxy margin: min %.3e ; #(<=1e-9): %d ; #steps>=500: %d' % (min(x[2] for x in eu), sum(x[2] <= 1e-9 for x in eu), sum(x[3] >= 500 for x in eu)))
        print('   non-eutactic:', [x for x in eu if not x[2] > 1e-9][:10])
