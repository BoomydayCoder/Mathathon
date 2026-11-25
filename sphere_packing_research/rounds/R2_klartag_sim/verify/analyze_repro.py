#!/usr/bin/env python3
"""Summarize the reproduction + independent re-verification (repro_*.jsonl)."""
import json, glob, math
def kappa(n): return math.pi ** (n / 2) / math.gamma(n / 2 + 1)
best = {2: (2, 3), 3: (2, 4), 4: (2, 4), 5: (2, 4), 6: (2, 3), 7: (2, 2), 8: (2, 1), 9: (4, 512), 10: (4, 768)}
def Vb(n): m, dd = best[n]; return kappa(n) * m ** (n / 2) / math.sqrt(dd)
tot = 0; bad = []
for f in sorted(glob.glob('repro_*.jsonl')):
    R = [json.loads(l) for l in open(f) if l.strip()]
    err = [r for r in R if 'error' in r]; R = [r for r in R if 'error' not in r]
    if not R: print(f, 'no runs', err); continue
    n = R[0]['n']; d = n * (n + 1) // 2; vb = Vb(n)
    nm = sum(r['lam1_match'] and r['V_reldiff'] < 1e-9 and r['nevents_match'] for r in R)
    minq = min(r['my_minq'] for r in R); maxq = max(r['my_minq'] for r in R)
    pairs_ok = sum(r['my_nmin'] == r['saved_npairs'] for r in R)
    perfect = sum(r['my_rank'] == d for r in R)
    Vok = sum(abs(r['my_V'] - r['Vsaved']) < 1e-9 * r['Vsaved'] for r in R)
    over = [r['i'] for r in R if r['my_Vhonest'] > vb * (1 + 1e-7)]
    boxn = sum(r['box_minq'] is not None for r in R)
    boxok = sum(r['box_minq'] is not None and abs(r['box_minq'] - r['my_minq']) < 1e-9 and r['box_nmin'] == r['my_nmin'] for r in R)
    vmax = max(r['my_V'] for r in R); vmaxsaved = max(r['Vsaved'] for r in R)
    tot += len(R)
    print(f"{f}: n={n} reproduced {len(R)} runs; bit-identical to saved (lam1, V rel<1e-9, #events): {nm}; own min A x.x in [{minq:.10f}, {maxq:.10f}]; "
          f"own #min pairs == saved: {pairs_ok}; perfect (rank d): {perfect}; own V == saved V: {Vok}; max own V {vmax:.4f} (saved max in sample {vmaxsaved:.4f}, best known {vb:.4f}); exceed best: {over}; "
          f"box-enum available {boxn}, agrees {boxok}; errors: {err}")
    for r in R:
        if not (r['lam1_match'] and r['V_reldiff'] < 1e-9 and r['my_nmin'] == r['saved_npairs'] and r['my_rank'] == d and r['my_minq'] > 1 - 1e-7): bad.append((f, r['i'], r['V_reldiff'], r['my_minq'], r['my_nmin'], r['saved_npairs'], r['my_rank']))
print('total verified forms', tot); print('non-conforming records:', bad)
