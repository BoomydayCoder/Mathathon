import json, glob, sys
rows = {}
for f in sorted(glob.glob('log_scan*.jsonl')):
    for line in open(f):
        line = line.strip()
        if not line.startswith('{'): continue
        try: d = json.loads(line)
        except Exception: continue
        key = (d['N'], d['p'])
        if key in rows and rows[key].get('n_roots') is not None and d.get('n_roots') is None: continue
        rows[key] = d
from sympy import isprime
out = ['| N (rank) | p | (p+1)/2 prime | k | roots (pairs by abs s) | norm-4 pairs | log eps | t-steps | fail | hits | min N(nu\') on circle [norm, s] | threshold | verdict |', '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for (N, p) in sorted(rows):
    d = rows[(N, p)]
    b = d.get('best')
    if d['fail']: verdict = 'INCOMPLETE'
    elif b is None: verdict = 'no vector below bound (?)'
    elif b[0] < d['thr']: verdict = '**VIOLATION**'
    elif b[0] == d['thr']: verdict = 'borderline only'
    else: verdict = 'ok'
    out.append('| %d (%d) | %d | %s | %d | %s | %s | %.2f | %d | %d | %d | %s | %d | %s |' % (N, 2 * N, p, 'yes' if isprime((p + 1) // 2) else 'no', d['k'],
        ('%s %s' % (d['n_roots'], d['roots_pairs_by_abs_s'])) if d.get('n_roots') is not None else 'n/c', d.get('norm4_pairs', 'n/c'), d['log_eps'], d['nt'], d['fail'], d['hits'],
        ('%d [%d, %d]' % (b[0], b[1], b[2])) if b else '-', d['thr'], verdict))
print('\n'.join(out))
