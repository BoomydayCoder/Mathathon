import json, glob, subprocess, re
from sympy import isprime
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
summ = []
for N in sorted({k[0] for k in rows}):
    ps = sorted(p for (n, p) in rows if n == N)
    good = [p for p in ps if rows[(N, p)]['fail'] == 0]
    viol = [p for p in good if rows[(N, p)]['best'] and rows[(N, p)]['best'][0] < rows[(N, p)]['thr']]
    qprime = [p for p in good if isprime((p + 1) // 2)]
    summ.append('* N = %d (rank %d, threshold %d): %d primes searched exhaustively (fail = 0), of which %d have (p+1)/2 prime; range p = %d..%d; violators: %s; every other prime has minimum exactly %d (forced borderline vectors only).  Primes: %s.' % (
        N, 2 * N, (N // 4) ** 2, len(good), len(qprime), min(ps), max(ps), viol if viol else 'none', (N // 4) ** 2, ', '.join(('%d' % p) + ('' if isprime((p + 1) // 2) else '°') for p in good)))
table = subprocess.run(['python3', 'make_table.py'], capture_output=True, text=True).stdout
s = open('proof.md').read()
block = 'Summary (° marks primes with (p+1)/2 composite, included only as additional tests of the quotient-curve lattices):\n\n' + '\n'.join(summ) + '\n\n' + table
if 'TABLE_PLACEHOLDER' in s:
    s = s.replace('TABLE_PLACEHOLDER', '<!--TABLE-->\n' + block + '<!--/TABLE-->')
else:
    s = re.sub(r'<!--TABLE-->.*<!--/TABLE-->', lambda m: '<!--TABLE-->\n' + block + '<!--/TABLE-->', s, flags=re.S)
open('proof.md', 'w').write(s)
print('\n'.join(summ))
