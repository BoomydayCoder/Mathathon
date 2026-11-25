import sys
src = open('v.py').read().split('# ---- split torus fixed lattices')[0]
sys.argv = [sys.argv[0], '61', '--nodet']
exec(src)
from math import gcd, cos, pi
a = g0; Pd = op((a, 0, 0, pow(a, p - 2, p))); pw = [I]
for k in range(1, 60): pw.append(pw[-1] @ Pd)
assert (pw[-1] @ Pd == I).all()
def ram(d, k): return int(round(sum(cos(2 * pi * j * k / d) for j in range(d) if gcd(j, d) == 1)))
atoms = {d: sum(ram(d, k) * pw[k] for k in range(60)) for d in (4, 12, 20, 60)}
P('atom dims', {d: int(np.trace(Et)) // 60 for d, Et in atoms.items()}, ' sum of atoms == 60 I:', bool((sum(atoms.values()) == 60 * I).all()))
for ds in ((12,), (20,), (12, 20), (60,)):
    try:
        certify(sum(atoms[d] for d in ds), 60, 'split C_60 eigen-orders %s' % (ds,))
    except Exception as ex:
        P('FAILED', ds, repr(ex))
P('done')
