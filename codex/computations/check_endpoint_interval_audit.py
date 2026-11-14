"""Independent exact-fraction tests of the endpoint interval primitives.
The enclosure proof is in endpoint-independent-review.md; these are diagnostics.
"""
import importlib.util
from pathlib import Path
from fractions import Fraction as Q
from random import Random
path=Path(__file__).with_name('certify_equilibrium_endpoints.py')
spec=importlib.util.spec_from_file_location('endpoint_certificate',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
I,J,S=module.I,module.J,module.S
rng=Random(131)
def encloses(interval,values):
    assert Q(interval.lo,S)<=min(values)
    assert Q(interval.hi,S)>=max(values)
for _ in range(1000):
    aa=sorted(Q(rng.randrange(-100,101),rng.randrange(1,98)) for _ in range(2))
    bb=sorted(Q(rng.randrange(-100,101),rng.randrange(1,98)) for _ in range(2))
    a,b=I.range(*aa),I.range(*bb)
    encloses(a+b,[x+y for x in aa for y in bb])
    encloses(a-b,[x-y for x in aa for y in bb])
    encloses(a*b,[x*y for x in aa for y in bb])
    if bb[0]>0 or bb[1]<0:
        encloses(a/b,[x/y for x in aa for y in bb])
    if aa[0]>=0:
        r=a.sqrt()
        assert Q(r.lo,S)**2<=aa[0] and Q(r.hi,S)**2>=aa[1]
for left,right in [(Q(0),Q(1,100)),(Q(1,5),Q(1,4)),(Q(2),Q(21,10))]:
    x=J(I.range(left,right),True)
    reciprocal=1/(1+x*x)
    for point in (left,(left+right)/2,right):
        exact=(5*point**4-10*point**2+1)/(1+point**2)**5
        encloses(reciprocal.c[4],[exact])
    linear=((x+2)*(x+2)).sqrt()
    encloses(linear.c[0],[left+2,right+2])
    encloses(linear.c[1],[Q(1)])
    for coefficient in linear.c[2:]:encloses(coefficient,[Q(0)])
print('PASS 1000 exact-fraction interval-operation checks and representative whole-cell derivative jets')
