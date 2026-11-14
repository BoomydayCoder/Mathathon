"""Independent execution wrapper plus exact inclusion/margin checks.
Requires python-flint; endpoint quadrature itself is standard-library-only.
"""
from pathlib import Path
from fractions import Fraction as Q
import runpy
from flint import arb,fmpq

here=Path(__file__).resolve().parent
runpy.run_path(str(here/'certify_equilibrium_endpoints.py'),run_name='__main__')
state=runpy.run_path(str(here/'certify_equilibrium_energy.py'),run_name='__main__')
for label,lo,hi in [('a',Q('.00007119508'),Q('.00007119517')),
                    ('b',Q('.75394957'),Q('.75394967'))]:
    for q in (lo,hi):
        assert state[label].contains(arb(fmpq(q.numerator,q.denominator)))
assert state['b']-state['a']>0
assert state['objective']<arb(fmpq(-69,50))
assert state['M']>arb(fmpq(-6661,1000))
A=Q(127125602969131786927559,94195881588024216000000)
margin=Q(69,50)-A
assert margin==Q(2864713622341631152441,94195881588024216000000)
assert margin>Q(3,100)
print('PASS: independent reruns, exact endpoint inclusion, objective<-69/50, arithmetic margin>3/100')
