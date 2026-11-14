"""Exact exploration, not proof of eventual signs or growth rates."""
from fractions import Fraction as F
from math import comb, log
import json
from pathlib import Path

def moments(n):
    a3=b3=a5=b5=F(0)
    for k in range(n+1):
        c0=(comb(n,k)**3*comb(n+k,k)*comb(2*n-k,n-k))**2
        den=[j-k for j in range(n+1) if j!=k]
        num=[-k-j for j in range(1,n+1)]+[-k+n+j for j in range(1,n+1)]
        def h(seq,r): return sum((F(1,z**r) for z in seq),F(0))
        l1=2*h(num,1)-6*h(den,1)
        l2=-h(num,2)+3*h(den,2)
        l3=F(2,3)*h(num,3)-2*h(den,3)
        c1=c0*l1;c2=c0*(l2+l1*l1/2);c3=c0*(l3+l1*l2+l1**3/6)
        v=F(n,2)-k
        a3+=2*(c0+3*v*c1+3*v*v*c2+v**3*c3)
        b3+=2*(c2+v*c3)
        a5+=2*(3*v*v*c0+v**3*c1)
        b5+=2*(c0+v*c1)
    return a3,b3,a5,b5

rows=[]
for n in range(1,41):
    a3,b3,a5,b5=moments(n)
    assert b3!=0
    lam=-a3/b3
    assert lam>0
    z5=a5+lam*b5
    assert z5<0
    rows.append(dict(n=n,lambda_exact=str(lam),lambda_over_n_squared=float(lam/n**2),log_lambda_denominator_over_n=log(lam.denominator)/n,zeta5_coefficient_exact=str(z5)))
Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS finite observations n=1..40: b3 nonzero, lambda positive, zeta5 coefficient negative')
for n in [1,2,3,4,5,10,20,40]:
    row=rows[n-1]
    print(n, row['lambda_over_n_squared'],row['log_lambda_denominator_over_n'])
