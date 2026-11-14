"""Exact certificates for sharp cone slope gaps n=1..5."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
from weighted_family_extension import triples

def coefficient(q,ts,index):return sum((q[j]*ts[j][index] for j in range(3)),F(0))
def run(n,L=None):
    ts=triples(n)
    if L is None:L=F(9*n*n,4)
    ell=[ts[0][1]-2*L*ts[1][1]+L*L*ts[2][1],ts[1][1]-L*ts[2][1],ts[2][1]]
    assert ell[1] and ell[0]*ell[2]>0 and ell[0]*ell[1]<0
    u=-ell[0]/ell[1];v=-ell[2]/ell[1]
    assert u>0 and v>0
    qa=[F(1),u-2*L,L*L-u*L]
    qc=[F(0),v,1-v*L]
    for q in [qa,qc]:assert coefficient(q,ts,1)==0
    Z1=coefficient(qa,ts,0);R1=coefficient(qa,ts,2)
    Z2=coefficient(qc,ts,0);R2=coefficient(qc,ts,2)
    assert Z1>0 and Z2<0
    lower=R1/Z1;upper=R2/Z2
    assert lower<upper
    tt=n+1 if L<(F(3*n,2)+1)**2 else n+2;zz=(F(2*tt+n,2))**2;yy=zz-L
    H=F(factorial(n)**2)
    for j in range(1,n+1):H*=((tt-j)*(tt+n+j))**2
    for j in range(n+1):H/=F((tt+j)**6)
    w=(2*tt+n)*H
    lower_gap=w*(yy*yy+u*yy)/Z1
    upper_gap=w*(1+v*yy)/(-Z2)
    assert lower_gap>0 and upper_gap>0
    return dict(n=n,L=str(L),ell=list(map(str,ell)),u=str(u),v=str(v),qa=list(map(str,qa)),qc=list(map(str,qc)),Z1=str(Z1),R1=str(R1),Z2=str(Z2),R2=str(R2),lower=str(lower),upper=str(upper),width=str(upper-lower),first_term_lower_gap=str(lower_gap),first_term_upper_gap=str(upper_gap),primitive_error_lower_bound=str(min(lower_gap,upper_gap)))
if __name__=='__main__':
    rows=[run(n) for n in range(1,6)]
    for row in rows:
        n=row['n']
        row['discrete_domain']=run(n,(F(3*n,2)+1)**2)
    assert rows[0]['u']=='88/37' and rows[0]['v']=='78/185'
    assert rows[0]['lower']=='659/640' and rows[0]['upper']=='557/536'
    Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('PASS exact cone sign conditions, extreme ray coefficient signs, sharp rational slope intervals, and positive rational gap certificates n=1..5')
    for row in rows:print(row['n'],row['lower'],row['upper'],float(F(row['width'])))
