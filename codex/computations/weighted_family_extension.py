"""Three-dimensional positive weights; exact lattice search without zeta-based selection."""
from fractions import Fraction as F
from math import comb,gcd,lcm,factorial
from pathlib import Path
import json
import mpmath as mp

def triples(n):
    ans=[[F(0),F(0),F(0)] for _ in range(3)]
    unwanted=[{j:F(0) for j in [1,2,4,6]} for _ in range(3)]
    for k in range(n+1):
        num=[-k-j for j in range(1,n+1)]+[-k+n+j for j in range(1,n+1)]
        den=[j-k for j in range(n+1) if j!=k]
        ell=[F(0)]+[F((-1)**(r+1),r)*(2*sum((F(1,z**r) for z in num),F(0))-6*sum((F(1,z**r) for z in den),F(0))) for r in range(1,6)]
        c=[F((comb(n,k)**3*comb(n+k,k)*comb(2*n-k,n-k))**2)]
        for j in range(1,6):c.append(sum((r*ell[r]*c[j-r] for r in range(1,j+1)),F(0))/j)
        v=F(n,2)-k
        hs=[F(0)]+[sum((F(1,l**j) for l in range(1,k+1)),F(0)) for j in range(1,7)]
        for idx,m in enumerate([5,3,1]):
            cs={6-j:sum((2*comb(m,i)*v**(m-i)*c[j-i] for i in range(min(m,j)+1)),F(0)) for j in range(6)}
            for j in [1,2,4,6]: unwanted[idx][j]+=cs[j]
            ans[idx][0]+=cs[5];ans[idx][1]+=cs[3]
            ans[idx][2]+=sum((cs[j]*hs[j] for j in range(1,7)),F(0))
    assert all(value==0 for d in unwanted for value in d.values())
    return ans

def bezout(a,b):
    if b==0:return (1 if a>0 else -1,0,abs(a))
    u,v,g=bezout(b,a%b)
    return v,u-(a//b)*v,g

def dot(u,v):return sum(a*b for a,b in zip(u,v))
def gauss(u,v):
    while True:
        if dot(v,v)<dot(u,u):u,v=v,u
        q=round(F(dot(u,v),dot(u,u)))
        if q==0:return u,v
        v=tuple(b-q*a for a,b in zip(u,v))

def kernel(row):
    a,b,c=row
    u,v,g=bezout(b,c);h=gcd(a,g)
    x=(g//h,-a//h*u,-a//h*v)
    y=(0,c//g,-b//g)
    assert dot(row,x)==dot(row,y)==0
    return gauss(x,y)

def positive(q,n):
    C,D,E=q;z0=(F(3*n,2)+1)**2
    if C<0:return False
    if C==0:return D>=0 and D*z0+E>=0 and (D or E)
    if -F(D,2*C)<=z0:return C*z0*z0+D*z0+E>=0
    return 4*C*E-D*D>=0

def primitive(Z,R):
    d=lcm(Z.denominator,R.denominator);a=int(d*Z);b=int(d*R)
    g=gcd(abs(a),abs(b));sgn=1 if a>=0 else -1
    return sgn*a//g,sgn*b//g

def run(n):
    ts=triples(n)
    d=lcm(*(r[1].denominator for r in ts));row=[int(d*r[1]) for r in ts]
    g=gcd(*row);row=[a//g for a in row]
    u,v=kernel(row)
    seen=set();options=[]
    for a in range(-4,5):
        for b in range(-4,5):
            if not(a or b):continue
            q=tuple(a*x+b*y for x,y in zip(u,v));g=gcd(*q)
            q=tuple(x//g for x in q)
            first=next(x for x in q if x)
            if first<0:q=tuple(-x for x in q)
            if q in seen or not positive(q,n):continue
            seen.add(q)
            Z=sum((F(q[i])*ts[i][0] for i in range(3)),F(0));R=sum((F(q[i])*ts[i][2] for i in range(3)),F(0))
            assert sum((F(q[i])*ts[i][1] for i in range(3)),F(0))==0
            if not Z:continue
            U,V=primitive(Z,R)
            options.append((q,U,V))
    assert options
    chosen={
        'small_weight':min(options,key=lambda o:(dot(o[0],o[0]),o[0])),
        'small_pair':min(options,key=lambda o:(max(abs(o[1]),abs(o[2])),abs(o[1]),o[0]))}
    results={}
    for name,(q,U,V) in chosen.items():
        with mp.workdps(max(len(str(U)),len(str(V)))+100):
            err=mp.mpf(U)*mp.zeta(5)-V
            results[name]=dict(weight=list(map(str,q)),primitive_U=str(U),primitive_V=str(V),log_error_over_n=mp.nstr(mp.log(abs(err))/n,30),error_sign=int(mp.sign(err)),log_U_over_n=mp.nstr(mp.log(abs(U))/n,30))
    nodes=[(F(2*t+n,2))**2 for t in range(n+1,n+4)]
    ratios=[]
    for i,t in enumerate(range(n+1,n+4)):
        others=[nodes[j] for j in range(3) if j!=i]
        denom=(nodes[i]-others[0])*(nodes[i]-others[1])
        L=[1/denom,-sum(others)/denom,others[0]*others[1]/denom]
        ZL=sum((L[j]*ts[j][0] for j in range(3)),F(0))
        H=F(factorial(n)**2)
        for j in range(1,n+1):H*=((t-j)*(t+n+j))**2
        for j in range(n+1):H/=F((t+j)**6)
        ratios.append(abs(ZL)/((2*t+n)*H))
    gap=1/max(ratios)
    det=ts[0][0]*(ts[1][1]*ts[2][2]-ts[1][2]*ts[2][1])-ts[0][1]*(ts[1][0]*ts[2][2]-ts[1][2]*ts[2][0])+ts[0][2]*(ts[1][0]*ts[2][1]-ts[1][1]*ts[2][0])
    assert det!=0
    return dict(n=n,fixed_n_rational_positivity_gap=str(gap),coefficient_map_determinant=str(det),zeta3_row=list(map(str,row)),reduced_kernel_basis=[list(map(str,u)),list(map(str,v))],positive_candidates=len(options),results=results,coefficient_triples=[[str(x) for x in row] for row in ts])

if __name__=='__main__':
    rows=[run(n) for n in range(1,13)]
    S=sum((F(1,k**5) for k in range(1,21)),F(0))
    lo=12544*(S+F(1,4*21**4))-13007
    hi=12544*(S+F(1,4*20**4))-13007
    assert 0<lo<hi<1
    assert sum((F(1,k**5) for k in range(1,5)),F(0))>F(23063,22272)
    assert sum((F(1,k**5) for k in range(1,4)),F(0))+F(1,324)<F(557,536)
    print('PASS exact n2 small-form enclosure and n1 positivity gap')
    Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')
    for row in rows:
        print(row['n'],row['positive_candidates'],[(k,v['log_error_over_n'],v['weight']) for k,v in row['results'].items()])
