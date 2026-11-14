"""Audit exact first-integral identity and finite digit-polynomial shapes."""
from math import comb
import sympy as s


def mul(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%p
    while len(out)>1 and out[-1]==0:out.pop()
    return out


def div_delta(a,p):
    out=a[:]
    q=[0]*max(1,len(a)-2)
    for i in range(len(a)-1,1,-1):
        u=out[i]%p
        q[i-2]=u
        out[i]=(out[i]-u)%p
        out[i-1]=(out[i-1]+34*u)%p
        out[i-2]=(out[i-2]-u)%p
    if any(out[:2]):return None
    while len(q)>1 and q[-1]==0:q.pop()
    return q


def square_root(a,p):
    if (len(a)-1)%2:return None
    d=(len(a)-1)//2
    assert a[0]==1
    q=[1]
    for n in range(1,d+1):
        q.append((a[n]-sum(q[j]*q[n-j] for j in range(1,n)))*pow(2,-1,p)%p)
    return q if mul(q,q,p)==a else None


def main():
    t=s.symbols('t')
    y=s.Function('y')(t)
    delta=1-34*t+t*t
    a=t*t*delta;b=2*t*(1-51*t+2*t*t);c=t*(t-10)
    I=a*(2*y*s.diff(y,t,2)-s.diff(y,t)**2)+b*y*s.diff(y,t)+c*y*y
    Ldiv=t*t*delta*s.diff(y,t,3)+3*t*(1-51*t+2*t*t)*s.diff(y,t,2)+(1-112*t+7*t*t)*s.diff(y,t)+(t-5)*y
    assert s.expand(s.diff(I,t)-2*y*Ldiv)==0
    assert s.expand(b-s.diff(a,t))==0
    assert s.expand(b-t*t*s.diff(delta,t)-2*t*delta)==0
    count=0; eps0=[];eps1=[]
    for p in list(s.primerange(5,102)):
        A=[sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))%p for n in range(p)]
        assert A[0]==A[-1]==1 and A==A[::-1]
        roots=[]
        for epsilon in (0,1):
            R=A if epsilon==0 else div_delta(A,p)
            if R is None:continue
            Q=square_root(R,p)
            if Q is not None:roots.append((epsilon,Q))
        assert len(roots)==1
        epsilon,Q=roots[0]
        assert len(Q)-1==(p-1-2*epsilon)//2
        poly=sum(v*t**j for j,v in enumerate(Q))
        assert s.degree(s.gcd(s.Poly(poly,t,modulus=p),s.Poly(s.diff(poly,t),t,modulus=p)))==0
        assert s.degree(s.gcd(s.Poly(poly,t,modulus=p),s.Poly(delta,t,modulus=p)))==0
        # Experimental check only: the symbol rule is not used in the proof.
        assert epsilon==(1-int(s.legendre_symbol(-6,p)))//2
        (eps1 if epsilon else eps0).append(p)
        count+=1
    print('Exact first-integral derivative identity passed.')
    print(f'Square/Delta-square, squarefreeness and coprimality checked at {count} primes 5..101.')
    print('epsilon=0:',eps0)
    print('epsilon=1:',eps1)
    print('Legendre-symbol rule agrees in this finite range; not proved here.')

if __name__=='__main__':main()
