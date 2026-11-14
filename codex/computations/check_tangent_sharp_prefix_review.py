"""Independent sharp-prefix audit: ODE convolution, digit carries, finite-field delta.
No source checker or project helper is imported; standard library only.
"""
from math import comb, gcd


def factor(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out


def vp(n,p):
    assert n
    v=0
    while n%p==0:n//=p;v+=1
    return v


def carry_count(x,y,p):
    carry=0;count=0
    while x or y or carry:
        carry=(x%p+y%p+carry)//p
        count+=carry;x//=p;y//=p
    return count


def tangent_from_ode(N):
    # y'=1+y^2 implies y^(n+1)(0)=1_(n=0)+sum binom(n,k)y^k(0)y^(n-k)(0).
    values=[0]
    for n in range(N):
        values.append(int(n==0)+sum(comb(n,k)*values[k]*values[n-k] for k in range(n+1)))
    return values


def delta(poly,p):
    out={}
    for degree,c in poly.items():
        if degree:
            for j in (degree-1,degree+1):out[j]=(out.get(j,0)+degree*c)%p
    return {j:c for j,c in out.items() if c}


def main():
    N=401;T=tangent_from_ode(N)
    assert [T[r] for r in (1,3,5,7,9)]==[1,2,16,272,7936]
    assert all(T[r]==0 for r in range(0,N+1,2))
    counts={'prefixes':0,'exact_binomial_carries':0,'witnesses':0,'operator_monomials':0,'operator_prime_powers':0}
    for s in range(3,N+1,2):
        fs=factor(s);Q=max(p**h for p,h in fs.items());E=(Q-1)//2
        running=0;first=None
        for e in range(E+1):
            A=T[2*e+1]*comb(s+2*e,s-1)
            running=gcd(running,A)
            if running==1 and first is None:first=e
        assert first==E
        assert (Q==s)==(len(fs)==1)
        for p,h in fs.items():
            assert p!=2
            P=p**h;u=s//P;critical=(P-1)//2
            for e in range(critical+1):
                r=2*e+1;b=comb(s+r-1,s-1)
                expected=h-vp(r,p)
                assert vp(b,p)==expected==carry_count(s-1,r,p)
                if e<critical:assert b%p==0
                counts['exact_binomial_carries']+=1
            b=comb(s+P-1,P)
            scalar=(-1)**(h*(p-1)//2)%p
            assert b%p==u%p and T[P]%p==scalar
            assert (b*T[P])%p!=0
            counts['witnesses']+=1
        counts['prefixes']+=1
    for p in (3,5,7,11,13,17,19):
        alpha=(-1)**((p-1)//2)%p
        for degree in range(2*p+1):
            start={degree:1};result=start
            for _ in range(p):result=delta(result,p)
            expected={j:alpha*c%p for j,c in delta(start,p).items()}
            expected={j:c for j,c in expected.items() if c}
            assert result==expected
            counts['operator_monomials']+=1
        h=1;P=p
        while P<=N:
            result={1:1}
            for _ in range(P):result=delta(result,p)
            scalar=pow(alpha,h,p)
            assert result=={0:scalar,2:scalar}
            assert result[0]==T[P]%p
            counts['operator_prime_powers']+=1
            h+=1;P*=p
    # Non-prime-power and prime-power sharpness examples.
    assert (max(p**h for p,h in factor(15).items())-1)//2==2
    assert (max(p**h for p,h in factor(9).items())-1)//2==4
    print('PASS independent sharp-prefix audit:',counts)

if __name__=='__main__':main()
