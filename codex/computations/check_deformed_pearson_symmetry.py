"""Exact checks of Pearson reciprocity and generic associated recurrences.

All arithmetic is rational. Z_sigma is a formal independent symbol,
represented by an ordered pair (its coefficient, the constant term).
"""
from fractions import Fraction as F
from math import comb, factorial
import check_higher_order_independent as A


def compose(p,q):
    out=[F(0)]
    for c in reversed(p):out=A.add(A.mul(out,q),[c])
    return out


def shifted(p,c):return compose(p,[F(c),F(1)])


def divide_linear(p,s):
    p=list(p)
    if len(p)==1:
        assert p[0]==0
        return [F(0)]
    q=[F(0)]*(len(p)-1)
    for j in range(len(p)-1,0,-1):
        q[j-1]=p[j]
        p[j-1]+=s*p[j]
    assert p[0]==0
    return A.trim(q)


def fadd(*terms):
    return tuple(sum(t[j] for t in terms) for j in range(2))


def fscale(c,t):return tuple(c*x for x in t)


parameters=[F(0),F(1,4),F(1,9),F(1,5),F(-1),F(2,3)]
counts=dict(commutator=0,pearson=0,associated=0,grid=0,adjacent=0,diagonal=0,kernel=0)
for sigma in parameters:
    ell=[sum((F(2*comb(d+1,2*h),d+1)*sigma**(h-1)*A.B[d+1-2*h]
              for h in range(1,(d+1)//2+1)),F(0)) for d in range(31)]
    def LL(p):return sum(c*ell[j] for j,c in enumerate(p))
    def theta(p):
        out=[F(0)]*max(1,len(p)-1)
        for d,c in enumerate(p):
            for j in range(d):out[d-1-j]+=c*ell[j]
        return A.trim(out)
    v=[F(0),-sigma,F(0),F(1)]
    vp=shifted(v,1)
    def L(p):
        return A.add(A.mul(vp,A.add(shifted(p,1),A.scale(p,-1))),
                     A.mul(v,A.add(shifted(p,-1),A.scale(p,-1))))
    def vv(j):return F(j)*(j*j-sigma)
    assert ell[:2]==[0,1]
    for d in range(19):
        p=[F(0)]*d+[F(1)]
        assert LL(compose(p,[F(-1),F(-1)]))==-LL(p)
        lhs=A.add(L(theta(p)),A.scale(theta(L(p)),-1))
        rhs=A.add(A.scale(A.add(shifted(p,-1),A.scale(shifted(p,1),-1)),2),[LL(p)])
        assert lhs==rhs
        counts["commutator"]+=1
        if d<=14:
            vh=A.mul(v,p)
            assert LL(A.add(shifted(vh,1),A.scale(vh,-1)))==0
            counts["pearson"]+=1

    # Formal Laurent identity v(y+1)[S(y+1)-S(y)]=2.
    delta=[F(0)]*25
    for k in range(1,25):
        shifted_coefficient=-sum(ell[j]*(-1)**(k-j-1)*comb(k-1,k-j-1)
                                  for j in range(k))
        delta[k]=shifted_coefficient+ell[k-1]
    for r in range(21):
        result=delta[r+3]+3*delta[r+2]+(3-sigma)*delta[r+1]+(1-sigma)*delta[r]
        assert result==(2 if r==0 else 0)

    for x in [F(1,3),F(2),F(7,2)]:
        for y in [F(-2,3),F(5,4),F(-4)]:
            d=x-y
            if d in (0,1,-1):continue
            vx=lambda q:q*(q*q-sigma)
            adj=vx(x)/(d-1)+vx(x+1)/(d+1)-(vx(x)+vx(x+1))/d
            dual=vx(y+1)*(1/(d-1)-1/d)+vx(y)*(1/(d+1)-1/d)
            assert adj-dual==-1
            counts["kernel"]+=1

    moments=[]
    power=[F(1)]
    for k in range(8):
        moments.append(LL(A.mul([F(1,2),F(1)],power)))
        power=A.mul(power,[F(-1,4),F(-1),F(-1)])
    def mu(p):return sum(c*moments[j] for j,c in enumerate(p))
    DP=[[F(1)]];dn=[F(1)]
    for k in range(1,7):
        DP.append(A.mul(DP[-1],[F((2*k-1)**2,4),F(1)]))
        dn.append(dn[-1]*(k*k-sigma))
    Phi=[];Pi=[];Theta=[]
    for n in range(7):
        p=[F(0)]
        for k in range(n+1):
            p=A.add(p,A.scale(DP[k],F((-1)**k*comb(n,k)*comb(n+k,k))/dn[k]))
        Phi.append(p)
        pp=compose(p,[F(-1,4),F(-1),F(-1)])
        Pi.append(pp);Theta.append(theta(pp))
        assert compose(pp,[F(-1),F(-1)])==pp and LL(pp)==0
        assert L(pp)==A.scale(A.mul([F(1),F(2)],pp),2*n*(n+1))

    tail=[F(0)]
    for m in range(1,7):tail.append(tail[-1]+1/vv(m))
    FF={}
    for n in range(7):
        for m in range(7):
            sample=-F((2*m+1)**2,4)
            K=A.value(Pi[n],F(m))
            assert A.value(Phi[n],sample)==K
            direct=mu(divide_linear(A.add([K],A.scale(Phi[n],-1)),sample))
            associated=A.value(Theta[n],F(m))
            assert direct==associated
            FF[n,m]=(K,-K*tail[m]-associated/2)
            counts["associated"]+=1
    for n in range(7):
        for m in range(7):assert FF[n,m]==FF[m,n]
        assert FF[n,0]==(F(1),-tail[n])
    for n in range(6):
        for m in range(6):
            ln=fadd(fscale(vv(n+1),FF[n+1,m]),
                    fscale(-vv(n+1)-vv(n)-2*(2*n+1)*m*(m+1),FF[n,m]),
                    fscale(vv(n),FF[n-1,m]) if n else (F(0),F(0)))
            lm=fadd(fscale(vv(m+1),FF[n,m+1]),
                    fscale(-vv(m+1)-vv(m)-2*(2*m+1)*n*(n+1),FF[n,m]),
                    fscale(vv(m),FF[n,m-1]) if m else (F(0),F(0)))
            assert ln==(F(0),F(-int(n==0)))
            assert lm==(F(0),F(-int(m==0)))
            counts["grid"]+=2
    for n in range(1,7):
        c=2*(3*n*n-sigma)/(n*n-sigma)
        assert fadd(FF[n,n],FF[n-1,n-1],fscale(-c,FF[n,n-1]))==(0,0)
        counts["adjacent"]+=1
    for n in range(6):
        alpha=(n+1)*((n+1)**2-sigma)**2/(2*(3*(n+1)**2-sigma))
        gamma=n*(n*n-sigma)**2/(2*(3*n*n-sigma)) if n else F(0)
        beta=(2*n+1)*(3*n*(n+1)+1-sigma)-alpha-gamma
        residual=fadd(fscale(alpha,FF[n+1,n+1]),fscale(-beta,FF[n,n]),
                      fscale(gamma,FF[n-1,n-1]) if n else (F(0),F(0)))
        assert residual==(0,-int(n==0))
        counts["diagonal"]+=1
    assert -FF[1,1][1]==2*(3-sigma)/(1-sigma)**2

print("PASS: exact polynomial Pearson/second-kind commutators, including non-even test polynomials")
print("PASS: formal Laurent normalization and rational adjoint-kernel constant -1")
print("PASS: direct associated numerators, both sourced grids, symmetry, adjacent balance")
print("PASS: full diagonal recurrence tests and the nonhomogeneous index-zero source")
print(counts)
