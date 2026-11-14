#!/usr/bin/env python3
"""Independent exact audit checks for the original Schur/Apéry identity.

No imports from the author's moment or Heun/RII checkers. Moments come from
derivatives of the complete pole polynomial; endpoint algebra is independently
derived from monomial leading coefficients and evaluation at x=1.
"""
import sympy as s

h,k,M=s.symbols('h k M')
x,z,t=s.symbols('x z t')
n=h-1
# Derive the top two D_M coefficients from the Laurent expansion of A_M at infinity.
a2=s.Rational(1,4)
a1=(5-2*M)/8
a0=(7-10*M)/16
am1=(1-14*M)/32
binprod=lambda r:s.prod(2*k-j for j in range(r))/s.factorial(r)
lead=s.factor(2*a2*binprod(2)+2*a1*binprod(1))
sublead=s.factor(2*a2*binprod(4)+2*a1*binprod(3)+2*a0*binprod(2)+2*am1*binprod(1))
assert s.expand(lead-k*(k+2-M))==0
kap=1-2*n*n
eta=s.factor((sublead.subs(k,n)+sublead.subs(k,n+1)-lead.subs(k,n)).subs(M,2*h))
assert s.expand(eta-(1-2*h-2*h**3*(h+1)))==0
yy=s.factor((1-2*h-eta)/(4*h))
zz=s.factor(yy*(3*(eta+kap)+4*(1-2*h))/(16*(1-2*h)))
assert s.factor(yy-h*h*(h+1)/2)==0
assert s.factor(zz-h*h*(h+1)**2*(3*h**3+3*h-2)/(16*(2*h-1)))==0
print('Independent Laurent-coefficient eigenvalue and direct x=0,1 endpoint values: exact.',flush=True)

Y=lambda j:yy.subs(h,j)
Z=lambda j:zz.subs(h,j)
fac=lambda u:(u-(2*h-1)**2)*(u-(2*h)**2)
mat=s.Matrix([[0,1,-fac(0)],
              [Y(h),Y(h),-fac(1)*Y(h-1)],
              [4*Z(h),Z(h),-fac(4)*Z(h-1)]])
det=s.factor(mat.det())
dh=12*h*h-15*h+5
assert s.factor(det-s.Rational(3,8)*h**6*(h+1)**4*dh)==0
sol=[s.factor(q) for q in mat.inv()*s.Matrix([1,Y(h+1),Z(h+1)])]
a,b,c=sol
ph=204*h**5+51*h**4-134*h**3-43*h*h+30*h+12
assert s.factor(a-ph/(h*h*(h+1)**2*(2*h+1)*dh))==0
assert s.factor(c-(h-1)*(12*(h+1)**2-15*(h+1)+5)/(4*h*h*(h+1)*(2*h-1)*(2*h+1)*dh))==0
alpha=s.factor(a*h*h*(2*h+2)*(2*h+1)/(h+1)**2)
gamma=s.factor(c*(h-1)**2*(2*h-1)*(2*h)*(2*h+1)*(2*h+2)/(h+1)**2)
u,v=s.symbols('u v')
K=lambda j:(2*j+1)*(17*j*j+17*j+5)
up=(K(h)*u-h**3*v)/(h+1)**3
um2=(K(h-1)*v-h**3*u)/(h-1)**3
assert s.factor(5*up-u-alpha*(5*u-v)+gamma*(5*v-um2))==0
print('Free-h endpoint determinant, independently solved coefficients, and factorial/Apéry elimination: exact.',flush=True)

N=6
AA=[s.Integer(1),s.Integer(5)]
BB=[s.Integer(0),s.Integer(6)]
for j in range(1,N):
    AA.append(s.cancel((K(j)*AA[-1]-j**3*AA[-2])/(j+1)**3))
    BB.append(s.cancel((K(j)*BB[-1]-j**3*BB[-2])/(j+1)**3))
RR={}
SS={}
qq={}
for H in range(1,N+1):
    D=s.prod(t+j*j for j in range(1,2*H+1))
    beta=lambda j:-s.Integer(j*j)*s.harmonic(j,3)-s.Rational(1,2)+s.Rational(1,2*j)
    const=[]
    for e in range(2*H-1):
        weights=[s.Rational((-j*j)**e,s.diff(D,t).subs(t,-j*j)) for j in range(1,2*H+1)]
        xc=sum(w*j*j for j,w in enumerate(weights,1))
        assert xc==(-1 if e==2*H-2 else 0)
        const.append(sum(w*beta(j) for j,w in enumerate(weights,1)))
    G=s.Matrix(H,H,lambda i,j:const[i+j])
    for size in range(1,H+1):
        assert G[:size,:size].det()>0
    if H==1:
        pc=[s.Integer(1)]
        q=G[0,0]
    else:
        pc=list(G[:-1,:-1].inv()*(-G[:-1,-1]))+[s.Integer(1)]
        q=s.factor(G[-1,-1]+(G[-1,:-1]*s.Matrix(pc[:-1]))[0])
    Q=s.Poly((-1)**(H-1)*sum(pc[j]*(-z)**j for j in range(H)),z)
    assert Q.LC()==1 and Q.eval(0)>0
    assert Q.count_roots(-s.oo,0)==H-1
    R=s.expand(Q.as_expr()/Q.eval(0))
    Rpoly=s.Poly(R,z)
    S=s.interpolate([(j*j,beta(j)*R.subs(z,j*j)) for j in range(1,H+2)],z)
    for j in range(1,2*H+1):
        assert s.cancel(S.subs(z,j*j)-beta(j)*R.subs(z,j*j))==0
    assert s.Poly(S,z).degree()<=H
    assert s.Poly(S,z).coeff_monomial(z**H)==-q*Rpoly.LC()
    assert R.subs(z,1)==Y(H) and R.subs(z,4)==Z(H)
    RR[H]=R
    SS[H]=S
    qq[H]=q
    assert q==s.cancel((5*BB[H]-BB[H-1])/(5*AA[H]-AA[H-1]))
    assert Rpoly.LC()==s.Rational(H*H,12*s.factorial(2*H))*(5*AA[H]-AA[H-1])
    if H>=3:
        j=H-1
        old1,old2=(2*j-1)**2,(2*j)**2
        ratio1=RR[H].subs(z,old1)/RR[j].subs(z,old1)
        ratio2=RR[H].subs(z,old2)/RR[j].subs(z,old2)
        ai=s.cancel((ratio2-ratio1)/(old2-old1))
        bi=s.cancel(ratio1-ai*old1)
        quotient,remainder=s.div(s.Poly(RR[H]-(ai*z+bi)*RR[j],z),
                                s.Poly((z-old1)*(z-old2),z))
        assert remainder.is_zero
        ci=-quotient.eval(0)
        assert s.expand(quotient.as_expr()+ci*RR[j-1])==0
        assert ai==a.subs(h,j) and bi==b.subs(h,j) and ci==c.subs(h,j)
        assert s.expand(SS[H]-(ai*z+bi)*SS[j]+ci*(z-old1)*(z-old2)*SS[j-1])==0
    print('h=',H,'direct pole Gram, root signs, interpolation, original-node RII, and scalar factors: passed.',flush=True)
assert qq[1]==s.Rational(5,4) and qq[2]==s.Rational(577,480)
assert RR[1]==1 and RR[2]==1+5*z
print('All independent exact audit checks passed.',flush=True)
