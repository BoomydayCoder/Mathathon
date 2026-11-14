"""Exact checks supplementing the all-M boundary Pearson proof.

The endpoint/Apéry tests at the bottom are explicitly conjectural.
"""
import sympy as s
from explore_rank_one_subspaces import beta_moments
x,t,Z=s.symbols('x t Z')
def ell(M,p):
    den=s.prod(t+j*j for j in range(1,M+1))
    pp=sum(c*(-t)**(k//2) for (k,),c in s.Poly(s.expand(p),x).terms())
    pol,rem=s.div(pp,den,t)
    ans=sum(c*(-1)**k*s.Rational(2*k+3,2)*s.bernoulli(2*k+2)
            for (k,),c in s.Poly(pol,t).terms())
    for j in range(1,M+1):
        w=s.cancel(rem.subs(t,-j*j)/s.diff(den,t).subs(t,-j*j))
        ans+=w*(j*j*(Z-s.harmonic(j,3))-s.Rational(1,2)+s.Rational(1,2*j))
    return s.factor(ans)
def Dop(M,p):
    aa=(x+1)**3*(x-M)/(2*x*(2*x+1))
    return s.cancel(aa*(p.subs(x,x+1)-p)+aa.subs(x,-x)*(p.subs(x,x-1)-p))
count=0
for M in range(1,11):
    cc=s.Rational(M,2*s.factorial(M)**2)
    for k in range(1,7):
        phi=s.prod(x*x-j*j for j in range(1,k+1))
        old=s.prod(x*x-j*j for j in range(1,k))
        assert s.cancel(Dop(M,phi)-k*(k+2-M)*phi-k*(k+1)**2*(k-M)*old)==0
        assert ell(M,Dop(M,phi))==cc*(phi.subs(x,0)-phi.subs(x,1))
        count+=1
    for r,j in [(0,1),(1,2),(2,3),(1,4),(3,5)]:
        p=x**(2*r)+s.Rational(2,3);q=x**(2*j)-s.Rational(5,7)
        assert ell(M,p*Dop(M,q)-q*Dop(M,p))==-cc*(p.subs(x,0)*q.subs(x,1)-p.subs(x,1)*q.subs(x,0))
        count+=1
    for k in range(M-2):
        q=s.prod(x*x-j*j for j in range(1,k+1))/((-1)**k*s.factorial(k)*s.factorial(k+1))
        formula=cc*s.Rational(k+1,(M-k-2)*(M-k-1))*sum(s.Rational(M-r-1,r*(r+1)**2) for r in range(k+1,M-1))
        closed=cc*(M-k-2-M*(k+1)*(s.harmonic(M-1,2)-s.harmonic(k+1,2)))/((M-k-2)*(M-k-1))
        assert s.factor(formula-closed)==0 and ell(M,q)==formula
        count+=1
print('Passed',count,'exact all-M identity checks.',flush=True)
A=[s.Integer(1),s.Integer(5)]
for n in range(1,9):A.append(s.cancel(((2*n+1)*(17*n*n+17*n+5)*A[-1]-n**3*A[-2])/(n+1)**3))
for h in range(1,10):
    _,mom=beta_moments(3,h)
    G=s.Matrix(h,h,lambda i,j:s.Rational(mom[i+j].numerator,mom[i+j].denominator))
    v=[s.Integer(1)] if h==1 else list(-G[:h-1,:h-1].inv()*G[:h-1,h-1])+[s.Integer(1)]
    # Construct Q_h(x^2)=(-1)^(h-1) p_h(-x^2) explicitly.
    p=s.expand((-1)**(h-1)*sum(a*(-x*x)**k for k,a in enumerate(v)))
    cc=s.Rational(h,s.factorial(2*h)**2)
    for k in range(h-1):
        q=x**(2*k)
        assert ell(2*h,(Dop(2*h,p)+(h-1)**2*p)*q)==cc*(p.subs(x,0)*q.subs(x,1)-p.subs(x,1)*q.subs(x,0))
    ah=5*A[h]-A[h-1]
    assert p.subs(x,0)==s.Rational(12*s.factorial(2*h),h*h)/ah
    assert p.subs(x,1)/p.subs(x,0)==s.Rational(h*h*(h+1),2)
    if h>=2:
        rr=p/p.subs(x,0)
        a1=s.cancel((rr.subs(x,2)-rr.subs(x,1))/3)
        assert a1==s.Rational(h*h*(h+1)*(h-1)*(h**3+2*h*h+3*h-2),16*(2*h-1))
    print('h',h,'inhomogeneous equation passed; finite endpoint patterns passed.',flush=True)
print('Endpoint patterns and the original all-h Apéry identity remain conjectural.')
