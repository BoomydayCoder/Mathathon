"""Nine bounded subspace experiments. Real energies are NOT certified."""
from fractions import Fraction as F
from math import log
from explore_odd_parameters import outer_integral,equilibrium

assert outer_integral(5,6,F(3,40))==F(129101,96000)
assert outer_integral(7,6,F(3,40))==F(27071,16000)
assert outer_integral(5,6,F(3,40),refined=True)==F(127751,96000)
s=7;r=6;m=r//2
for alpha in (F(1,20),F(3,40),F(1,10)):
    for lam in (F(3,5),F(3,4),F(9,10)):
        I=outer_integral(s,r,alpha,lam)
        af=float(alpha);lf=float(lam)
        A=float(I)+(m*m+m+1)/144+lf*(s/3+4*r*af/3+2*r*af*log(af))
        a,b,E=equilibrium(r,af,lf)
        print('alpha',alpha,'lambda',lam,'outer',I,'arithmetic',round(A,9),'real',round(E,9),'total',round(A+E,9),flush=True)

# One motivated follow-up: the best preceding point, with r=s+1.
alpha=F(1,10);lam=F(3,4)
for r in (6,8):
    m=r//2;af=float(alpha);lf=float(lam)
    I=outer_integral(s,r,alpha,lam,refined=True)
    A=float(I)+(m*m+m+1)/144+lf*(s/3+4*r*af/3+2*r*af*log(af))
    _,_,E=equilibrium(r,af,lf)
    print('refined alpha',alpha,'lambda',lam,'r',r,'outer',I,'total',round(A+E,9),flush=True)
