"""Tiny exact certificate for the shifted-trial exponential improvement."""
from fractions import Fraction as F
from math import factorial

# Machin's identity plus the first two alternating terms.
pi_lower=16*(F(1,5)-F(1,5)**3/3)-F(4,239)
assert pi_lower>F(157,50)

# This implies3pi/2>21/sqrt20 and hence v>sqrt20 for
# v-atan(v)=pi, because atan(1/sqrt20)<1/sqrt20.
assert F(21**2,20)<(F(3,2)*F(157,50))**2

# The exponential series gives e<8/3+5/96=87/32<11/4.
assert F(8,3)+F(5,96)==F(87,32)<F(11,4)

a=F(3,20)
b=F(3,5)
middle=F(1,441)*(1-a*a/(b*b))**2
assert middle==F(25,12544)
inner=a**4*F(11,4)**4/16

# At x=b, use atan(3/10)<3/10 and pi>157/50.
q=F(323,250)
exp_lower=sum(q**j/F(factorial(j)) for j in range(4))
tail=F(81,11881)/exp_lower
assert inner<middle and tail<middle and middle<F(1,500)

# Optional precise monomial-trial root bracket, still all rational.
def atan_interval(x,n=16):
    center=sum((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n))
    next_sum=center+(-1)**n*x**(2*n+1)/F(2*n+1)
    return min(center,next_sum),max(center,next_sum)

atan5=atan_interval(F(1,5))
atan239=atan_interval(F(1,239))
pi=(16*atan5[0]-4*atan239[1],16*atan5[1]-4*atan239[0])
vlo,vhi=F(4493409,1000000),F(4493410,1000000)
alo=atan_interval(1/vlo)
ahi=atan_interval(1/vhi)
assert vlo-F(3,2)*pi[0]+alo[1]<0
assert vhi-F(3,2)*pi[1]+ahi[0]>0
rho=(1/(1+vhi*vhi)**2,1/(1+vlo*vlo)**2)
assert rho[1]<F(1,441)
print('PASS all exact trial bounds.')
print('Shifted uniform rate bound:',middle,'<1/500.')
print('Monomial rate bracket, decimal display:',float(rho[0]),float(rho[1]))
