"""Exact recurrence/growth checks and illustrative Apéry asymptotic values."""
from fractions import Fraction as F
from math import comb
import sympy as S
import mpmath as mp

m,n=S.symbols('m n')
b=34*n**3+51*n*n+27*n+5
assert S.expand(5*(b-n**3/S.Integer(5)-30*(n+1)**3).subs(n,m+12))==19*m**3+489*m*m+3213*m+847
assert S.expand(5*b-n**3-25*(n+1)**3)==144*n**3+180*n*n+60*n
Q=F(135,4)
assert S.cancel(b/(n+1)**3-(34-51/(n+1)+(27*n+22)/(n+1)**3))==0
assert F(169,5)-F(51,1020)==Q
assert Q**4/F(2**6*3**9)==F(16875,16384)>1
h=[j-j//2-2*(j//3)+j//6 for j in range(6)]
assert h==[0,1,1,0,0,1]
A=[sum(comb(j,k)**2*comb(j+k,k)**2 for k in range(j+1)) for j in range(32)]
B=[F(0),F(6)]
for j in range(1,31):
    bj=34*j**3+51*j*j+27*j+5
    assert (j+1)**3*A[j+1]-bj*A[j]+j**3*A[j-1]==0
    B.append((bj*B[-1]-j**3*B[-2])/(j+1)**3)
    assert A[j]*B[j+1]-A[j+1]*B[j]==F(6,(j+1)**3)
print('PASS exact recurrence, Casoratian, growth polynomials, floor table and logarithm-free comparison')
mp.mp.dps=80
lam=(1+mp.sqrt(2))**4
CA=(1+mp.sqrt(2))**2/(2**mp.mpf('2.25')*mp.pi**mp.mpf('1.5'))
print('n, A_n*n^(3/2)/lambda^n, ratio to predicted constant')
for j in (10,30,100,300,1000):
    a=sum(comb(j,k)**2*comb(j+k,k)**2 for k in range(j+1))
    value=mp.mpf(a)*j**mp.mpf('1.5')/lam**j
    assert mp.mpf(a)>=mp.exp(-3)*lam**j/(4*(j+1)**4)
    print(j,mp.nstr(value,18),mp.nstr(value/CA,15))
print('Predicted constant:',mp.nstr(CA,20))
print('Perturbation experiment: b, first nonintegral index through 12')
for parameter in range(1,21):
    U=[F(1),F(parameter)]
    first=None
    for j in range(1,12):
        U.append(((2*j+1)*(17*j*j+17*j+parameter)*U[-1]-j**3*U[-2])/(j+1)**3)
        if U[-1].denominator!=1 and first is None:
            first=j+1
    print(parameter,first)
