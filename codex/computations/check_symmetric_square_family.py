"""Exact symbolic proof check for arbitrary recurrence parameters."""
import sympy as S
t,a,b=S.symbols('t a b'); u=S.Function('u')(t)
def th(f):return t*S.diff(f,t)
inner=a*th(th(u))+a*th(u)+b*u
L=th(th(th(u)))-t*(2*th(inner)+inner)+t*t*(th(th(th(u)))+3*th(th(u))+3*th(u)+u)
lead=S.expand(L).coeff(S.diff(u,t,3))
c=[S.cancel(S.expand(L).coeff(S.diff(u,t,j))/lead) for j in range(3)]
p=(1-3*a*t+2*t*t)/(t*(1-2*a*t+t*t))
q=(t-2*b)/(4*t*(1-2*a*t+t*t))
assert S.cancel(c[2]-3*p)==0
assert S.cancel(c[1]-S.diff(p,t)-2*p*p-4*q)==0
assert S.cancel(c[0]-2*S.diff(q,t)-4*p*q)==0
assert S.cancel((3*(2*a+b)*b-1)/8).subs({a:2,b:1})==S.Rational(7,4)
for n in range(1,20):
 assert (n+1)**3==(2*n+1)*(n*(n+1)+1)-n**3
print('PASS arbitrary-parameter symmetric-square identity and contrasting examples')
