"""Exact algebra supporting zudilin-explained.md (requires SymPy).

The recurrence certificates are verified with n,t indeterminate, so these are
polynomial identity checks rather than finite samples. Analysis, convergence,
and denominator arguments are proved in the exposition, not by this script.
Run: python codex/computations/check_zudilin.py
"""
import sympy as s

n,t=s.symbols('n t')
mid=(2*n+1)*(17*n*n+17*n+5)
short=lambda z:4*(2*n+1)*(-2*z*z+z+(2*n+1)**2)
short_up=(t-n-1)**2/(t+n+1)**2
short_down=(t+n)**2/(t-n)**2
short_shift=t**4/((t-n)**2*(t+n+1)**2)
assert s.cancel((n+1)**3*short_up-mid+n**3*short_down
                -short(t+1)*short_shift+short(t))==0
print('PASS: arithmetic rational-function recurrence certificate')

P=(-t**6-(8*n-1)*t**5+(4*n**2+27*n+5)*t**4
   +2*n*(67*n**2+71*n+15)*t**3
   +(358*n**4+339*n**3+76*n**2-7*n-3)*t**2
   +(384*n**5+396*n**4+97*n**3-29*n**2-17*n-2)*t
   +n*(153*n**5+183*n**4+50*n**3-30*n**2-22*n-4))
cert=P/((2*t+n)*(t+2*n-1)*(t+2*n))
up=((n+1)**2*(2*t+n+1)*(t-n-1)*(t+2*n+1)*(t+2*n+2)
    /((2*t+n)*(t+n+1)**5))
down=((2*t+n-1)*(t+n)**5
      /(n**2*(2*t+n)*(t-n)*(t+2*n-1)*(t+2*n)))
shift=((2*t+n+2)*t**5*(t+2*n+1)
       /((2*t+n)*(t-n)*(t+n+1)**5))
assert s.cancel((n+1)**3*up-mid+n**3*down
                -cert.subs(t,t+1)*shift+cert)==0
print('PASS: Ball rational-function recurrence certificate')

x=s.symbols('x')
root2=s.sqrt(2)
alpha=(root2-1)**4
critical_poly=x**5*(x+2)-(x-1)*(x+1)**5
factored=-(x+s.Rational(1,2))*(2*(x+s.Rational(1,2))**4
             -5*(x+s.Rational(1,2))**2-s.Rational(7,8))
assert s.expand(critical_poly-factored)==0
# The critical point is the positive root of x²+x=1+sqrt(2).
# This remainder proves the stated value of exp(f) at that point.
value_poly=(x+2)**2*(x-1)-alpha*(x+1)**5
assert s.simplify(s.rem(value_poly,x*x+x-(1+root2),x))==0
print('PASS: saddle critical polynomial and exact maximum value')

# Ball n=1 partial fractions and telescoping initial values.
ball_one=(2*t+1)*(t-1)*(t+2)/(t**4*(t+1)**4)
ball_one_partial=(-2/t**4+2/(t+1)**4+5/t**3+5/(t+1)**3
                  -5/t**2+5/(t+1)**2)
assert s.cancel(ball_one-ball_one_partial)==0
print('PASS: Ball n=1 partial fractions')

# h(x+30)=h(x) because the sum of the shifts is zero. Each floor is
# constant between adjacent integers, so the table covers every real x.
assert 30-15-10-6+1==0
weights=[m-m//2-m//3-m//5+m//30 for m in range(30)]
assert all(v in (0,1) for v in weights)
assert weights[1:6]==[1]*5
# Exact sufficient exponent inequality from the exposition. It uses
# 1+sqrt(2)>12/5, since 2>49/25.
assert 2*s.Integer(25)>49
assert 5**115 < 2**158 * 3**73
print('PASS: lcm floor staircase and sufficient exponential inequality')

# Simplified analytic tail bound: write t/n=y>1. The inequality
# (2t+n)/(t+n)<=2 is just 2t+n<=2(t+n), and t+2n<=3t.
# These are proved in prose; no finite sampling is used to assert them.
print('All exact checks passed. SymPy',s.__version__)
