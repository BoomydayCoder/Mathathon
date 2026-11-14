"""Finite rational-polynomial audit of endpoint congruences; not a proof."""
import sympy as s
x=s.symbols('x')
for a,b in [(1,1),(2,3)]:
    # r=1 assumed rational pi=a/b: nonzero subset sums +/- i a/b.
    g=b*b*x*x+a*a
    c=b*b
    m=2
    k=2
    for p in [3,5,7]:
        if k*c*a*a%p==0:
            continue
        f=s.Poly(c**(m*p)*x**(p-1)*g**p/s.factorial(p-1),x)
        cur=f
        F=f
        while cur.degree()>0:
            cur=cur.diff()
            F=F+cur
        F0=F.eval(0)
        assert F0.is_integer
        assert (F0-c**(m*p)*(a*a)**p)%p==0
        remainder=F.rem(s.Poly(g,x))
        root_sum=2*remainder.eval(0)
        assert root_sum.is_integer and root_sum%p==0
        N=k*F0+root_sum
        assert N.is_integer and N%p!=0
assert 3**18*5**135<2**347
print('PASS: finite pi-power endpoint congruences; exact zeta2 lcm-rate inequality')
