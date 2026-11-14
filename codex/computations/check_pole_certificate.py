"""Exact symbolic checks of the independent pole-matching derivation.
Requires SymPy. Symbolic identities are stronger than finite numerical tests.
The mathematical proof is in problem131/simplification.md.
"""
import sympy as s
n,t=s.symbols('n t')
a,b,c,A,B,C=s.symbols('a b c A B C')
q=lambda x:a*x*x+b*x+c
D=(t-n)**2*(t+n+1)**2
polynomial=s.Poly(s.expand(
    A*(t-n-1)**2*(t-n)**2+B*D+C*(t+n)**2*(t+n+1)**2
    -q(t+1)*t**4+q(t)*D),t)
M,_=s.linear_eq_to_matrix(polynomial.all_coeffs(),[A,B,C,a,b,c])
nullspace=M.nullspace()
assert len(nullspace)==1
normalized=[s.factor(x/nullspace[0][0]*(n+1)**3) for x in nullspace[0]]
expected=[(n+1)**3,-(2*n+1)*(17*n*n+17*n+5),n**3,
          -8*(2*n+1),4*(2*n+1),4*(2*n+1)**3]
assert all(s.factor(x-y)==0 for x,y in zip(normalized,expected))
assert s.expand(polynomial.as_expr().subs(dict(zip([A,B,C,a,b,c],expected))))==0
pole_equations=[n*(2*n+1)*s.diff(q(t),t).subs(t,n+1)+(4*n+3)*q(n+1),
                (n+1)*(2*n+1)*s.diff(q(t),t).subs(t,-n)-(4*n+1)*q(-n)]
P,_=s.linear_eq_to_matrix(pole_equations,[a,b,c])
assert s.factor(P[:,[0,2]].det())==-3*(2*n+1)**3
assert all(s.factor(e.subs({a:-2,b:1,c:(2*n+1)**2}))==0 for e in pole_equations)
print('PASS: both pole constraints, rank determinant, unique normalized certificate, full symbolic identity')
print('Recurrence/certificate coefficients:',normalized)
