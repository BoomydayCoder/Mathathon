"""Independent exact audit of the reciprocal real region and cube map.
No numerical integration and no geometric theorem is used.
"""
import sympy as S
x,y,z,a,b,c,u,v,w=S.symbols('x y z a b c u v w')
V=S.Matrix([x,y,z]); C=(1+x)*(1+y); Q=C+x*y*z
L=C*(1+z)*Q/(x*y*z)
T=S.Matrix([-z/(1+z),-C/Q,-x/(1+x)])
def sub(f,vs,vals): return f.subs(dict(zip(vs,vals)),simultaneous=True)
def zero(f): assert S.cancel(f)==0, S.factor(f)
zero(sub(L,V,T)*L-1)
for e in sub(T,V,T)-V: zero(e)
# Full ordinary determinant, independently of the sparse log-Jacobian formula.
zero(T.jacobian(V).applyfunc(S.cancel).det(method='domain-ge')*x*y*z/S.prod(T)-1/L)
F=-T
R= a*b*c-(1-a)*(1-b)
Inv=S.Matrix([c/(1-c),(1-a)*(1-b)/R,a/(1-a)])
for e in sub(F,V,Inv)-S.Matrix([a,b,c]):zero(e)
for e in sub(Inv,(a,b,c),F)-V:zero(e)
Phi=(1-a)*(1-b)*(1-c)*R/(a*b*c)
zero(sub(Phi,(a,b,c),F)-1/L)
zero(F.jacobian(V).applyfunc(S.cancel).det(method='domain-ge')*x*y*z/S.prod(F)-1/L)
D=1-w+u*v*w
cube=S.Matrix([w,(1-w)/D,u])
zero(sub(Phi,(a,b,c),cube)-u*(1-u)*v*(1-v)*w*(1-w)/D)
zero(cube.jacobian([u,v,w]).applyfunc(S.cancel).det(method='domain-ge')-u*w*(1-w)/D**2)
zero(cube.jacobian([u,v,w]).applyfunc(S.cancel).det(method='domain-ge')/S.prod(cube)-1/D)
zero(sub((1-a)*(1-b)/(a*b*c),(a,b,c),cube)-v)
point=[1/S.sqrt(2),1/S.sqrt(2),1+S.sqrt(2)]
lam=17+12*S.sqrt(2); rho=17-12*S.sqrt(2)
assert S.simplify(sub(L,V,point)-lam)==0
for q in V:assert S.simplify(sub(S.diff(L,q),V,point))==0
image=[1/S.sqrt(2),1/S.sqrt(2),S.sqrt(2)-1]
for e in sub(F,V,point)-S.Matrix(image):assert S.simplify(e)==0
assert S.simplify(sub(Phi,(a,b,c),image)-rho)==0
# Positive Laurent coefficients and exponent support spanning 3-space imply
# strict convexity of L(exp(X),exp(Y),exp(Z)), hence global unique minimum.
exponents=[]
for term in S.Add.make_args(S.expand(L)):
    p=term.as_powers_dict(); e=[int(p.get(q,0)) for q in V]
    coeff=S.cancel(term/S.prod(q**k for q,k in zip(V,e)))
    assert coeff>0
    exponents.append(e)
assert S.Matrix(exponents).rank()==3
# The source-map counterexample, independent direct rational substitution.
phi=(1-x)*(1-y)*(1-z)*((1-x)*(1-y)-x*y*z)/(-x*y*z)
printed=S.Matrix([z/(1-z),-(1-x)*(1-y)/((1-x)*(1-y)-x*y*z),x/(1-x)])
assert sub(sub(phi,V,printed),V,[-1,-1,-1])==S.Rational(2349,40)
assert sub(phi,V,[-1,-1,-1])==40
print('PASS: full Jacobians, involution, both real inverses, cube pullback, critical point, strict convexity support, printed-map counterexample')
