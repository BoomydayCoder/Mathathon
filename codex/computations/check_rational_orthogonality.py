"""Exact finite verification of rational Gram formulas and Legendre bridge."""
import sympy as S
z=S.symbols('z')

def coefficients(f,poles):
    out=[]
    for a,order in poles:
        for p in range(1,order+1):
            c=S.residue(f*(z+a)**(p-1),z,-a)
            if c:out.append((a,p,c))
    assert S.cancel(f-sum(c/(z+a)**p for a,p,c in out))==0
    return out

def gram_entry(cs,ds):
    return S.factor(sum(c*d*S.binomial(p+q-2,p-1)/(a+b)**(p+q-1) for a,p,c in cs for b,q,d in ds))

basis=[]
for n in range(5):
    f=S.prod((z-S.Rational(2*j-1,2))/(z+S.Rational(2*j-1,2)) for j in range(1,n+1))/(z+S.Rational(2*n+1,2))
    mellin=(-1)**n*S.prod(z+S.Rational(1,2)-j for j in range(1,n+1))/S.prod(z+S.Rational(1,2)+j for j in range(n+1))
    assert S.cancel(mellin-(-1)**n*f)==0
    basis.append(coefficients(f,[(S.Rational(2*j-1,2),1) for j in range(1,n+2)]))
for i in range(5):
 for j in range(5):
    assert gram_entry(basis[i],basis[j])==(S.Rational(1,2*i+1) if i==j else 0)
# A proper damped comparison: K=2,N=0,h=2,r=0,a=1,L=1.
g=[1/((z+1)**2*(z+2)),-z*z/((z+1)**2*(z+2))]
c=[coefficients(f,[(S.Integer(1),2),(S.Integer(2),1)]) for f in g]
H=S.Matrix(2,2,lambda i,j:gram_entry(c[i],c[j]))
assert H.det()>0
print('PASS half-integer rational orthogonality and exact Mellin correspondence n=0..4')
print('Damped comparison Gram:',H,'det=',H.det())
