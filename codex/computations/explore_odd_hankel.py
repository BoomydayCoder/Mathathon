"""Small exact Hankel experiments. Numeric evaluations are evidence, not bounds."""
from pathlib import Path
import json
import sympy as S
import mpmath as mp

t,X=S.symbols('t X')
mp.mp.dps=150

def functional_parts(s,K,N,r,d):
    # polynomial quotient plus residues after exact cancellation
    num=S.Poly(S.prod(t+j*j for j in range(1,N+1))**r*t**d,t)
    den=S.Poly(S.prod(t+j*j for j in range(1,K+1)),t)
    q,rem=S.div(num,den)
    const=S.Rational(0);lead=S.Rational(0)
    for (e,),coeff in q.terms():
        const+=coeff*(-1)**e*S.bernoulli(2*e+2)*S.factorial(2*e+s)/(S.factorial(s-1)*S.factorial(2*e+2))
    for j in range(1,K+1):
        res=S.cancel(rem.as_expr().subs(t,-j*j)/S.diff(den.as_expr(),t).subs(t,-j*j))
        lead+=res*j**(s-1)
        const+=res*(-j**(s-1)*sum(S.Rational(1,v**s) for v in range(1,j+1))-S.Rational(1,s-1)+S.Rational(1,2*j))
    return lead,const

def run(s,K,N,r):
    h=K-N
    parts=[functional_parts(s,K,N,r,d) for d in range(2*h-1)]
    mat=S.Matrix(h,h,lambda i,j:parts[i+j][0]*X+parts[i+j][1])
    poly=S.Poly(mat.det(method='domain-ge'),X)
    denominator,integral=poly.clear_denoms(convert=True)
    content,primitive=integral.primitive()
    # Rational scaling clears denominators AND removes integer content.
    coefficients=[int(c) for c in primitive.all_coeffs()]
    expected=(-1)**(h*(h-1)//2)*S.prod(j**(s-1)*S.prod(i*i-j*j for i in range(1,N+1))**(r-1) for j in range(N+1,K+1))
    assert poly.degree()==h and poly.LC()==expected
    value=mp.polyval(coefficients,mp.zeta(s))
    assert value>0
    return dict(s=s,K=K,N=N,r=r,h=h,degree=poly.degree(),coefficients=coefficients,
                normalization=str(S.Rational(denominator,content)),
                log_abs_value=str(mp.log(value)),neg_log_value_per_degree=str(-mp.log(value)/h),
                max_coefficient_digits=max(len(str(abs(c))) for c in coefficients))

if __name__=='__main__':
    rows=[]
    for K in range(1,7):
        for N in [0]+([1] if K>=2 else []):
            for s in [3,5,7]:
                row=run(s,K,N,6);rows.append(row)
                print(s,K,N,row['max_coefficient_digits'],float(row['neg_log_value_per_degree']),flush=True)
    target=Path(__file__).resolve().parent/'odd_hankel_results.json'
    target.write_text(json.dumps(dict(precision=mp.mp.dps,rows=rows),indent=2)+'\n')
