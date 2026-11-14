"""Construct the z=0 boundary rational function of the lifted certificate."""
import sympy as s,json
from pathlib import Path
x,y,t=s.symbols('x y t');p=(1+x)*(1+y);ss=x*y
D=s.expand((ss-t*p*(ss+p))**2-4*t*t*p**3*ss)
U=s.Poly.from_dict({(i,j,k):s.Rational(c) for i,j,k,c in json.loads(Path(__file__).with_name('geometric_certificate_U.json').read_text())},(x,y,t)).as_expr();V=U.xreplace({x:y,y:x})
b=ss-t*p*(ss+p)
Dx=s.diff(D,x);Dy=s.diff(D,y);Dt=s.diff(D,t)
T1N=2*ss*D**2
T2N=-3*ss*t*Dt*D
T3N=-4*ss*(t*Dt+t*t*s.diff(D,t,2))*D+s.Rational(15,2)*ss*t*t*Dt**2
TLN=s.expand((1-34*t+t*t)*T3N+(-51*t+3*t*t)*T2N+(-27*t+3*t*t)*T1N)
Bx=s.expand(b*(p+3*x*s.diff(p,x))-2*x*p*s.diff(b,x));By=Bx.xreplace({x:y,y:x})
P=s.expand(p*TLN-(U/x)*Bx-(V/y)*By)
W=s.cancel(P/(p*D**3))
if __name__=='__main__':
 print('boundary numerator terms',len(s.Poly(s.fraction(W)[0],x,y,t).terms()),flush=True)
 print('numerator factor:',s.factor(s.fraction(W)[0]),flush=True)
 print('denominator factor:',s.factor(s.fraction(W)[1]),flush=True)
 Path(__file__).with_name('geometric_boundary_W0.txt').write_text(str(s.factor(W))+'\n')
