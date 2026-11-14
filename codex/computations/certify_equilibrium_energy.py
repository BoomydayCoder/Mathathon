"""Rigorous scalar energy enclosure for the proposed exact equilibrium.
Requires python-flint. Endpoint box is supplied by independently certified
face signs in certify_equilibrium_endpoints.py. No sampled values are proofs.
"""
from flint import arb,acb,fmpq,ctx
from pathlib import Path
import json,time
ctx.prec=160
start=time.time()
a=arb('0.000071195125','0.000000000045')
b=arb('0.75394962','0.00000005')
alpha=arb(fmpq(3,40));lam=arb(fmpq(37,40));pi=arb.pi();d=b-a
A=acb(a);B=acb(b);D=B-A
nan=acb('nan')

def V(x,analytic=False):
 if analytic and not x.real>0:return nan
 w=x.sqrt(analytic=analytic)
 def term(u):
  v=acb(u)/w
  if analytic and not v.real>0:return nan
  return u*(x+u*u).log(analytic=analytic)-2*u+2*w*v.atan()
 return 2*pi*w+term(arb(1))-6*term(alpha)

# Integral0^T du/[(u²+x)sqrt((u²+a)(u²+b))], Carlson form.
# The RF arguments are positive real constants; RJ's fourth argument has Re>0.
pre=[]
for T in [None,arb(1),alpha]:
 if T is None:
  sin=acb(1);cos2=acb(0)
 else:
  sin=acb(T)/(A+T*T).sqrt();cos2=A/(A+T*T)
 sin2=sin*sin;second=1-(1-A/B)*sin2
 RF=acb.elliptic_rf(cos2,second,1)
 pre.append((sin,cos2,sin2,second,RF))

def Jpart(x,entry):
 sin,c2,s2,second,rf=entry
 fourth=c2+(A/x)*s2
 if not fourth.real>0:return nan
 rj=acb.elliptic_rj(c2,second,1,fourth)
 return (sin*rf-(A/x)*sin**3*rj/3)/(x*B.sqrt())

def J(x):
 if not x.real>0:return nan
 return 2*Jpart(x,pre[0])+Jpart(x,pre[1])-6*Jpart(x,pre[2])

def density_integrand(theta,analytic):
 sn=theta.sin();cs=theta.cos();x=A+D*sn*sn
 if not x.real>0:return nan
 return -D*D/pi*sn*sn*cs*cs*V(x,analytic)*J(x)

# Split near the hard edge to keep complex integration disks narrow.
knots=[acb(0),acb('0.01'),acb('0.03'),acb('0.1'),acb('0.3'),acb('0.7'),acb(pi/2)]
parts=[]
for lo,hi in zip(knots,knots[1:]):
 val=acb.integral(density_integrand,lo,hi,abs_tol=arb('0.000002'),rel_tol=arb('0.000002'),eval_limit=20000,depth_limit=30)
 print('segment',lo,hi,val,flush=True)
 assert val.is_finite()
 parts.append(val)
Iv=sum(parts,acb(0))
# Closed plateau constant, derived from endpoint equations.
V0=-2+12*alpha*(1-alpha.log())
def kap(u):return (4*u*u/((u*u+a).sqrt()+(u*u+b).sqrt())**2).log()
M=2*lam*(d/4).log()-V0+kap(arb(1))-6*alpha*kap(alpha)-4*lam
Cstar=-2*lam+12*alpha*lam*(1-alpha.log())+3*lam*lam-2*lam*lam*(2*lam).log()
objective=lam*M/2-Iv.real/2+Cstar
print('M',M,'intVrho',Iv,'Cstar',Cstar,'objective',objective,flush=True)
assert Iv.imag.contains(0)
assert objective < arb('-1.38')
result=dict(endpoint_a=str(a),endpoint_b=str(b),M=str(M),integral_V_rho=str(Iv),Cstar=str(Cstar),objective=str(objective),verified_upper_target='-1.38',segments=[str(x) for x in parts],elapsed_seconds=time.time()-start,precision_bits=ctx.prec)
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS rigorous energy upper bound < -1.38')
