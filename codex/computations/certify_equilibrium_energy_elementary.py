"""Independent virial/potential cross-check of the Carlson energy certificate."""
from flint import arb,acb,fmpq,ctx
from pathlib import Path
import json,time
ctx.prec=160
start=time.time()
a=arb('0.000071195125','0.000000000045')
b=arb('0.75394962','0.00000005')
alpha=arb(fmpq(3,40))
lam=arb(fmpq(37,40))
pi=arb.pi()
A=acb(a)
B=acb(b)
D=B-A
nan=acb('nan')
V0=-2+12*alpha*(1-alpha.log())
def kap(u):
    return (4*u*u/((u*u+a).sqrt()+(u*u+b).sqrt())**2).log()
M=2*lam*((b-a)/4).log()-V0+kap(arb(1))-6*alpha*kap(alpha)-4*lam
Cstar=-2*lam+12*alpha*lam*(1-alpha.log())+3*lam*lam-2*lam*lam*(2*lam).log()

def Dt(u):
    return ((acb(u)+A).sqrt()+(acb(u)+B).sqrt())**2/4
def tau(u):
    return D/(4*Dt(u))
t1=tau(1)
ta=tau(alpha*alpha)

def combo_u(u,analytic):
    sa=(u*u+A).sqrt(analytic=analytic)
    sb=(u*u+B).sqrt(analytic=analytic)
    t=D/(sa+sb)**2
    return (1-t1*t).log(analytic=analytic)-6*alpha*(1-ta*t).log(analytic=analytic)

def combo_tail(v,analytic):
    sa=(1+A*v*v).sqrt(analytic=analytic)
    sb=(1+B*v*v).sqrt(analytic=analytic)
    z0=D/(sa+sb)**2
    def term(ts):
        z=ts*z0*v*v
        if analytic and not z.real<1:
            return nan
        return -ts*z0*z.hypgeom_2f1(1,1,2)
    return term(t1)-6*alpha*term(ta)

small=acb.integral(combo_u,0,alpha,abs_tol=arb('0.000001'),rel_tol=arb('0.000001'),eval_limit=10000)
medium=acb.integral(combo_u,0,1,abs_tol=arb('0.000001'),rel_tol=arb('0.000001'),eval_limit=10000)
large=acb.integral(combo_tail,0,1,abs_tol=arb('0.000001'),rel_tol=arb('0.000001'),eval_limit=10000)
Iv2=2*lam*lam+lam*(Dt(1).log()-6*alpha*Dt(alpha*alpha).log())+lam*(-2+12*alpha)+3*medium-6*small+2*large
obj2=lam*M/2-Iv2.real/2+Cstar
assert Iv2.imag.contains(0)
assert obj2 < arb('-1.38')
# Compare with the previously saved first route's outward-rounded real interval.
assert (arb('1.916','0.000305')-Iv2.real).contains(0)
print('Independent potential/virial check:',Iv2,'objective',obj2,flush=True)
second=dict(integral_V_rho=str(Iv2),objective=str(obj2),elementary_integrals=dict(small=str(small),medium=str(medium),tail=str(large)),elapsed_seconds=time.time()-start)
Path(__file__).with_suffix('.json').write_text(json.dumps(second,indent=2)+'\n')
print('PASS independent elementary integral certificate < -1.38')
