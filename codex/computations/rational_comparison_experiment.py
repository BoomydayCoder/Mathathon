"""Bounded exact n=1 rational Hardy Gram comparator, with Arb scalar certificate.
Requires python-flint (installed in /workspace/scratch/odd-zeta-python for this run).
No asymptotic claim. The reported bound includes the explicit S_K, not the
prime-arithmetic multiplier m_(K,M), whose asymptotic theorem has huge thresholds.
"""
from flint import fmpq,fmpz_poly,fmpq_mat,arb,ctx
from math import factorial,comb,prod
from pathlib import Path
import json,hashlib,time

ctx.prec=192
K,N,h,r,L,a=40,3,37,6,51,41
start=time.time()
# Canceled base rational function DN^5/[Dtail*(t+a²)^L].
num=fmpz_poly([1]);den=fmpz_poly([1])
for j in range(1,N+1):num*=fmpz_poly([j*j-a*a,1])**(r-1)
for j in range(N+1,K+1):den*=fmpz_poly([j*j-a*a,1])
bs=[]
for k in range(L):
 value=fmpq(num[k])-sum((fmpq(den[j])*bs[k-j] for j in range(1,min(k,den.degree())+1)),fmpq(0))
 bs.append(value/den[0])
# J_p=(1/pi) integral_0^inf (y²+a²)^(-p)dy.
J=[None]+[fmpq(comb(2*p-2,p-1),2*4**(p-1)*a**(2*p-1)) for p in range(1,L+1)]
Q=[sum((bs[i]*J[L-i-j] for i in range(L-j)),fmpq(0)) for j in range(L)]
simple=[]
for j in range(N+1,K+1):
 numerator=prod(i*i-j*j for i in range(1,N+1))**(r-1)
 denominator=prod(i*i-j*j for i in range(N+1,K+1) if i!=j)*(a*a-j*j)**L
 simple.append((j,fmpq(numerator,denominator)))
mom=[]
for d in range(2*h-1):
 value=sum((c*(-j*j)**d/(2*j) for j,c in simple),fmpq(0))
 value+=sum((fmpq(comb(d,j))*(-a*a)**(d-j)*Q[j] for j in range(min(d,L-1)+1)),fmpq(0))
 assert value>0
 mom.append(value)
print('73 exact rational moments done in',round(time.time()-start,2),'seconds',flush=True)
H=fmpq_mat([[mom[i+j] for j in range(h)] for i in range(h)])
det=H.det()
assert det>0
print('Exact37x37 determinant done in',round(time.time()-start,2),'seconds',flush=True)
# Rigorous scalar certificate: w5(y)*(1+y²/a²)^L < 7/20 for all y>=0.
pi=arb.pi();c5=(2*pi)**4/12
bound=arb(fmpq(7,20));delta=fmpq(1,100);Y=5
small=(2*pi*arb(delta)).exp()/pi*(1+arb(delta)**2/a**2)**L
assert small<bound
# Interval bounds use monotonic factors separately, avoiding interval dependency issues.
maxupper=small.upper();cells=0
for i in range(100,50000):
 lo=arb(fmpq(i,10000));hi=arb(fmpq(i+1,10000))
 q=(-2*pi*lo).exp()
 upper=c5*hi**5*q*(1+11*q+11*q*q+q**3)/(1-q)**5*(1+hi**2/a**2)**L
 assert upper<bound,(i,str(upper))
 if upper.upper()>maxupper:maxupper=upper.upper()
 cells+=1
# For y>=5, sum ell^4 exp(-2pi ell y) <= 2 exp(-2pi y).
q=(-2*pi).exp()
assert (1+11*q+11*q*q+q**3)/(1-q)**5<2
beta=2*pi-fmpq(L,a)
assert beta>0 and 5*beta>5
# log(1+u²)<=u supplies damping; y^5exp(-beta*y) decreases past5.
tail=2*c5*arb(Y)**5*(-beta*Y).exp()
assert tail<bound
C=fmpq(7*a**(2*L),20)
# Explicit determinant normalization from Fauzan equation2.5.
S=fmpq(factorial(K)**(2*h)*4**(h-1),factorial(N)**(12*h)*prod(factorial(2*i)**2 for i in range(1,h)))
logdet=arb(det).log();logupper=h*(pi.log()+arb(C).log())+logdet
logscaled=arb(S).log()+logupper
logbestlower=arb(S).log()+2*L*h*arb(a).log()+logdet
assert logbestlower>2015
raw=Path(__file__).with_suffix('.exact.txt')
raw.write_text('det_H='+str(det)+'\nS_K='+str(S)+'\nC_bound='+str(C)+'\n')
result=dict(parameters=dict(K=K,N=N,h=h,r=r,L=L,a=a),scalar_normalized_bound='7/20',interval_cells=cells,largest_interval_upper=str(maxupper),small_region_upper=str(small),tail_upper=str(tail),log_det_H=str(logdet),log_raw_determinant_upper=str(logupper),log_S_normalized_determinant_upper=str(logscaled),log_best_scalar_comparison_lower=str(logbestlower),exact_data_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),elapsed_seconds=time.time()-start)
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
