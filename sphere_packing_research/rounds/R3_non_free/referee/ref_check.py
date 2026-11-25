#!/usr/bin/env python3
"""Referee checks for R3_non_free (independent re-derivations). Light."""
import math, itertools
import numpy as np
from mpmath import mp, mpf, loggamma, log, pi, exp, sqrt
mp.dps = 30
def logvol(m, R): return (mpf(m)/2)*log(pi) + m*log(R) - loggamma(mpf(m)/2+1)
def h(x): return 0.0 if x<=0 or x>=1 else -x*math.log(x)-(1-x)*math.log(1-x)

print("(A) Prop 2: exact T_U vs asymptotic; and exact 'poly' factor in V = T^{1/f} (f tau^2)^{n/2} * P")
for (D,sd,k,tau2,V) in [(4,2,50,2.0,1e4),(12,2,30,6.0,1e3),(40,2,45,20.0,1e3),(200,8,60,25.0,1e6)]:
    n=D*k; f=mpf(sd)/D; fn=sd*k
    logR=(log(V)+loggamma(mpf(n)/2+1)-(mpf(n)/2)*log(pi))/n; R=exp(logR)
    cU=mpf(tau2)**(mpf(sd)/2)
    exact=exp(logvol(fn,R)-k*log(cU))
    approx=mpf(V)**f*(pi*n)**(f/2)/sqrt(pi*f*n)*((1/f)/tau2)**(mpf(fn)/2)
    # the factor P := V / (T^{1/f} (f tau^2)^{n/2})  -- notes call it poly(n)
    P = mpf(V)/(exact**(1/f)*(f*tau2)**(mpf(n)/2))
    Ppred = (pi*f*n)**(1/(2*f))/sqrt(pi*n)
    print(f"  D={D} dimU={sd} k={k} n={n} f={float(f):.3f}: exact/approx={float(exact/approx):.5f}; log10 P={float(log(P)/log(10)):.2f} (pred {(float(log(Ppred)/log(10))):.2f}); log10 n^3={3*math.log10(n):.2f}")

print("(B) borderline term T_U = f^{-1/2} pi^{(f-1)/2} V^f n^{-(1-f)/2}: growth for V=n, n^2")
for f in [0.5, 1/3, 0.25]:
    for a in [1,2]:
        print(f"  f={f:.3f}, V=n^{a}: T_U ~ n^{a*f-(1-f)/2:.3f}")

print("(C) entropy windows, own bisection")
def ok(t,q):  # a/b=t
    a=t; b=1.0; m=a+b*(q-1)
    return m*h(a/m) <= 2*min(a,b)*math.log(q)+1e-12
for q in [2,3,5,7,11,101,10007]:
    ts=[i/10000 for i in range(1,40000)]
    good=[t for t in ts if ok(t,q)]
    if good: print(f"  q={q}: a/b in [{min(good):.4f},{max(good):.4f}]  rho in [{(min(good)-1)/(min(good)+q-1):.4f},{(max(good)-1)/(max(good)+q-1):.4f}]  e/q={math.e/q:.4f}")
    else: print(f"  q={q}: none")
print("  q=2 exact: only a=b:", [ (a,b) for a in range(1,30) for b in range(1,30) if (a+b)*h(a/(a+b)) <= 2*min(a,b)*math.log(2)+1e-12 and a!=b])

print("(D) SL2(p), p=1 mod 4, nonsplit torus prime q | (p+1)/2, (a,b) ~ (N-1,N), N=(p+1)/(2q); split torus a=b")
from sympy import primerange, factorint
bad=[]
for p in primerange(5,2000):
    if p%4!=1: continue
    for q in factorint((p+1)//2):
        N=(p+1)//(2*q)
        if N==1: continue
        a=4*(N-1); b=4*N; m=a+b*(q-1)
        assert m==2*(p-1)
        if not (m*h(a/m) <= 2*min(a,b)*math.log(q)+1e-12): bad.append((p,q,N))
    for q in factorint((p-1)//2):
        if q==2: continue
        a=b=2*(p-1)//q; m=a+b*(q-1); assert m==2*(p-1)
        if not (m*h(a/m) <= 2*min(a,b)*math.log(q)+1e-12): bad.append((p,q,'split'))
print("  failures for p<2000:", bad)

print("(E) SL2(p), p=3 mod 4: (i) involution in degree (p-1)/2: Delta; (ii) unipotent (a,b)=(2,1)")
for p in [7,11,19,23,47,103,1003]:
    m=p-1; a=(p-1)//2+1; b=m-a
    D1=(m/2)*h(a/m)-min(a,b)*math.log(2)
    m2=p+1; D2=(m2/2)*h(2/m2)-math.log(p)
    print(f"  p={p}: (i) Delta={D1:.4f} (log2={math.log(2):.4f});  (ii) Delta={D2:.4f} (1-log2={1-math.log(2):.4f})")

print("(F) Unramified example: Q-span vs K-span of G = mu_6 x| C_2 on E=Q(sqrt-3,sqrt5)")
# basis of E over Q: 1, s3, s5, s3*s5  (s3=sqrt(-3), s5=sqrt5). mult by zeta6=(1+s3)/2 ; Galois tau: s3->-s3, s5->-s5 (fixes s3*s5)
def mul_s3():  # matrix of multiplication by s3 on basis (1,s3,s5,s3s5)
    M=np.zeros((4,4)); M[1,0]=1; M[0,1]=-3; M[3,2]=1; M[2,3]=-3; return M
S3=mul_s3(); Z6=(np.eye(4)+S3)/2
TAU=np.diag([1.,-1.,-1.,1.])
els=[]; 
for i in range(6):
    Zi=np.linalg.matrix_power(Z6,i)
    els.append(Zi); els.append(Zi@TAU)
rk=np.linalg.matrix_rank(np.array([e.flatten() for e in els]))
tr=[np.trace(e) for e in els]
print("  |G|=",len(els)," dim_Q span(G) =",rk,"  (M_2(K) would have dim 8);  sum rho^2 =",sum((t/4)**2 for t in tr), " |mu_E|/[E:Q] =",6/4)

print("(G) D4 over Z[sqrt-2]: all K-lines super-stable? min over w of covol(w O_K) with saturation check")
def qmul(a,b):
    a0,a1,a2,a3=a;b0,b1,b2,b3=b
    return (a0*b0-a1*b1-a2*b2-a3*b3,a0*b1+a1*b0+a2*b3-a3*b2,a0*b2-a1*b3+a2*b0+a3*b1,a0*b3+a1*b2-a2*b1+a3*b0)
s=(0,1,1,0); mn=None
for w in itertools.product(range(-3,4),repeat=4):
    if sum(w)%2 or not any(w): continue
    sw=qmul(s,w); c=math.sqrt(sum(x*x for x in w)*sum(x*x for x in sw))
    mn=c if mn is None else min(mn,c)
print("  min covol(w O_K) =",mn," -> tau^2_min <= ",mn/math.sqrt(2)," (c_W=2, tau^2=c_U/sqrt(c_W)); O_K PID so every line lattice is w O_K, c_U=sqrt2|w|^2>=2sqrt2 => tau^2>=2 exactly")

print("(H) Blichfeldt check p=3,5,7,11 and group identity SL2(5)xC7")
for p in [3,5,7,11]:
    m=p*(p-1); rd=p**((p-2)/(p-1)); bl=float((2/pi)*exp(2*loggamma(mpf(2)+mpf(m)/2)/m))
    print(f"  p={p}: need gamma>= {m/rd:.3f}; Blichfeldt gamma_{m} <= {bl:.3f}")
print("  SL2(5)xC7: g = 15*7/6 =",15*7/6, " = 840/(12*4) =",840/48)
