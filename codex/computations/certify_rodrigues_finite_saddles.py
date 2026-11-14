#!/usr/bin/env python3
"""Exact rational Rouche/argument/modulus certificate for multiplicities 2..5.
Sympy approximate roots propose centers ONLY; every accepted claim is then
checked by exact rational inequalities. Prints the certificates for review.
"""
from fractions import Fraction as F
from math import factorial, floor
import json
from pathlib import Path
import sympy as s


def atan_small(x):
    if x<0:
        lo,hi=atan_small(-x);return -hi,-lo
    assert 0<=x<=F(1,2)
    total=F(0)
    for k in range(50):total+=(-1)**k*x**(2*k+1)/(2*k+1)
    remainder=x**101/101
    return total,total+remainder

p1=atan_small(F(1,5));p2=atan_small(F(1,239))
PI=(16*p1[0]-4*p2[1],16*p1[1]-4*p2[0])
assert F(3141592653589793238,10**18)<PI[0]<PI[1]<F(3141592653589793239,10**18)

def atan(x):
    assert x>=0
    if x>1:
        lo,hi=atan(1/x);return PI[0]/2-hi,PI[1]/2-lo
    if x>F(1,2):
        lo,hi=atan_small((x-1)/(x+1));return PI[0]/4+lo,PI[1]/4+hi
    return atan_small(x)

def atan_interval(lo,hi):
    return atan(lo)[0],atan(hi)[1]

def complex_l1(v):
    v=s.expand(v);return abs(F(s.re(v)))+abs(F(s.im(v)))

def complex_lower(v):
    v=s.expand(v);return max(abs(F(s.re(v))),abs(F(s.im(v))))

x=s.symbols('x');radius=F(1,10**7);records=[]
for m in range(2,6):
    roots=[]
    for eps in [1,-1]:
        poly=s.Poly((m-x)*(1+x)**m-eps*x**(m+1),x)
        degree=poly.degree(); centers=[]
        for approximate in s.nroots(poly,n=40,maxsteps=300):
            re,im=approximate.as_real_imag()
            ar=F(int(s.floor(re*10**10+s.Rational(1,2))),10**10)
            ai=F(int(s.floor(im*10**10+s.Rational(1,2))),10**10)
            center=s.Rational(ar.numerator,ar.denominator)+s.I*s.Rational(ai.numerator,ai.denominator)
            main=complex_lower(poly.diff().eval(center))*radius
            rest=complex_l1(poly.eval(center))
            for k in range(2,degree+1):
                rest+=complex_l1(poly.diff((x,k)).eval(center))*radius**k/factorial(k)
            assert rest<main,('Rouche failed',m,eps,ar,ai)
            centers.append((ar,ai))
            roots.append({'eps':eps,'a':ar,'b':ai,'rouche_margin':main-rest})
        assert len(centers)==degree
        for i in range(degree):
            for j in range(i):
                # Stronger coordinate-separation condition suffices for disjoint disks.
                assert max(abs(centers[i][0]-centers[j][0]),abs(centers[i][1]-centers[j][1]))>2*radius
    assert len(roots)==2*m
    components=[]
    for root in roots:
        # Keep t=-x with Im t>=0, excluding positive x (the large period).
        if root['a']>0:continue
        if root['b']>0:continue
        a=-root['a'];b=-root['b'];al=a-radius;ah=a+radius
        bl=max(F(0),b-radius);bh=b+radius
        assert 0<al<ah<1
        if b==0:
            # A real-centered isolating disk has exactly one root, hence it is real.
            bl=bh=F(0); j=0
        else:
            A=atan_interval(bl/ah,bh/al)
            B=atan_interval(bl/(1-al),bh/(1-ah))
            D=atan_interval(bl/(ah+m),bh/(al+m))
            anglelo=(2*m+2)*A[0]+2*m*B[0]-2*D[1]
            anglehi=(2*m+2)*A[1]+2*m*B[1]-2*D[0]
            # Polynomial equation gives real part of psi'=0 and Im psi' an
            # integral multiple of2pi. A strict interval isolates that integer.
            k=round(float((anglelo+anglehi)/(2*(PI[0]+PI[1]))))
            assert 1<=k<=m-1
            assert anglelo>(2*k-1)*PI[1] and anglehi<(2*k+1)*PI[0]
            j=-k
        assert root['eps']==(-1)**(j+m+1)
        numeratorlo=(1-ah)**2+bl**2; numeratorhi=(1-al)**2+bh**2
        denominatorlo=(al+m)**2+bl**2; denominatorhi=(ah+m)**2+bh**2
        moduluslo=(m*m*numeratorlo/denominatorhi)**m
        modulushi=(m*m*numeratorhi/denominatorlo)**m
        assert 0<moduluslo<modulushi<1
        nonreal=False
        if b:
            B=atan_interval(bl/(1-al),bh/(1-ah))
            D=atan_interval(bl/(ah+m),bh/(al+m))
            # -arg(eta)/pi modulo2 equals a positive explicitly bounded angle.
            phase_lo=2*m*(B[0]+D[0])/PI[1]
            phase_hi=2*m*(B[1]+D[1])/PI[0]
            assert floor(phase_lo)==floor(phase_hi) and phase_lo>floor(phase_lo)
            nonreal=True
        components.append({'j':j,'t_real_center':str(a),'t_imag_center':str(b),'radius':str(radius),
                           'modulus_lower':str(moduluslo),'modulus_upper':str(modulushi),
                           'modulus_decimal_interval':[float(moduluslo),float(modulushi)],
                           'nonreal_value':nonreal,'rouche_margin':str(root['rouche_margin'])})
    assert sorted(row['j'] for row in components)==list(range(1-m,1))
    dominant=next(row for row in components if row['j']==1-m)
    for row in components:
        if row is not dominant:assert F(dominant['modulus_lower'])>F(row['modulus_upper'])
    assert dominant['nonreal_value']
    records.append({'m':m,'components_upper_half':components,'dominant_j':1-m})
    print('PASS m=',m,'all2m roots isolated; all components assigned; unique dominant pair',
          dominant['modulus_decimal_interval'])

out=Path(__file__).resolve().parent/'rodrigues_finite_saddles_certificate.json'
out.write_text(json.dumps({'method':'exact rational Rouche, alternating atan/Machin pi, rational modulus bounds',
                           'records':records},indent=2)+'\n')
print('Saved',out)
