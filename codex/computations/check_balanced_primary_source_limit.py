#!/usr/bin/env python3
"""Finite raising-operator check of the Ismail--Rains--Stanton degeneration."""
from pathlib import Path
from functools import lru_cache
import json
import mpmath as mp
base=Path(__file__).parent
scope={}
exec((base/'check_balanced_rodrigues_integrality.py').read_text().split('rows=[];conditions=0')[0],scope)
mp.mp.dps=90
rows=[]
for m,n in [(1,1),(1,3),(2,1),(2,3),(3,2),(3,4),(4,3)]:
 cf=scope['integral_rodrigues'](m,n)
 x=mp.mpf('1.3')
 target=mp.factorial(n)*sum(c*x**(2*e) for e,c in enumerate(cf))
 errors=[]
 for eps in ['0.01','0.001','0.0001']:
  q=1-mp.mpf(eps)
  @lru_cache(None)
  def inside(s,h):
   if s==n:return mp.mpf(1)
   xx=x+mp.mpf(h)/2
   cs=lambda z:-q**(-mp.mpf('.5')-2*z)*(1-q**(mp.mpf(s)/2+z))*(1-q**(1+mp.mpf(s)/2+z))**(2*m+1)
   top=cs(xx)*inside(s+1,h+1)-cs(-xx)*inside(s+1,h-1)
   return 2*top/((q**mp.mpf('.5')-q**(-mp.mpf('.5')))*(q**xx-q**(-xx)))
  value=inside(0,0)/((-2)**n*(1-q)**(2*m*n))
  errors.append(abs(value/target-1))
 assert errors[2]<errors[1]/5 and errors[1]<errors[0]/5
 assert errors[2]<mp.mpf('.01')
 rows.append(dict(m=m,n=n,relative_errors=[mp.nstr(e,20) for e in errors]))
 print('PASS finite source degeneration:',m,n,[mp.nstr(e,6) for e in errors])
(base/'balanced_primary_source_limit_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
