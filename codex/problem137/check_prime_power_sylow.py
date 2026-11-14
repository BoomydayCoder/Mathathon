"""Full Sylow p-group orbit checks and leading-unit checks; standard library."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb, factorial
from check_tree_orbits import compositions, digits


def parameters(p,e,n):
    q=p**e;L=(q-1)*n+1;ds=digits(L,p);D=(sum(ds)-1)//(p-1)
    ls=[sum(d for j,d in enumerate(ds) if (i+j)%e==0) for i in range(e)]
    ks=[]
    for i in range(e):
        num=sum(p**r*ls[(i-r)%e] for r in range(e))-p**i
        assert num%(q-1)==0
        ks.append(num//(q-1))
    assert all(k>=0 for k in ks) and sum(ks)==D
    assert all(ks[i]+ls[i]==int(i==0)+p*ks[(i-1)%e] for i in range(e))
    top=1;den=1
    for k in ks:top*=factorial(p*k);den*=factorial(k)
    for d in ds:den*=factorial(d)
    assert top%den==0
    N=top//den
    assert N%p**D==0
    M=N//p**D
    unit=(-1)**D
    for d in ds:unit=unit*pow(factorial(d),-1,p)%p
    assert M%p==unit
    return D,M,unit


@lru_cache(None)
def slot_orbit(t,p,e):
    if e==0:return t,1
    b=len(t)//p
    sub=[slot_orbit(t[i*b:(i+1)*b],p,e-1) for i in range(p)]
    keys=tuple(x[0] for x in sub)
    rotations={keys[k:]+keys[:k] for k in range(p)}
    canonical=sum(min(rotations),())
    size=len(rotations)
    for _,s in sub:size*=s
    return canonical,size


def valuation_power(x,p):
    d=0
    while x%p==0:d+=1;x//=p
    assert x==1
    return d


def enumerate_sylow(p,e,N):
    q=p**e
    by_size=[[0]];exponents={0:0};next_id=1
    for n in range(1,N+1):
        canonical={}
        for sizes in compositions(n-1,q):
            for children in product(*(by_size[k] for k in sizes)):
                key,root_size=slot_orbit(children,p,e)
                exp=sum(exponents[v] for v in children)+valuation_power(root_size,p)
                if key in canonical:assert canonical[key]==exp
                else:canonical[key]=exp
        ids=[];hist=Counter()
        for key,exp in canonical.items():
            ids.append(next_id);exponents[next_id]=exp;next_id+=1;hist[exp]+=1
        by_size.append(ids)
        total=sum(count*p**exp for exp,count in hist.items())
        assert total==comb(q*n,n)//((q-1)*n+1)
        D,M,unit=parameters(p,e,n)
        assert min(hist)==D and hist[D]==M
        print(f'p={p}, e={e}, n={n}: D={D}, minimum orbit count={M}, all {len(ids)} orbit types verified')


def main():
    checks=0
    for p in (2,3,5,7):
        for e in (1,2,3,4):
            q=p**e
            for n in range(201):
                D,M,unit=parameters(p,e,n)
                C=comb(q*n,n)//((q-1)*n+1)
                assert C%p**D==0 and C//p**D%p==unit
                checks+=1
    print(f'Passed {checks} exact prime-power-arity valuation and leading-unit checks.')
    for p,e,N in ((2,1,10),(3,1,6),(2,2,8),(2,3,5),(3,2,4)):
        enumerate_sylow(p,e,N)
    print('All full Sylow orbit counts match the typed-skeleton factorial formula.')

if __name__=='__main__':main()
