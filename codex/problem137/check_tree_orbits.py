"""Enumerate child-orbit necklaces and check p-group valuation statements.
Full trees are not expanded. Orbit types are recursive canonical child tuples.
"""
from collections import Counter
from itertools import product
from math import comb, factorial


def compositions(total, slots):
    if slots==1:
        yield (total,)
    else:
        for a in range(total+1):
            for rest in compositions(total-a,slots-1):yield (a,)+rest


def vp(x,p):
    d=0
    while x%p==0:d+=1;x//=p
    return d


def digits(x,p):
    out=[]
    while x:out.append(x%p);x//=p
    return out or [0]


def cyclic_images(t):
    return {t[k:]+t[:k] for k in range(len(t))}


def sylow4_images(t):
    out=set()
    for swap,flip0,flip1 in product(range(2),repeat=3):
        blocks=[t[:2],t[2:]]
        if flip0:blocks[0]=blocks[0][::-1]
        if flip1:blocks[1]=blocks[1][::-1]
        if swap:blocks.reverse()
        out.add(blocks[0]+blocks[1])
    return out


def enumerate_orbits(q,p,N,images):
    by_size=[[0]]
    exponent={0:0}
    next_id=1
    histograms=[Counter({0:1})]
    for n in range(1,N+1):
        canonical={}
        for sizes in compositions(n-1,q):
            for children in product(*(by_size[k] for k in sizes)):
                imgs=images(children)
                canonical.setdefault(min(imgs),sum(exponent[v] for v in children)+vp(len(imgs),p))
        ids=[];hist=Counter()
        for childtuple,d in canonical.items():
            ids.append(next_id);exponent[next_id]=d;next_id+=1;hist[d]+=1
        by_size.append(ids);histograms.append(hist)
        total=sum(p**d*count for d,count in hist.items())
        assert total==comb(q*n,n)//((q-1)*n+1)
        fixed=hist.get(0,0)
        r=0
        while r<n:r=q*r+1
        assert fixed==int(r==n)
        if q==p:
            ds=digits((p-1)*n+1,p);D=(sum(ds)-1)//(p-1)
            assert min(hist)==D
            M=factorial(p*D)//(p**D*factorial(D))
            for d in ds:M//=factorial(d)
            # Use exact product denominator as a separate check.
            den=p**D*factorial(D)
            for d in ds:den*=factorial(d)
            M=factorial(p*D)//den
            assert hist[D]==M and M%p!=0
        print(f'q={q}, n={n}, trees={total}, min exponent={min(hist)}, min orbit count={hist[min(hist)]}, valuation={vp(total,p)}')
    return histograms


def main():
    enumerate_orbits(2,2,12,cyclic_images)
    enumerate_orbits(3,3,7,cyclic_images)
    enumerate_orbits(4,2,8,cyclic_images)
    enumerate_orbits(4,2,8,sylow4_images)
    # Large cyclic-q4 counterexample can be counted without full enumeration.
    q,n,p=4,28,2
    value=comb(q*n,n)//((q-1)*n+1)
    assert vp(value,p)==3
    perfect_sizes=[0,1,5,21]
    tuples=[t for t in product(perfect_sizes,repeat=4) if sum(t)==27]
    orbits={min(cyclic_images(t)) for t in tuples}
    assert len(tuples)==24 and len(orbits)==6
    assert all(len(cyclic_images(t))==4 for t in tuples)
    sylow_orbits={min(sylow4_images(t)) for t in tuples}
    assert len(sylow_orbits)==3 and all(len(sylow4_images(t))==8 for t in tuples)
    print('q4 n28: cyclic minimum = 6 orbits of size4; Sylow gives3 orbits of size8; exact valuation3.')
    print('All recursive orbit enumerations and prime-arity minimum-orbit formulas passed.')

if __name__=='__main__':main()
