"""Combinatorial excess-two profile formula and independent orbit/raw-tree checks."""
from collections import Counter
from functools import lru_cache
from itertools import product, combinations, combinations_with_replacement
from math import comb
from check_tree_orbits import compositions, digits
from check_prime_power_sylow import slot_orbit, valuation_power
from check_next_digit_orbits import plane_profile, predicted


def excess_two(p,e,n):
    assert p>=5
    ds=digits((p**e-1)*n+1,p);D=(sum(ds)-1)//(p-1)
    N0,_=plane_profile(p,e,ds)
    M2=0;details=[]
    for j,k in combinations_with_replacement(range(len(ds)-1),2):
        c=[0]*len(ds);c[j]+=1;c[k]+=1
        a=[d+p*c[i]-(c[i-1] if i else 0) for i,d in enumerate(ds)]
        if min(a)<0:continue
        assert sum(a)==sum(ds)+2*(p-1)
        Na,_=plane_profile(p,e,a)
        marked=0
        for h,count in enumerate(a[:-1]):
            if count<p:continue
            b=a[:];b[h]-=p;b[h+1]+=1
            Nb,_=plane_profile(p,e,b)
            marked+=(a[h+1]+1)*Nb
        z=comb(ds[j+1],2) if j==k else ds[j+1]*ds[k+1]
        pairs=z*N0
        F=Na-marked+pairs
        assert F>=0 and F%p**(D+2)==0
        M2+=F//p**(D+2)
        details.append((a,Na,marked,pairs,F))
    return M2,details


def enumerate_orbits(p,e,N):
    q=p**e;by_size=[[0]];exp={0:0};nextid=1
    for n in range(1,N+1):
        orbitkeys={}
        for sizes in compositions(n-1,q):
            for children in product(*(by_size[k] for k in sizes)):
                key,rootsize=slot_orbit(children,p,e)
                d=sum(exp[x] for x in children)+valuation_power(rootsize,p)
                if key in orbitkeys:assert orbitkeys[key]==d
                else:orbitkeys[key]=d
        ids=[];hist=Counter()
        for key,d in orbitkeys.items():ids.append(nextid);exp[nextid]=d;nextid+=1;hist[d]+=1
        by_size.append(ids)
        D,M0,M1,_=predicted(p,e,n);M2,_=excess_two(p,e,n)
        assert hist[D]==M0 and hist.get(D+1,0)==M1 and hist.get(D+2,0)==M2
        assert sum(count*p**d for d,count in hist.items())==comb(q*n,n)//((q-1)*n+1)
        print(f'p={p},e={e},n={n}: M0,M1,M2={M0},{M1},{M2}; complete histogram checked')


@lru_cache(None)
def shapes(p,K):
    if K==0:return (None,)
    out=[]
    for ks in compositions(K-1,p):out.extend(product(*(shapes(p,k) for k in ks)))
    return tuple(out)


def label_shape(shape, labels):
    it=iter(labels)
    def rec(t):return next(it) if t is None else tuple(rec(x) for x in t)
    return rec(shape)


def visible_stars(t):
    if isinstance(t,int):return 0
    return int(all(isinstance(x,int) and x==t[0] for x in t))+sum(visible_stars(x) for x in t)


def raw_enumeration(p,K,ones):
    leaves=(p-1)*K+1;hist=Counter()
    for t in shapes(p,K):
        for positions in combinations(range(leaves),ones):
            labels=[0]*leaves
            for i in positions:labels[i]=1
            hist[visible_stars(label_shape(t,labels))]+=1
    return hist


def main():
    checked=0
    for p in (5,7,11):
        for e in (1,2,3):
            q=p**e
            for n in range(201):
                D,M0,M1,_=predicted(p,e,n)
                M2,_=excess_two(p,e,n)
                C=comb(q*n,n)//((q-1)*n+1)
                assert C//p**D%p**3==(M0+p*M1+p*p*M2)%p**3
                checked+=1
    print(f'Passed {checked} normalized-mod-p3/profile-count checks.')
    assert excess_two(5,1,4)==(1,[([12,1],455,360,30,125)])
    assert excess_two(5,1,6)[0]==25 and excess_two(5,2,1)[0]==0
    h=raw_enumeration(5,3,1)
    assert h==Counter({0:125,1:300,2:30}),h
    nested=raw_enumeration(5,2,4)
    assert nested==Counter({0:625,1:5}),nested
    print('Actual raw labeled skeleton enumeration: 455 and630 cases match all star buckets.')
    for p,e,N in ((5,1,7),(5,2,4),(7,1,5),(7,2,3)):
        enumerate_orbits(p,e,N)
    # Negative boundary: p3,n7 has a repeated-nonperfect excess-two orbit.
    N0=plane_profile(3,1,[0,2,1])[0]
    f0=plane_profile(3,1,[6,0,1])[0]-plane_profile(3,1,[3,1,1])[0]+N0
    f1=plane_profile(3,1,[3,4,0])[0]-5*plane_profile(3,1,[0,5,0])[0]-plane_profile(3,1,[3,1,1])[0]+2*N0
    assert (f0+f1)//27==14
    actual=(comb(21,7)//15)//3%27
    two_uncarry_prediction=(1+3*6+9*14)%27
    assert actual==19 and two_uncarry_prediction==10 # One omitted orbit contributes9.
    print('All exact next-next orbit counts and modulo-p3 checks passed.')

if __name__=='__main__':main()
