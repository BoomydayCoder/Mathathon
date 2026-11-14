"""Full visible-star IE recipe for all k<p-1; exact independent comparisons."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb
from check_tree_orbits import compositions, digits
from check_prime_power_sylow import slot_orbit, valuation_power
from check_next_digit_orbits import plane_profile, predicted
from check_excess_two_orbits import excess_two


def carry_vectors(k,length):
    if length==0:
        if k==0:yield ()
        return
    yield from compositions(k,length)


def orbit_counts(p,e,n,maxk=None):
    ds=digits((p**e-1)*n+1,p);h=len(ds);D=(sum(ds)-1)//(p-1)
    if maxk is None:maxk=p-2
    @lru_cache(None)
    def N(a):return plane_profile(p,e,list(a))[0]
    counts=[]
    for k in range(maxk+1):
        total=0
        for shortc in carry_vectors(k,h-1):
            c=shortc+(0,)
            a=tuple(ds[i]+p*c[i]-(c[i-1] if i else 0) for i in range(h))
            if min(a)<0:continue
            assert sum(a)==sum(ds)+(p-1)*k
            bounds=[x//p for x in a]
            assert all(bounds[i]<=c[i] for i in range(h)) and bounds[-1]==0
            F=0
            for b in product(*(range(r+1) for r in bounds)):
                ap=tuple(a[i]-p*b[i]+(b[i-1] if i else 0) for i in range(h))
                assert min(ap)>=0
                multiplicity=1
                for j in range(h-1):
                    assert ap[j+1]>=b[j]
                    multiplicity*=comb(ap[j+1],b[j])
                F+=(-1)**sum(b)*N(ap)*multiplicity
            assert F>=0 and F%p**(D+k)==0
            total+=F//p**(D+k)
        counts.append(total)
    return D,counts


def enumerate_sylow(p,e,N):
    q=p**e;by_size=[[0]];exp={0:0};nextid=1
    for n in range(1,N+1):
        keys={}
        for sizes in compositions(n-1,q):
            for children in product(*(by_size[k] for k in sizes)):
                key,size=slot_orbit(children,p,e)
                d=sum(exp[x] for x in children)+valuation_power(size,p)
                if key in keys:assert keys[key]==d
                else:keys[key]=d
        ids=[];hist=Counter()
        for key,d in keys.items():ids.append(nextid);exp[nextid]=d;nextid+=1;hist[d]+=1
        by_size.append(ids)
        D,Ms=orbit_counts(p,e,n)
        assert all(hist.get(D+k,0)==M for k,M in enumerate(Ms))
        assert sum(count*p**d for d,count in hist.items())==comb(q*n,n)//((q-1)*n+1)
        print(f'p={p},e={e},n={n}: all {len(Ms)} claimed low-excess orbit buckets match')


def main():
    configs=((2,1,100),(2,2,100),(3,1,100),(3,2,100),
             (5,1,60),(5,2,60),(5,3,40),(7,1,25),(7,2,25),
             (11,1,10),(11,2,10))
    checked=0
    for p,e,N in configs:
        q=p**e;mod=p**(p-1)
        for n in range(N+1):
            D,Ms=orbit_counts(p,e,n)
            C=comb(q*n,n)//((q-1)*n+1)
            assert C//p**D%mod==sum(p**k*M for k,M in enumerate(Ms))%mod
            if p>=3:
                oldD,M0,M1,_=predicted(p,e,n)
                assert oldD==D and Ms[:2]==[M0,M1]
            if p>=5:assert Ms[2]==excess_two(p,e,n)[0]
            checked+=1
    print(f'Passed {checked} complete-low-excess/mod-p^(p-1) checks and exact k1/k2 comparisons.')
    for p,e,N in ((2,1,8),(3,1,7),(5,1,7),(5,2,4),(7,1,5),(7,2,3)):
        enumerate_sylow(p,e,N)
    # The first excluded repeated-block orbit produces an actual wrong next digit.
    for p in (3,5,7):
        n=2*p+1;q=p
        D,starfree=orbit_counts(p,1,n,maxk=p-1)
        assert D==1
        actual=(comb(q*n,n)//((q-1)*n+1))//p**D
        approximation=sum(p**k*M for k,M in enumerate(starfree))
        assert (actual-approximation)%p**p==p**(p-1)
        print(f'p={p},q=p,n={n}: extending naively to k=p-1 misses exactly one boundary orbit modulo next digit')
    print('All general IE, exact orbit buckets, and sharp-boundary checks passed.')

if __name__=='__main__':main()
