"""Odd-prime second-unit-digit checks and independent next-orbit enumeration."""
from collections import Counter
from itertools import product
from math import comb, factorial
from check_tree_orbits import compositions, digits
from check_prime_power_sylow import slot_orbit, valuation_power


def plane_profile(p,e,a):
    ls=[sum(x for j,x in enumerate(a) if (i+j)%e==0) for i in range(e)]
    q=p**e
    ks=[]
    for i in range(e):
        num=sum(p**r*ls[(i-r)%e] for r in range(e))-p**i
        assert num%(q-1)==0
        ks.append(num//(q-1))
    assert all(k>=0 for k in ks)
    num=1;den=1
    for k in ks:num*=factorial(p*k);den*=factorial(k)
    for x in a:den*=factorial(x)
    assert num%den==0
    return num//den,ks


def predicted(p,e,n):
    L=(p**e-1)*n+1;ds=digits(L,p);D=(sum(ds)-1)//(p-1)
    N0,ks=plane_profile(p,e,ds)
    assert N0%p**D==0
    M=N0//p**D
    M1=0
    H=[0]
    for r in range(1,p):H.append((H[-1]+pow(r,-1,p))%p)
    S=sum(ds[j+1]*H[ds[j]] for j in range(len(ds)-1))%p
    for j in range(len(ds)-1):
        c=ds[j+1]
        if not c:continue
        a=ds[:];a[j]+=p;a[j+1]-=1
        Na,ka=plane_profile(p,e,a)
        t=(-j-1)%e
        assert ka==[k+int(i==t) for i,k in enumerate(ks)]
        difference=Na-c*N0
        assert difference>=0 and difference%p**(D+1)==0
        M1+=difference//p**(D+1)
    assert M1%p==(-M*S)%p
    mod=p*p;den=1
    for d in ds:den=den*factorial(d)%mod
    rhs=pow(factorial(p-1),D,mod)*pow(den,-1,mod)*(1-p*S)%mod
    assert M%mod==pow(factorial(p-1),D,mod)*pow(den,-1,mod)%mod
    assert (M+p*M1)%mod==rhs
    return D,M,M1,rhs


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
        D,M,M1,rhs=predicted(p,e,n)
        assert hist[D]==M and hist.get(D+1,0)==M1
        assert sum(count*p**d for d,count in hist.items())==comb(q*n,n)//((q-1)*n+1)
        print(f'p={p},e={e},n={n}: minimum orbits={M}, next orbits={M1}, exact hist verified')


def main():
    checked=0
    for p in (3,5,7):
        for e in (1,2,3):
            q=p**e
            for n in range(501):
                D,M,M1,rhs=predicted(p,e,n)
                C=comb(q*n,n)//((q-1)*n+1)
                assert C//p**D%(p*p)==rhs
                checked+=1
    print(f'Passed {checked} second-unit-digit, profile-count and next-orbit residue checks.')
    for p,e,N in ((3,1,9),(3,2,5),(5,1,6),(5,2,3)):
        enumerate_orbits(p,e,N)
    assert (comb(4,2)//3)//2==1 # Actual normalized Catalan value at p2,n2.
    assert (1-2*1)%4==3 # Blind odd-prime formula would instead predict3.
    print('All independent next-size orbit histograms agree with the uncarry-minus-star formula.')

if __name__=='__main__':main()
