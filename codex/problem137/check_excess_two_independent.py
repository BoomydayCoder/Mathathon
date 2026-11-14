"""Independent plane-tree/star/orbit audit; imports no project helper code."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import factorial,comb

def compositions(n,k):
    if k==1:
        yield (n,)
        return
    for a in range(n+1):
        for rest in compositions(n-a,k-1):
            yield (a,)+rest

@lru_cache(None)
def trees(p,n):
    if n==0:
        return (None,)
    out=[]
    for sizes in compositions(n-1,p):
        out.extend(product(*(trees(p,k) for k in sizes)))
    return tuple(out)

def profile_count(p,e,a):
    leaves=[sum(a[j] for j in range(len(a)) if (j+t)%e==0) for t in range(e)]
    ks=[]
    for t in range(e):
        numerator=sum(p**r*leaves[(t-r)%e] for r in range(e))-p**t
        assert numerator%(p**e-1)==0
        ks.append(numerator//(p**e-1))
    assert min(ks)>=0
    numerator=1
    denominator=1
    for k in ks:
        numerator*=factorial(p*k)
        denominator*=factorial(k)
    for count in a:
        denominator*=factorial(count)
    assert numerator%denominator==0
    return numerator//denominator

def normalize(p,a):
    digits=[]
    carries=[]
    carry=0
    for count in a:
        digit=(count+carry)%p
        carry=(count+carry)//p
        digits.append(digit)
        carries.append(carry)
    while carry:
        digits.append(carry%p)
        carry//=p
        carries.append(carry)
    return digits,carries

def describe_shape(tree):
    depths=[]
    stars=[]
    def rec(t,depth):
        if t is None:
            depths.append(depth)
            return
        start=len(depths)
        for child in t:
            rec(child,depth+1)
        if all(child is None for child in t):
            stars.append(tuple(range(start,len(depths))))
    rec(tree,0)
    return depths,stars

def assignments(a,depths,e):
    counts=list(a)
    labels=[]
    def rec(i):
        if i==len(depths):
            yield tuple(labels)
            return
        for label,count in enumerate(counts):
            if count and (label+depths[i])%e==0:
                counts[label]-=1
                labels.append(label)
                yield from rec(i+1)
                labels.pop()
                counts[label]+=1
    yield from rec(0)

def check_raw(p,e,a):
    K=(sum(a)-1)//(p-1)
    digits,carries=normalize(p,a)
    assert sum(carries)==2
    a=list(a)+[0]*(len(digits)-len(a))
    histogram=Counter()
    marked=Counter()
    pairs=Counter()
    for tree in trees(p,K):
        depths,positions=describe_shape(tree)
        for labels in assignments(a,depths,e):
            visible=[labels[pos[0]] for pos in positions
                     if len({labels[i] for i in pos})==1]
            assert len(visible)<=2
            histogram[len(visible)]+=1
            marked.update(visible)
            if len(visible)==2:
                pairs[tuple(sorted(visible))]+=1
    assert sum(histogram.values())==profile_count(p,e,a)
    predicted_marked=Counter()
    for j,count in enumerate(a[:-1]):
        if count>=p:
            b=a[:]
            b[j]-=p
            b[j+1]+=1
            predicted_marked[j]=(a[j+1]+1)*profile_count(p,e,b)
    assert marked==predicted_marked
    indices=[j for j,c in enumerate(carries) for _ in range(c)]
    j,k=indices
    z=comb(digits[j+1],2) if j==k else digits[j+1]*digits[k+1]
    expected_pairs=z*profile_count(p,e,digits)
    assert histogram[2]==expected_pairs
    assert pairs==Counter({(j,k):expected_pairs}) if expected_pairs else not pairs
    free=profile_count(p,e,a)-sum(predicted_marked.values())+expected_pairs
    assert free==histogram[0]
    D=(sum(digits)-1)//(p-1)
    assert free%p**(D+2)==0
    print('PASS raw',p,e,a,'buckets',dict(sorted(histogram.items())))

@lru_cache(None)
def orbit(tree):
    if tree is None:
        return (),0
    children=[orbit(t) for t in tree]
    keys=tuple(k for k,d in children)
    exponent=sum(d for k,d in children)+(0 if len(set(keys))==1 else 1)
    rotations=[keys[i:]+keys[:i] for i in range(len(keys))]
    return min(rotations),exponent

def canonical_profile(tree,p):
    def rec(t):
        if t is None:
            return 0,Counter({0:1})
        children=[rec(x) for x in t]
        heights=[h for h,c in children]
        if all(h is not None and h==heights[0] for h in heights):
            h=heights[0]+1
            return h,Counter({h:1})
        out=Counter()
        for h,c in children:
            out.update(c)
        return None,out
    return rec(tree)[1]

def check_ternary_counterexample():
    p=3
    leaf=None
    perfect=(leaf,leaf,leaf)
    nonperfect=(perfect,leaf,leaf)
    exceptional=(nonperfect,nonperfect,nonperfect)
    key,exponent=orbit(exceptional)
    assert exponent==3 and canonical_profile(exceptional,p)==Counter({0:6,1:3})
    groups=Counter()
    exponents={}
    for t in trees(3,7):
        k,d=orbit(t)
        groups[k]+=1
        exponents[k]=d
    assert len(trees(3,7))==7752
    assert all(groups[k]==3**d for k,d in exponents.items())
    hist=Counter(exponents.values())
    assert hist[1]==1 and hist[2]==6 and hist[3]==15
    assert groups[key]==27
    # Exactly one excess-two orbit has the three-carry canonical profile.
    exceptional_keys=set()
    for t in trees(3,7):
        k,d=orbit(t)
        if d==3 and sum(canonical_profile(t,3).values())==9:
            exceptional_keys.add(k)
    assert exceptional_keys=={key}
    print('PASS ternary full enumeration: 7752 trees, M0/M1/M2=1/6/15; omitted orbit size27.')

def main():
    # Same-level stars, adjacent uncarries (with and without a pair),
    # separated levels, genuinely cyclic depth restrictions, and root cases.
    cases=[
        (5,1,[12,1]),
        (5,1,[5,4,0]),
        (5,2,[5,4,0]),
        (5,1,[8,5,0]),
        (5,1,[8,0,5,0]),
        (5,1,[11,0,2]),
        (5,2,[5,1,5,2]),
        (7,1,[7,6,0]),
        (7,2,[7,6,0]),
    ]
    for case in cases:
        check_raw(*case)
    check_ternary_counterexample()
    print('All independent excess-two checks passed.')

if __name__=='__main__':
    main()
