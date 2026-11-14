"""Independent high-order star and boundary tests.

Counts visible stars by an ordered-child convolution, not by the proposed
inclusion-exclusion formula, source helper functions, or tree enumeration.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb,factorial

@lru_cache(None)
def star_distribution(p,K,ones):
    """Plane p-trees with K internals, zero or one leaf of label1,
    all other leaves label0; polynomial coefficients count visible stars."""
    assert ones in (0,1)
    if K==0:
        return ((0,1),)
    states={(0,0,0):1}
    for child in range(p):
        nxt=Counter()
        for (used,marked,stars),count in states.items():
            for n in range(K-used):
                for one in range(ones-marked+1):
                    for childstars,childcount in star_distribution(p,n,one):
                        nxt[(used+n,marked+one,stars+childstars)]+=count*childcount
        states=nxt
    result=Counter()
    rootstar=int(K==1 and ones==0)
    for (used,marked,stars),count in states.items():
        if used==K-1 and marked==ones:
            result[stars+rootstar]+=count
    return tuple(sorted(result.items()))

def check_high_orders():
    profiles=0
    selections=0
    for p in (2,3,5,7,11):
        for k in range(p-1):
            K=k+1
            a0=(p-1)*K
            histogram=dict(star_distribution(p,K,1))
            assert sum(histogram.values())==factorial(p*K)//(factorial(K)*factorial(a0))
            # Digit profile (p-1-k,k+1), carry c0=k, D=1.
            assert a0//p==k
            alternating=0
            for b in range(k+1):
                actual=sum(count*comb(stars,b) for stars,count in histogram.items() if stars>=b)
                ap0=a0-p*b
                ap1=1+b
                Kp=K-b
                raw=factorial(p*Kp)//(factorial(Kp)*factorial(ap0)*factorial(ap1))
                predicted=comb(ap1,b)*raw
                assert actual==predicted
                alternating+=(-1)**b*predicted
                selections+=1
            assert alternating==histogram.get(0,0)
            assert alternating%p**(1+k)==0
            profiles+=1
    print('PASS convolution:',profiles,'profiles;',selections,'marked subset orders, through order9.')

def visible(tree,path=()):
    if isinstance(tree,int):
        return []
    out=[]
    if all(isinstance(x,int) and x==tree[0] for x in tree):
        out.append((path,tree[0]))
    for i,child in enumerate(tree):
        out.extend(visible(child,path+(i,)))
    return out

def collapse(tree,selected,path=()):
    if path in selected:
        return selected[path]+1
    if isinstance(tree,int):
        return tree
    return tuple(collapse(child,selected,path+(i,)) for i,child in enumerate(tree))

def expand(tree,selected,p,path=()):
    if path in selected:
        assert tree==selected[path]+1
        return (selected[path],)*p
    if isinstance(tree,int):
        return tree
    return tuple(expand(child,selected,p,path+(i,)) for i,child in enumerate(tree))

def profile(tree):
    if isinstance(tree,int):
        return Counter({tree:1})
    out=Counter()
    for child in tree:
        out.update(profile(child))
    return out

def depth_valid(tree,e,depth=0):
    if isinstance(tree,int):
        return (tree+depth)%e==0
    return all(depth_valid(child,e,depth+1) for child in tree)

def check_simultaneous_maps():
    p=5
    examples=[
        (1,((0,)*p,(0,)*p,(1,)*p,2,2)), # Adjacent marked heights, three stars.
        (2,((0,)*p,(0,)*p,(2,)*p,1,1)), # Distinct heights with depth restrictions.
        (3,(((0,)*p,1,1,1,1),2,2,2,2)), # Collapse creates a new ancestor star.
    ]
    subsets=0
    for e,tree in examples:
        assert depth_valid(tree,e)
        stars=visible(tree)
        original=profile(tree)
        for size in range(len(stars)+1):
            for subset in combinations(stars,size):
                selected=dict(subset)
                reduced=collapse(tree,selected)
                assert depth_valid(reduced,e)
                assert expand(reduced,selected,p)==tree
                counts=Counter(label for path,label in subset)
                expected=original.copy()
                for label,count in counts.items():
                    expected[label]-=p*count
                    expected[label+1]+=count
                assert +expected==profile(reduced)
                assert all(profile(reduced)[label+1]>=count for label,count in counts.items())
                subsets+=1
    print('PASS simultaneous collapse/expand:',subsets,'subsets; adjacent, separated and nested-merger cases.')

def normalize(p,counts):
    length=max(counts)+2
    carry=0
    digits=[]
    carries=[]
    for j in range(length):
        value=counts[j]+carry
        digits.append(value%p)
        carry=value//p
        carries.append(carry)
    assert carry==0
    return digits,carries

def check_boundary_every_q():
    checked=0
    for p in (2,3,5,7,11):
        for e in range(1,6):
            if e==1:
                low,high=0,1
                expected_n=2*p+1
            else:
                low,high=e-2,2*e-2
                expected_n=p**(e-1)+1
            # A mixed root with p-1 perfect low branches and one perfect
            # high branch has orbit exponent1 and digit bound1.
            child_leaves=(p-1)*p**low+p**high
            ds=[]
            value=child_leaves
            while value:
                ds.append(value%p);value//=p
            assert (sum(ds)-1)//(p-1)==1
            child_exponent=1
            # The outer root repeats this nonperfect orbit type p times.
            exponent=p*child_exponent
            labels=Counter({low:p*(p-1),high:p})
            digits,carries=normalize(p,labels)
            D=(sum(digits)-1)//(p-1)
            assert D==1 and exponent-D==p-1
            assert sum(carries)==p
            assert (sum(labels.values())-1)//(p-1)==p+1
            L=p*child_leaves
            assert (L-1)%(p**e-1)==0
            assert (L-1)//(p**e-1)==expected_n
            depths={2+low,2+high}
            assert all(depth%e==0 for depth in depths)
            if e>=2:
                assert depths=={e,2*e}
            checked+=1
    print('PASS sharp repeated-block boundary:',checked,'prime-power arities, including p2.')

def main():
    check_high_orders()
    check_simultaneous_maps()
    check_boundary_every_q()
    print('All independent general low-excess checks passed.')

if __name__=='__main__':
    main()
