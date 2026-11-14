"""Independent depth-weighted tree recurrence checks; standard library only."""
from math import factorial


def conv(a,b,N):
    out=[0]*(N+1)
    for i,x in enumerate(a):
        if not x:continue
        for j,y in enumerate(b[:N+1-i]):
            if y:out[i+j]+=x*y
    return out


def counts(q,N,b):
    f=[1]+[0]*N
    for h in range(N-1,-1,-1):
        power=[1]+[0]*N
        for _ in range(q):power=conv(power,f,N)
        f=[1]+[b(h)*power[n-1] for n in range(1,N+1)]
    return f


def main():
    checked=0
    for p,e,N in ((2,1,20),(2,2,12),(2,3,8),(3,1,14),(3,2,8),(5,1,8)):
        q=p**e
        for u in range(1,p):
            for sign in (-1,1):
                b=lambda h,u=u,p=p,sign=sign:u+sign*p*(h*h+2*h+1)
                ws=counts(q,N,b)
                for n,w in enumerate(ws):
                    L=(q-1)*n+1;ds=[]
                    while L:ds.append(L%p);L//=p
                    D=(sum(ds)-1)//(p-1)
                    unit=pow(u,n,p)*(-1)**D
                    for d in ds:unit=unit*pow(factorial(d),-1,p)%p
                    assert w%p**D==0
                    assert w//p**D%p==unit
                    checked+=1
    assert counts(2,2,lambda h:1+2*h)[2]==6
    assert 1**2+1*3==4
    print(f'Passed {checked} weighted-depth valuation/leading-unit checks.')
    print('Ordinary-depth vs left-edge example: 6 versus4 at n2, weights1,3.')

if __name__=='__main__':main()
