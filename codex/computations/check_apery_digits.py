"""Independent direct-binomial checks of digit arithmetic; standard library only."""
from math import comb


def B(n, a=2, b=2, t=1):
    return sum(comb(n, k)**a * comb(n+k, k)**b * t**k
               for k in range(n+1))


def mul(x, y, p, length):
    out = [0]*length
    for i, u in enumerate(x):
        for j, v in enumerate(y[:max(0, length-i)]):
            out[i+j] = (out[i+j] + u*v) % p
    return out


def main():
    ps = (2, 3, 5, 7, 11, 13, 17, 19)
    N = 400
    A = [B(n) for n in range(N+1)]
    cases = 0
    for p in ps:
        digits = [A[r] % p for r in range(p)]
        print('p =', p, 'digit residues =', digits)
        for n in range(N+1):
            assert A[n] % p == A[n//p]*A[n % p] % p
            q, product = n, 1
            while q:
                product = product * digits[q % p] % p
                q //= p
            assert A[n] % p == product
            cases += 1
        L = 130
        f = [x % p for x in A[:L]]
        # Verify P_p F^{p-1}=1 as truncated series.
        power = [1] + [0]*(L-1)
        for _ in range(p-1):
            power = mul(power, f, p, L)
        assert mul(digits, power, p, L) == [1]+[0]*(L-1)
    for a, b in ((1,1), (1,2), (2,1), (2,2), (3,1), (1,3)):
        for t in (-2, 0, 1, 2, 5):
            values = [B(n,a,b,t) for n in range(100)]
            for p in ps[:5]:
                for n, v in enumerate(values):
                    assert v % p == values[n//p]*values[n%p] % p
                    cases += 1
    polynomial_cases = 0
    for a, b in ((1,1), (1,2), (2,1), (2,2), (3,1)):
        for p in ps[:5]:
            for r in range(p):
                lhs = [comb(p-1-r,k)**a * comb(p-1-r+k,k)**b % p
                       if k <= p-1-r else 0 for k in range(p)]
                rhs = [comb(r,k)**b * comb(r+k,k)**a * (-1)**((a+b)*k) % p
                       if k <= r else 0 for k in range(p)]
                assert lhs == rhs
                polynomial_cases += 1
                for n in range(4):
                    top = p*n+r
                    lhs = [comb(top,k)**a * comb(top+k,k)**b % p
                           for k in range(top+1)]
                    rhs = [0]*(top+1)
                    for j in range(n+1):
                        for k in range(r+1):
                            rhs[p*j+k] = (comb(n,j)**a * comb(n+j,j)**b
                                          * comb(r,k)**a * comb(r+k,k)**b) % p
                    assert lhs == rhs
                    polynomial_cases += 1
    print(f'Passed {polynomial_cases} full polynomial factorization/reflection checks.')
    assert all(x % 2 == 1 for x in A)
    assert all(x % 3 == (-1)**n % 3 for n,x in enumerate(A))
    assert [A[r] % 5 for r in range(5)] == [1,0,3,0,1]
    print(f'Passed {cases} numeric Lucas/digit checks and 8 series identities.')


if __name__ == '__main__':
    main()
