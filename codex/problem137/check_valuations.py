"""Exact finite verification of the independently proved forced-carry formula.
Run: python codex/problem137/check_valuations.py
No external dependencies. Finite checks are not a substitute for the proof.
"""
from math import comb


def vp(x, p):
    assert x > 0
    out = 0
    while x % p == 0:
        out += 1
        x //= p
    return out


def carries(x, y, p):
    out, incoming, pos = [], 0, 0
    while x or y or incoming:
        incoming = (x % p + y % p + incoming) // p
        if incoming:
            out.append(pos)
        x //= p
        y //= p
        pos += 1
    return out


def transducer(m, n, p):
    u, c, t, score = 0, 0, True, 0
    while n or u or c or t:
        d, n = n % p, n // p
        e, u = (m*d + u) % p, (m*d + u) // p
        c = (d + e + c) // p
        t = t and e == p - 1
        assert c - int(t) >= 0
        score += c - int(t)
    return score


def main():
    checked = 0
    for m in range(1, 21):
        for n in range(301):
            b = comb((m+1)*n, n)
            C, rem = divmod(b, m*n+1)
            assert rem == 0
            for p in (2, 3, 5, 7, 11, 13):
                a, cs = vp(m*n+1, p), carries(n, m*n, p)
                assert all(j in cs for j in range(a))
                value = vp(C, p)
                assert value == sum(j >= a for j in cs)
                assert value == transducer(m, n, p)
                checked += 1
    prime_power_checks = 0
    for p in (2, 3, 5, 7):
        for e in (1, 2, 3):
            q = p**e
            repunits = {0}
            r = 0
            while r <= 150:
                repunits.add(r)
                r = q*r + 1
            for n in range(151):
                C = comb(q*n,n)//((q-1)*n+1)
                assert (vp(C,p) == 0) == (n in repunits)
                x, ds = (q-1)*n, 0
                while x:
                    ds += x % p
                    x //= p
                assert vp(C,p) == ds//(p-1)-vp((q-1)*n+1,p)
                prime_power_checks += 1
    print(f'Passed {prime_power_checks} prime-power-arity classification checks.')
    print(f'Passed {checked} exact valuation/carry/transducer checks.')


if __name__ == '__main__':
    main()
