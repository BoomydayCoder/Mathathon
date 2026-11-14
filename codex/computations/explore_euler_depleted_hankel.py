"""Exact small determinants for original and Euler-depleted positive models.

All polynomials use the SAME variable X=zeta(s), substituting
eta(s)=(1-2^(1-s))*X in the depleted models. Numerical evaluations are
diagnostics, not certified asymptotic bounds.
"""
from pathlib import Path
import json
import sympy as S
import mpmath as mp

t, X = S.symbols('t X')
mp.mp.dps = 250


def check_model():
    polynomial_checks = 0
    pole_checks = 0
    for s in range(3, 20, 2):
        h = s-1
        a = (h & -h).bit_length()-1
        uniform = 2**(a+1)
        for e in range(81):
            k = 2*e+2
            G = 2*(1-2**k)*S.bernoulli(k)
            assert G.is_Integer and int(G) % 2
            value = (-1)**(e+1)*G*S.binomial(k+s-2, s-2)/(2*h)
            assert (uniform*value).is_Integer and value > 0
            polynomial_checks += 1
        e_sharp = 2**(a-1)-1
        k = 2*e_sharp+2
        value = ((-1)**e_sharp*(2**k-1)*S.bernoulli(k)
                 *S.factorial(2*e_sharp+s)/(S.factorial(s-1)*S.factorial(k)))
        assert S.denom(value) == uniform
        # Compute the scaled difference using Hurwitz tails independently.
        for j in range(1, 41):
            if j % 2 == 0:
                hz_lead = S.Integer(1)
                hz_const = -sum(S.Rational(1, k**s) for k in range(1, j//2))
            else:
                hz_lead = S.Integer(2**s-1)
                hz_const = -2**s*sum(S.Rational(1, (2*k-1)**s)
                                      for k in range(1, (j+1)//2))
            left_a = j**(s-1)*(S.Rational(1, 2**(s-1))*hz_lead-1)
            left_b = j**(s-1)*(S.Rational(1, 2**(s-1))*hz_const
                       +sum(S.Rational(1, k**s) for k in range(1, j)))
            left_b -= S.Rational(1, 2*j)
            right_a = (-1)**(j-1)*j**(s-1)*(1-S.Rational(1, 2**(s-1)))
            right_b = -(-1)**(j-1)*j**(s-1)*sum(
                S.Rational((-1)**(k-1), k**s) for k in range(1, j+1))
            right_b += S.Rational(1, 2*j)
            assert left_a == right_a and left_b == right_b
            pole_checks += 1
    print(f'PASS: {polynomial_checks} dyadic moments, sharp constants, '
          f'{pole_checks} exact Hurwitz pole identities.', flush=True)


def parts(s, K, N, family):
    nodes = list(range(1, K+1)) if family != 'depleted_odd' else list(range(1, 2*K, 2))
    dn = S.Poly(S.prod(t+j*j for j in nodes[:N]), t)
    dk = S.Poly(S.prod(t+j*j for j in nodes), t)
    numerator = dn**6
    h = K-N
    out = []
    for d in range(2*h-1):
        q = S.quo(numerator*S.Poly(t**d, t), dk)
        lead = S.Rational(0)
        const = S.Rational(0)
        for (e,), coeff in q.terms():
            if not coeff:
                continue
            moment = ((-1)**e*S.bernoulli(2*e+2)
                      *S.factorial(2*e+s)/(S.factorial(s-1)*S.factorial(2*e+2)))
            if family != 'original':
                moment *= 2**(2*e+2)-1
            const += coeff*moment
        for j in nodes[N:]:
            residue = S.Rational(numerator.eval(-j*j)*(-j*j)**d,
                                 S.prod(v*v-j*j for v in nodes if v != j))
            if family == 'original':
                a = j**(s-1)
                b = -a*sum(S.Rational(1, k**s) for k in range(1, j+1))
                b += -S.Rational(1, s-1)+S.Rational(1, 2*j)
            else:
                sign = (-1)**(j-1)
                a_eta = sign*j**(s-1)
                a = a_eta*(1-S.Rational(1, 2**(s-1)))
                b = -a_eta*sum(S.Rational((-1)**(k-1), k**s)
                               for k in range(1, j+1))+S.Rational(1, 2*j)
            lead += residue*a
            const += residue*b
        out.append((lead, const))
    return out


def val2(q):
    numerator, denominator = map(int, S.fraction(q))
    def power(n):
        n = abs(n)
        return (n & -n).bit_length()-1
    return power(numerator)-power(denominator)


def run(s, K, N, family):
    h = K-N
    values = parts(s, K, N, family)
    matrix = S.Matrix(h, h, lambda i, j: values[i+j][0]*X+values[i+j][1])
    poly = S.Poly(matrix.det(method='domain-ge'), X)
    denominator, integer = poly.clear_denoms(convert=True)
    content, primitive = integer.primitive()
    scale = S.Rational(denominator, content)
    coeffs = [int(c) for c in primitive.all_coeffs()]
    evaluation = mp.polyval(coeffs, mp.zeta(s))
    assert poly.degree() == h and evaluation > 0
    result = dict(family=family, s=s, K=K, N=N, degree=h,
                  coefficient_denominator=str(denominator),
                  primitive_multiplier=str(scale), multiplier_v2=val2(scale),
                  coefficients=coeffs, log_primitive_value=str(mp.log(evaluation)),
                  coefficient_denominator_factors={str(p): int(e) for p, e in
                                                  S.factorint(denominator).items()})
    print(f'{family:14} s={s} K={K} N={N} v2(scale)={val2(scale):4} '
          f'log primitive value={float(mp.log(evaluation)): .6f}', flush=True)
    return result


if __name__ == '__main__':
    check_model()
    rows = []
    for K in range(1, 7):
        for N in [0]+([1] if K > 1 else []):
            for s in (5, 7):
                for family in ('original', 'depleted_all', 'depleted_odd'):
                    rows.append(run(s, K, N, family))
    Path(__file__).with_suffix('.json').write_text(
        json.dumps(dict(precision=mp.mp.dps, target_variable='zeta(s)', rows=rows),
                   indent=2)+'\n')
