"""Exact algebra checks for irrational critical phases; no floating-point tests."""
from math import factorial
import sympy as S

x, L = S.symbols('x L')
f = 2*x**3-3*x-2
shifted = S.Poly(f.subs(x, x+1).expand(), x)
assert shifted.as_expr() == 2*x**3+6*x**2+3*x-3
assert shifted.LC() % 3 != 0
assert all(c % 3 == 0 for c in shifted.all_coeffs()[1:])
assert shifted.TC() % 9 != 0
assert S.discriminant(f, x) == -216
value = 256*L**3-2036001*L**2-17376*L-256
assert S.resultant(f, L*(2-x)**4-16*(1+x)**4, x) == 16*value
assert S.discriminant(value, L) < 0
assert S.Poly(value, L).nth(1)/S.Poly(value, L).LC() == -S.Rational(17376, 256)
assert f.subs(x, 1) < 0 and f.subs(x, S.Rational(3, 2)) > 0

discriminant_checks = 0
for m in range(2, 13):
    for epsilon in (1, -1):
        F = S.Poly((m-x)*(1+x)**m-epsilon*x**(m+1), x)
        assert S.expand(x*S.diff(F.as_expr(), x)-(m+1)*F.as_expr()
                        +m*(m+1)*(1+x)**(m-1)) == 0
        if epsilon == 1:
            assert F.degree() == m+1 and F.LC() == -2
            expected = (-1)**(m*(m+3)//2)*2*S.Integer(m)**m*(m+1)**(m+1)
        else:
            assert F.degree() == m-1 and F.LC() == m*(m+1)//2
            expected = ((-1)**((m-1)*(m-2)//2)*2
                        *S.Integer(m)**(m-3)*(m+1)**(m-2))
        assert S.discriminant(F.as_expr(), x) == expected
        assert S.gcd(F, F.diff()).degree() == 0
        discriminant_checks += 1
print(f'PASS: shifted Eisenstein cubic, value resultant, and {discriminant_checks} discriminant/derivative identities.')

expected_groups = {
    (1, 2): (6, 3, False), (1, 3): (24, 12, True),
    (1, 4): (120, 60, True), (1, 5): (720, 360, True),
    (-1, 3): (2, 1, False), (-1, 4): (6, 3, False),
    (-1, 5): (24, 12, True), (-1, 6): (120, 60, True),
    (-1, 7): (120, 60, True),
}
for (epsilon, m), expected in expected_groups.items():
    F = S.Poly((m-x)*(1+x)**m-epsilon*x**(m+1), x)
    assert F.is_irreducible
    G, inside_alternating = S.polys.numberfields.galois_group(F)
    H = G.derived_subgroup()
    d = F.degree()
    double_transitive = H.is_transitive() and len(H.stabilizer(0).orbit(1)) == d-1
    actual = (G.order(), H.order(), double_transitive)
    assert actual == expected, (epsilon, m, actual, expected)
    print(dict(branch='plus' if epsilon == 1 else 'minus', m=m, degree=d,
               group_order=G.order(), commutator_order=H.order(),
               commutator_doubly_transitive=double_transitive), flush=True)

# The degree-six minus branch at m=7 is not S6 and contains i via its discriminant.
F7 = (7-x)*(1+x)**7+x**8
disc7 = S.discriminant(F7, x)
assert disc7 < 0 and S.sqrt(-disc7).is_Integer
assert expected_groups[-1, 7][0] != factorial(6)
print('PASS: all nine exact group checks; the m=7 minus branch is not S6.')
