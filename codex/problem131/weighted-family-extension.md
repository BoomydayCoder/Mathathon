# Three-dimensional positive weights: actual freedom and a finite improvement

Status: exact structural lemmas and finite computations. This extends the family in weighted-family-arithmetic.md and produces a better explicit small form at n=2, but no sequence proving irrationality. No novelty claim.

## 1. The extended family

Keep the sextic-pole function H_n(t) from orthogonality-design.md. Put v=t+n/2 and take

    W_(C,D,E)(t)=(2t+n)[C v^4+D v²+E] H_n(t),
    C,D,E rational, n>=1.

The three basis weights correspond to 2v^5,2v^3,2v. Their series converge, since H_n(t)=O(t^(-2n-6)) and the largest multiplier has degree5. Reflection t->-n-t makes each series a rational linear form in 1,zeta(3),zeta(5), with all other zeta coefficients canceled.

Write the three coefficient triples as

    (a_5,a_3,a_0), (b_5,b_3,b_0), (c_5,c_3,c_0),

where the convention is coefficient_5*zeta(5)+coefficient_3*zeta(3)-coefficient_0. Eliminating zeta(3) means

    a_3 C+b_3 D+c_3 E=0.                                (1)

This is one rational linear condition, leaving a two-dimensional space of isolated fifth-zeta forms. Fixing C=1 leaves one free parameter, but allowing integer triples makes the exact lattice formulation more informative.

## 2. Positivity is an exact quadratic condition

Let z=v². At positive integers t>n, z>=z_0=(3n/2+1)² and the remaining factors are positive. A sufficient condition for positive summands is

    q(z)=Cz²+Dz+E >=0 for every real z>=z_0,

with q not the zero polynomial. This condition is exactly characterized as follows:

- C<0: impossible.
- C=0: D>=0 and Dz_0+E>=0, with (D,E)!=(0,0).
- C>0 and -D/(2C)<=z_0: require q(z_0)>=0.
- C>0 and -D/(2C)>z_0: require 4CE-D²>=0.

This is just the minimum of a quadratic on a half-line. It is stronger than positivity only at the discrete summation arguments, but easy to check without numerical estimates. Equality at finitely many points is harmless; a nonzero nonnegative quadratic is positive at infinitely many later arguments, so the full sum is strictly positive.

## 3. The cancellation lattice is saturated and computable exactly

Clear denominators from (1) and divide by their gcd to get a primitive integer row (a,b,c). Suppose g=gcd(b,c)>0, h=gcd(a,g), and choose integers u,v with ub+vc=g. Then

    e_1=(g/h,-(a/h)u,-(a/h)v),
    e_2=(0,c/g,-b/g)                                    (2)

form a full integer basis for the solutions of aC+bD+cE=0.

Proof: every solution has aC divisible by g, hence C divisible by g/h. Subtract the corresponding multiple of e_1. The remaining solution has C=0, so (D,E) is an integer multiple of (c/g,-b/g). This proves saturation: no admissible integer weight is missed by the basis. QED.

The implementation applies exact two-dimensional Gauss reduction to this basis, then considers integer combinations s e_1+t e_2 with -4<=s,t<=4, removes common coefficient gcds and duplicates, and retains weights satisfying the positivity criterion. This is a small bounded search, not an exhaustive search over the full lattice.

Two selection rules use only rational/integer data:

1. `small_weight`: smallest squared Euclidean norm C²+D²+E²;
2. `small_pair`: smallest max(|U|,|V|) after reducing the isolated form to a coprime integer pair U zeta(5)-V.

Neither rule uses a numerical value of zeta(5) to select a weight. Only after selection are the actual errors evaluated.

## 4. An explicit improved form

For n=2 the cancellation row is

    -4805 C-2383 D+28935 E=0.

The positive integer weight

    (C,D,E)=(147,165,38)

satisfies this equation exactly and is positive for all z>=0. Its primitive fifth-zeta form is

    12544 zeta(5)-13007.

This lies strictly between zero and one. A completely rational certificate follows from

    S_M+1/[4(M+1)^4] < zeta(5) < S_M+1/(4M^4),
    S_M=sum_(k=1)^M 1/k^5,

which follows from integral comparison of the decreasing positive tail. At M=20, substituting into the displayed form yields the exact bounds

    0.22016396... < 12544 zeta(5)-13007 < 0.22363900... <1.

The decimals are merely summaries; the script checks the two bounding fractions exactly.

This is a concrete improvement over the previous two-dimensional family's n=2 primitive form, which had magnitude greater than21. It is not an irrationality proof: one small nonzero integer-coefficient linear form exists for every irrational number and for many rational numbers, and no tending-to-zero sequence has been established here.

## 5. Bounded search results

The exact coefficient maps have nonzero determinant for every n=1,...,12 tested. Selected results of the `small_pair` rule are:

| n | log absolute primitive error / n | comment |
|---|---|---|
| 1 | 0.25883654 | primitive pair (8,7) |
| 2 | -0.75307861 | explicit form above |
| 3 | 1.11660625 | error greater than1 |
| 4 | 1.39945265 | error greater than1 |
| 5 | 1.89898677 | error greater than1 |
| 6 | 1.58757598 | error greater than1 |
| 8 | 1.75972320 | error greater than1 |
| 10 | 1.69369460 | improvement over two-dimensional rate5.11627 at the same n |
| 12 | 2.16746132 | error greater than1 |

These logarithms are numerical evaluations after exact coefficient and gcd calculations. Except for the separately certified n=2 example, they are not presented as rational error enclosures. The selected forms still fail to show eventual decay. The search is too small to prove optimality or an asymptotic obstruction.

The weights do not remain small: selected coefficients have grown to many digits by n=10. Thus the extra dimension helps finite cancellation and gcds but has not produced a bounded-complexity multiplier formula or a proved improved denominator budget.

## 6. Why representability alone is not a breakthrough

**Lemma.** If the 3x3 rational coefficient map above is invertible, every rational pair (U,V) occurs as an isolated form U zeta(5)-V for some rational weight (C,D,E).

Proof: solve the invertible linear system with target triple (U,0,V). Clearing coefficient denominators gives an integer weight and scales the target pair; reducing the pair restores its primitive normalization. QED.

Therefore merely producing an isolated zeta(5) form in this three-dimensional space is arithmetically uninformative unless the construction also controls positivity, coefficient size, and decay. Using a known high-accuracy approximation to zeta(5) as the target of this system would simply encode that approximation, not derive an irrationality argument.

The positivity constraint is substantial. At n=1 the exact inverse map, with output (Z,0,R), is

    C=(8912Z-8576R)/1408,
    D=(3320Z-3392R)/1408,
    E=(133Z-152R)/1408.

Any nonnegative half-line weight must satisfy C>=0 and q(25/4)>=0. Since

    1408 q(25/4)=369008Z-356352R,

these necessary conditions imply

    if Z>0: R/Z<=23063/22272,
    if Z<0: R/Z>=557/536.                               (3)

The interval between these two rational bounds contains zeta(5). This can be proved without high-precision evaluation:

    sum_(k=1)^4 1/k^5 >23063/22272,
    sum_(k=1)^3 1/k^5+1/324 <557/536.

The first sum is below zeta(5), and the second expression is above it. Thus the positive-weight cone excludes all sufficiently accurate slopes near zeta(5) at n=1, despite algebraic representability of every slope without positivity.

More quantitatively, take the positive rational number delta equal to the minimum of the two gaps in the preceding rational inequalities. Every isolated form arising from a nonnegative half-line weight at n=1 and normalized to primitive U>0 then satisfies

    |U zeta(5)-V| >= delta U >= delta.

Hence varying coefficients indefinitely at this fixed n cannot force the forms to zero. This is a proved restriction for n=1; no uniform positive gap in n is asserted.

## 6a. A general rational sampling bound for fixed polynomial degree

The fixed-n gap is not peculiar to n=1. It follows from positivity and finite-dimensional interpolation.

**Lemma.** Fix n and polynomial degree d for which all weighted sums converge. Put z_t=(t+n/2)² and w_t=(2t+n)H_n(t)>0 for t>n. Suppose q is a rational polynomial of degree at most d with q(z_t)>=0 for all integer t>n. Let S(q)=sum_(t>n)w_t q(z_t), and let Z(q) be its fifth-zeta coefficient. Choose d+1 distinct summation indices t_i>n and their Lagrange basis polynomials L_i(z), so L_i(z_(t_j))=1 if i=j and 0 otherwise. Define

    M=max_i |Z(L_i)|/w_(t_i).

If the fifth-zeta coefficient functional is nonzero, M>0 and

    S(q)>=|Z(q)|/M.                                    (5)

All quantities defining M are rational and explicitly computable from partial fractions.

Proof: interpolation gives q=sum_i q(z_(t_i))L_i. Therefore

    |Z(q)|<=sum_i q(z_(t_i))|Z(L_i)|
            <=M sum_i w_(t_i)q(z_(t_i))<=M S(q).

The first inequality uses nonnegativity of the sampled values and the last uses nonnegativity of every summand. QED.

When the third-zeta coefficient is canceled and the resulting pair is normalized to primitive integers U,V with U>0, equation (5) gives

    |U zeta(5)-V|>=U/M>=1/M.

Thus no sequence obtained by varying only coefficients within a fixed n and fixed polynomial degree can have primitive linear forms tending to zero, as long as positivity is proved in this summandwise way. This is a uniform bound over all coefficient choices, not just the bounded lattice search. It allows arbitrary rational multipliers and gcd reduction.

For the quartic-in-v family, take d=2 and t_i=n+1,n+2,n+3. The lemma does not obstruct n tending to infinity: M depends on n and can grow rapidly. Its role is to rule out the misleading idea that the algebraic freedom at one fixed n permits arbitrarily good positive constructions.

## 7. Next mathematical direction

The useful remaining task is to derive formulas for the two boundary slopes of the positive cone as functions of n, then compare the arithmetic of those slopes with their distance to zeta(5). This is more informative than enlarging an unstructured coefficient search. A dimension-three construction offers genuine freedom, but its useful portion is the intersection of an exact cancellation lattice with a positivity cone, followed by primitive gcd normalization.

Run `python codex/computations/weighted_family_extension.py`. It writes exact coefficient triples, map determinants, reduced lattice bases, chosen integer weights and primitive pairs to the accompanying JSON. The script verifies the n=2 rational error enclosure and the n=1 positivity-gap inequalities exactly. Requires mpmath only for reported logarithms.

Independent audit: odd-zeta-independent-review.md section 9 independently expands the n=2 weighted rational function into principal parts, recovers the primitive pair (12544,13007), and checks the exact rational enclosure between zero and one. It also checks the n=1 inverse-map cone bounds, rational zeta brackets, saturated lattice basis, and quadratic positivity criterion. The full bounded search was not independently rerun. The subsequently added general sampling lemma is separately marked for review.
