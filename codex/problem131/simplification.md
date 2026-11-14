# Deriving Apéry's recurrence by matching poles

Status: proved algebraically below and checked by exact symbolic computation. This is a simplification of the *derivation and verification of Zudilin's first certificate*, not a new proof of irrationality. No claim of historical novelty. The full irrationality argument still needs the arithmetic lemma, positive small linear forms, and the identity with Ball's series.

Primary source: W. Zudilin, *An elementary proof of Apéry's theorem*, arXiv:math/0202159v1 (17 February 2002), Lemmas 2–3, formulas (5)–(6). Read locally as `codex/sources/zudilin-0202159.txt`.

## What is mysterious, and what can be removed?

Zudilin gives a quadratic polynomial produced by creative telescoping and then asks the reader to check a rational identity. Instead, we can **discover both that polynomial and the recurrence coefficients together**, using only the two possible poles and infinity. This is a concrete reduction in unexplained calculation. It connects recurrence discovery to the elementary principle that a rational function is determined by its principal parts and its polynomial part.

Fix an integer n >= 1 and put

    R_n(t) = ((t-1)...(t-n) / [t(t+1)...(t+n)])².

Seek constants A,B,C (depending on n) and a quadratic q(t)=a t²+b t+c such that

    A R_{n+1}(t)+B R_n(t)+C R_{n-1}(t)
        = q(t+1)R_n(t+1)-q(t)R_n(t).                 (1)

Why a quadratic? The eventual proof differentiates and sums this identity. Since R_n(t)=O(t^-2), a quadratic gives q(t)R_n(t)=constant+O(t^-1), whose derivative tends to zero. Thus quadratic polynomials are the largest polynomial class with the desired derivative boundary behavior. This does not prove every possible telescoper must be polynomial; it motivates this small, sufficient search space.

The shift ratios are immediate from the products:

    R_{n+1}/R_n = ((t-n-1)/(t+n+1))²,
    R_{n-1}/R_n = ((t+n)/(t-n))²,
    R_n(t+1)/R_n(t) = t⁴/[(t-n)²(t+n+1)²].

Divide (1) by R_n. Its desired form becomes

    A ((t-n-1)/(t+n+1))² + B + C ((t+n)/(t-n))²
      = q(t+1)t⁴/[(t-n)²(t+n+1)²]-q(t).            (2)

Both sides have possible double poles only at t=n and t=-n-1. Match their principal parts, then their finite value at infinity. This suffices for a rational identity: the difference has no finite poles and vanishes at infinity, so it is zero.

## First pole: t=n

Only the C-term on the left has a pole. On the right only the first term does. Match values and derivatives of their numerators after multiplying by (t-n)². The value equation is

    C = n² q(n+1) / [4(2n+1)²].                    (3)

The derivative equation, after using (3), is

    n(2n+1) q'(n+1) + (4n+3) q(n+1) = 0.         (4)

These equations can be obtained by ordinary product differentiation; no asymptotic theorem or summation algorithm is needed. Equivalently, when q(n+1) is nonzero, logarithmic differentiation of q(t+1)t⁴/(t+n+1)² gives the same condition. The polynomial form (4) remains valid without division by q.

## Second pole: t=-n-1

The corresponding value equation is

    A = (n+1)² q(-n) / [4(2n+1)²],                (5)

and the derivative equation is

    (n+1)(2n+1) q'(-n) - (4n+1)q(-n) = 0.        (6)

Now (4) and (6) are just two linear equations in a,b,c. Their coefficient matrix has rank two: the determinant of its a,c columns is

    -3(2n+1)³,

which is nonzero for n >= 1. Substitution verifies that their one-dimensional solution space is

    q(t) = lambda [-2t²+t+(2n+1)²].               (7)

Thus this quadratic is forced by the two pole conditions; it is not an arbitrary lucky guess.

Equations (3) and (5) now simplify to

    C = lambda n³/[4(2n+1)],
    A = lambda (n+1)³/[4(2n+1)].

The cubic coefficients of the recurrence have appeared naturally from the pole locations and pole orders.

## Infinity determines the middle coefficient

Expand only through the constant term:

    t⁴/[(t-n)²(t+n+1)²]
      = 1 - 2/t + [2n(n+1)+3]/t² + O(t^-3).

For q=a t²+b t+c this gives

    lim_{t→infinity} [q(t+1)t⁴/((t-n)²(t+n+1)²)-q(t)]
      = 2a n(n+1)-b.

With a=-2lambda and b=lambda, this is -lambda(2n+1)². Therefore

    A+B+C = -lambda(2n+1)².

Choose the convenient normalization lambda=4(2n+1). We obtain

    A=(n+1)³,    C=n³,
    B=-4(2n+1)³-(n+1)³-n³
     =-(2n+1)(17n²+17n+5),
    q(t)=4(2n+1)[-2t²+t+(2n+1)²].

This proves exactly Zudilin's identity (6). The coefficient 34 in its expanded middle term is explained by a simple combination of three cubic expressions, not introduced by fitting numerical data.

## From the identity to a recurrence for linear forms

Define F_n=-sum_{t=1}^infinity R_n'(t). Differentiate (1) and sum to N. The right side telescopes to

    (q R_n)'(N+1) - (q R_n)'(1).

The first term tends to zero because qR_n=constant+O(t^-1). The second term is zero because R_n has a double zero at t=1 for n>=1. The differentiated series converge absolutely. Hence

    (n+1)³ F_{n+1} -(2n+1)(17n²+17n+5)F_n+n³F_{n-1}=0.

This argument does **not** justify separately equating the rational and zeta(3) coefficients before irrationality is known. That separation would be circular. Obtain the coefficient recurrence by a separate residue argument, as Zudilin's final-page footnote suggests, or work directly with the F sequence.

## What has been simplified, and what remains?

- Proved: the short telescoping certificate and all three recurrence coefficients are derivable from two double-pole conditions and infinity.
- Proved: no nonzero linear or constant q works within this polynomial ansatz; (7) has degree two whenever lambda is nonzero.
- Exact computational check: the symbolic rational identity is identically zero, and the full coefficient system has one-dimensional nullspace over Q(n). Run `python codex/computations/check_pole_certificate.py`.
- Not proved here: a correspondingly short derivation of the larger Ball-series certificate; a recurrence-free equality between the two linear forms; a new irrationality proof; historical novelty.

## Next concrete simplification experiments

1. Apply the same pole-matching method to Ball's rational summand. Its adjacent-n ratios have more poles and zeros; derive the required denominator of its certificate before choosing the numerator ansatz. This may explain the degree-six polynomial even if it does not shorten it.
2. Compare with the integral route: integration by parts is continuous telescoping, whereas (1) is discrete telescoping. Ask whether a change of representation makes the certificate shorter while preserving denominator control.
3. Try the same method with modified numerator roots. Record whether a quadratic certificate still exists, rather than assuming every rational approximation construction has a second-order recurrence.

These are proposed next tasks, not completed results.
