# Independent review: the critical quarter is elementary on polynomials

Status: accepted. The full theorem in `newton-hilbert-critical-quarter.md`,
its original standard-library checker, and the separate x-coordinate checker
pass this independent audit. No spectral closure, completeness assertion,
measure uniqueness, or historical novelty claim is needed for maximality.

The square identity is exact. Expanding its first square gives the original
Dirichlet term, the cross integral ∫x(1−x)ff′, and one quarter of∫(1−x)f².
The cross integral is−(1/2)∫(1−2x)f² with zero endpoint terms. Adding the
last quarter of∫xf² leaves precisely∫(x−1/4)f². The last integral in the
square representation is strictly positive for every nonzero polynomial,
so positivity at the critical value and below is strict.

The central-binomial recurrence gives the telescoping identity coefficient by
coefficient in u=1−x. It does not require convergence of an infinite series.
The induction bound a_N²≤1/(N+1) is valid because the needed polynomial
difference is3N+2≥0. The separate use of the positive binomial series to bound
f_N≤x^(−1/2) is valid on0<x≤1; the exceptional endpoint x=0 has no effect on
the integral and each finite trial is a polynomial.

The exact first-square integral is

    (2N+1)²a_N²/[8(N+1)]≤1/2,

and the other integral is≤1/4. Thus q(f_N)≤3/4 uniformly, with no unjustified
passage of an infinite squared-norm limit through integration.

The norm bounds are also finite and exact. The full central-binomial
convolution has coefficient1 in each degree. Truncating the two factors atN
retains that complete coefficient for degrees≤N and a positive subsum≤1 for
the remaining degrees through2N. Integrating these finite polynomials gives

    H_(N+1)≤||f_N||²≤H_(2N+1).

Divergence of the harmonic series therefore makes the positive critical
Rayleigh quotients tend to0. Every δ>0 admits a finite N with
H_(N+1)>3/(4δ), giving q(f_N)−δ||f_N||²<0. This proves necessity for every
forbidden shift, not merely the four shifts in the sample checker. Scaling
by c>0 uses δ=(λ−c/4)/c exactly as stated.

All-polynomial positivity is therefore equivalent to λ≤c/4 in the specified
a₀=c>0 family. A nonnegative-support representing measure must have that
positivity, so larger shifts are excluded without using a spectral theorem.
The theorem does not assert an exhaustive classification for general a₀/c.
Existence/density/essential selfadjointness of the classical spectral model
remain separate facts; the elementary maximality argument neither assumes nor
establishes them.

The original checker was run unchanged and passed41 exact cases, four negative
shift witnesses, and mixed-sign forms. The independent checker
`check_newton_critical_quarter_independent.py` constructs the trial polynomials
in x rather than the original u coordinate and verifies fifteen direct
Sturm–Liouville forms, both norm bounds, coefficient convolutions, the induction
inequality, and all applicable shifted witnesses. It also passes.

The classification note has been updated to cite this elementary proof for
maximality, while retaining Ismail–Koelink as the source for the classical
operator/spectral-density identification.
