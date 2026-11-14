# Primary-source comparison: the centered Hurwitz generalized Stieltjes family

Status: bounded literature comparison, 2026-11-14. The primary sources below
were retrieved and relevant statements read. This is not an exhaustive priority
search. No historical novelty claim follows from the absence of an exact match
in this small set of papers.

Our target is

\[
 F_\alpha(z)=2\zeta(2\alpha+1,\sqrt z+1/2),\qquad \alpha>0,
\]

with a positive generalized Stieltjes representation of exact order alpha.
See [the real-order proof](real-order-positive-stieltjes-zeta.md) and
[the integer density proof](higher-order-positive-stieltjes-zeta.md).

## 1. What is already part of the established framework

D. Karp and E. Prilepkina, *Generalized Stieltjes functions and their exact
order*, Journal of Classical Analysis (2012), 53–74,
[DOI 10.7153/jca-01-07](https://doi.org/10.7153/jca-01-07).
The archived author manuscript is titled *Generalized Stieltjes transforms:
basic aspects*, [arXiv:1111.4271v2](https://arxiv.org/abs/1111.4271).
Archive: `codex/sources/karp-prilepkina-2012-generalized-stieltjes.pdf` and `.txt`.
Numbering below follows that manuscript.

- Definition (1) includes a nonnegative constant and a possible atom at zero,
  as well as the positive measure integral with kernel (z+t)^(-alpha).
- Theorem 3, equations (10)–(13), changes the representing order of the SAME
  function by a left Riemann–Liouville fractional integral of its measure.
- Theorem 4, equation (14), and Remark 3 give the corresponding right
  Kober–Erdelyi integral for the alternative kernel (1+xz)^(-alpha).
- Section 4, equation (24), defines the exact order as inf{alpha:f in S_alpha}.
  Theorem 15 proves attainment; Theorem 16 gives a measure criterion.
  Corollary 2 proves exactness for compactly supported representing measures.

These supply the correct terminology and the fractional-integral ancestry.
Our compact-mass asymptotic lower bound is an elementary way to prove exactness
in the particular finite-mass family; it is not a new definition or general
criterion. A nonnegative constant in their definition must vanish here because
F_alpha(z) tends to zero. An atom at zero is already covered by our compact-mass
argument, so their broader convention does not change our exact-order result.

There is a distinction to preserve: the standard order-changing representation
keeps f fixed. Our Beta product changes BOTH the centered Hurwitz exponent and
the measure. Thus Theorem 3 or 4 alone is not the displayed zeta identity.
Our actual family needs the additional centered square-root kernel calculation.
No centered Hurwitz formula was identified in this manuscript.

S. Koumandos and H. L. Pedersen, *On Generalized Stieltjes Functions*,
[DOI 10.1007/s00365-018-9440-8](https://doi.org/10.1007/s00365-018-9440-8),
Constructive Approximation 50, 129–144;
[arXiv:1706.00606](https://arxiv.org/abs/1706.00606).
Archive: `codex/sources/koumandos-pedersen-2018-generalized-stieltjes.pdf` and `.txt`.
Their introduction, equation (1), characterizes order-lambda transforms as
Laplace transforms of t^(lambda-1) phi(t), with phi completely monotone,
plus a nonnegative constant. Their principal results sharpen Sokal's real
variable characterization. This is directly relevant general theory, but the
read text does not supply our specific Hurwitz/square-root density.

## 2. An exact primary ancestor: the hyperbolic mixing distribution

Jim Pitman and Marc Yor, *Infinitely Divisible Laws Associated with Hyperbolic
Functions*, Canadian Journal of Mathematics 55(2) (2003), 292–330,
[DOI 10.4153/CJM-2003-014-x](https://doi.org/10.4153/CJM-2003-014-x).
Verified published PDF downloaded from
[Cambridge](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/S0008414X00033484).
Archive: `codex/sources/pitman-yor-2003-hyperbolic-published.pdf` and `.txt`.

Their equation (1), published page 292, states

\[
 E e^{-\theta^2 C_r/2}=(\cosh\theta)^{-r}.
\]

Their page 293 identifies C_1 as the first hitting time of {-1,1} by standard
Brownian motion. Consequently the exact function used in our positivity proof,

\[
 g(t)=\sech^2(\pi\sqrt t),
\]

is the Laplace transform of the already studied variable
Y=(pi²/2)C_2. This is an EXACT match to an ingredient, not merely a similar
hyperbolic formula. Independent increments identify C_2 with a sum of two
independent copies of C_1. Brownian scaling gives the equivalent sum of two
exit times from [-pi/sqrt(2),pi/sqrt(2)]. No Brownian machinery is needed for
our original Euler-product proof of complete monotonicity.

This source makes our Beta hierarchy conceptually transparent. Define the
common tilted mixing law

\[
 dP^*(Y)=\frac{\pi\sqrt\pi}{2}Y^{-1/2}dP(Y).
\]

Its normalization follows by Tonelli, not by an unverified special-function
moment identity:

\[
 E Y^{-1/2}=\frac1{\sqrt\pi}\int_0^\infty t^{-1/2}g(t)dt
 =\frac2{\pi\sqrt\pi}.
\]

If H_alpha is independent Gamma(alpha+1/2, rate 1), then our normalized
probability measure mu_alpha is the law of H_alpha/Y*.
For integer m this follows directly from the established density:

\[
 m\sigma_m(t)=mC_m\pi t^{m-1/2}E[Y^m e^{-tY}],\qquad
 mC_m\pi\Gamma(m+1/2)=\pi\sqrt\pi/2.
\]

For real alpha the Beta-product construction and elementary beta–gamma algebra
identify the same mixture. Explicitly, if B~Beta(alpha+1/2,m-alpha), independent
of H_m, then BH_m~Gamma(alpha+1/2,1). This explains why the parameters alpha+1/2
occur and why the coupling is independent of the selected larger integer m.
The gamma mixture and beta–gamma mechanism should be described as standard
probability constructions applied to this family, not as newly invented tools.

Pitman–Yor study Mellin transforms, zeta values, and hyperbolic distributions;
their Table 1 and subsequent sections contain substantial prior zeta structure.
The exact centered-Hurwitz generalized-Stieltjes theorem above was not
identified in the portions examined. A complete priority comparison would
need to trace these older distributional and special-function references.

## 3. Nearby results that must not be mistaken for the same theorem

C. Berg, S. Koumandos and H. L. Pedersen, *Nielsen's beta function and some
infinitely divisible distributions*,
[arXiv:1905.04131v2](https://arxiv.org/abs/1905.04131).
Archive: `codex/sources/berg-koumandos-pedersen-2019-nielsen-beta.pdf` and `.txt`.
Corollary 3.2 treats cosh(nu x)/cosh(x), its infinite divisibility, and Laplace
transforms expressed through Nielsen's beta function, belonging to S_2.
This is genuine hyperbolic/Stieltjes overlap, but its variable and target
function differ from F_alpha(z). Nielsen's beta is also not the Beta probability
law used in our product coupling.

D. Askitis and H. L. Pedersen, *Ratios of Entire functions and generalized
Stieltjes functions*, [arXiv:2107.00905](https://arxiv.org/abs/2107.00905).
Archive: `codex/sources/askitis-pedersen-2021-ratios-entire-stieltjes.pdf` and `.txt`.
Section 3, Proposition 3.1 and Corollary 3.2, treat derivatives of
log|f(a+i sqrt(x))|, with discrete Stieltjes measures at squared zeros.
This provides a particularly close-looking square-root/gamma precedent.
Its imaginary square root and log modulus are essential: it is not the real
centered argument sqrt(z)+1/2 in our Hurwitz function. Their multiple-zeta
section concerns Barnes functions, again a different family.

Also, the elementary identity zeta(s,z)=sum_(n>=0)(z+n)^(-s) is already a
positive order-s transform in z. It does not establish the order-(s-1)/2
result after the centered square-root substitution. General complete
monotonicity alone likewise does not establish exact generalized Stieltjes order.

## 4. Recheck of the integer additions

The minimal-order and Beta-product arguments added to the integer note were
rechecked separately from its earlier core peer review. No correction is needed.
This is a mathematical recheck by the note's originating agent, not a substitute
for a separate independent reviewer.

For minimal order, sum–integral comparison gives F_m(z)~z^(-m)/m. A nonzero
positive representing measure of order beta has finite positive mass c on
some compact interval [0,T], since finiteness of its transform controls compact
mass. Thus F_m(z)>=c(z+T)^(-beta), excluding beta<m. Infinite total mass is
allowed; no unjustified finite-mass assumption enters this lower bound.

For the Beta identity, put mu_m=m sigma_m and a=m+1/2. The density recurrence
is exactly

\[
 a x^{-a}\mu_{m+1}(x)=-D_x[x^{1-a}\mu_m(x)].
\]

Integration from t>0 to infinity gives

\[
 \mu_m(t)=a t^{a-1}\int_t^\infty x^{-a}\mu_{m+1}(x)dx,
\]

the product density for Beta(a,1) times an independent mu_(m+1) variable.
The infinity boundary vanishes exponentially. No boundary at zero is used;
indeed x^(1-a)mu_m(x) need not vanish there. Positivity, normalization, and the
product-density formula prove the distributional identity directly, without
assuming that equality of moments determines the measure.

## 5. Search scope and next steps

This bounded pass used direct primary PDF retrieval, arXiv title/abstract API
queries, Crossref bibliographic queries, and reference-following. Broad queries
for "Hurwitz zeta Stieltjes" are dominated by the unrelated Stieltjes constants;
those hits were not counted as evidence about generalized Stieltjes functions.
Google's returned challenge HTML was unusable. The publisher endpoint for
Karp–Prilepkina redirected in a loop, so its arXiv manuscript was used.
A guessed author-hosted Pitman URL returned a different document and was
removed; the verified published paper is the retained source.

The appropriate current wording is: **proved in this research notebook;
related general framework and an exact hyperbolic-distribution ancestor are
identified; priority of the complete centered-Hurwitz statement is unresolved.**
A next targeted bibliography pass should inspect the earlier Mellin-transform
references of Pitman–Yor and literature on centered polygamma transforms,
fractional Stieltjes inversion, and beta–gamma mixtures. It is unnecessary to
make any priority claim to use this representation pedagogically or to test
its rational-approximation arithmetic.
