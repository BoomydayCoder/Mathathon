# Primary-source comparison for the Rodrigues / Eulerian construction

**Date:** 2026-11-14 UTC. This bounded source investigation was performed after
our independent derivations. It finds exact specializations of established
constructions, not merely similarities. It does not establish priority for the
specific Laurent coordinates, critical-point identities, or numerical rates.
No novelty conclusion follows from not locating those particular formulas.

## 1. Main attribution result

The squared rational function, derivative sum and pure odd-zeta form are an
exact specialization of Nesterenko's 2003 general construction. His Barnes
representation theorem also contains the finite Fourier decomposition. The
Eulerian numbers give an explicit combinatorial name/formula for its
coefficients. The underlying Rodrigues polynomial and its confluent logarithmic
orthogonality occur in the classical Jacobi–Piñeiro framework, explicitly
including the confluent case in Smet–Van Assche's work.

Accordingly, describe our work as an independently developed, explicit study
of this specialization: a conceptual construction, exact Laurent geometry,
sharp contour estimates and reproducible certificates. Do not describe the
rational-function family, multiple orthogonality, or finite Fourier/Barnes
method as a historically new method.

## 2. Exact dictionary into Nesterenko's rational-function construction

Primary source: Yuri V. Nesterenko, *Integral identities and constructions of
approximations to zeta-values*, Journal de Théorie des Nombres de Bordeaux
**15** (2003), 535–550, [DOI 10.5802/jtnb.412](https://doi.org/10.5802/jtnb.412),
[primary PDF](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.412.pdf).
Relevant locations: pp.536–538, equations (3)–(10), Proposition1 and the
observation A₁(1)=0 when δ≥2.

Write our squared Mellin function as

\[
 R_{n,m}(t)=\kappa_{n,m}^2
 \frac{\Gamma(t)^{2m+2}}{\Gamma(t-n)^{2m}\Gamma(t+mn+1)^2},\qquad
 \kappa_{n,m}=\frac{(mn)!}{(n!)^m}.
\]

In Nesterenko's notation, choose the following parameters; M here replaces his
number of gamma factors so it is not confused with our multiplicity m:

| Source parameter | Our specialization |
|---|---|
| number M of gamma factors | 2m+2 |
| a₁,…,a_M | all n+1 |
| b₁,…,b₂ₘ | all 1 |
| b₂ₘ₊₁,b₂ₘ₊₂ | both (m+1)n+2 |
| source variable s | t−n−1 |
| derivative parameter μ | 2m |
| sum-shift parameter a in Proposition1 | n+1 |
| polylogarithm variable z | 1 |

The source's normalized polynomial bricks contribute (n!)^(−2m), and its two
reciprocal bricks contribute ((mn)!)². Its normalization γ is therefore exactly
κ², and its R(s) is Rₙ,ₘ(s+n+1), not just a proportional function. Its degree
difference is δ=2 and its maximum pole order is q=2.

Proposition1 now gives exactly

\[
 G_{2m}(1)=-\frac1{(2m-1)!}\sum_{\nu\ge1}R_{n,m}^{(2m-1)}(\nu)
          =\mathscr L_{n,m}.
\]

Only ζ(2m) and ζ(2m+1) can occur; δ=2 makes the ζ(2m) coefficient zero.
At a double pole t=−k the coefficient is
binom(mn,k)²binom(n+k,k)^(2m). The source's equation (10), whose binomial
multiplier is μ for double poles, therefore gives our exact coefficient

\[
 2m\sum_{k=0}^{mn}\binom{mn}{k}^{2}\binom{n+k}{k}^{2m}.
\]

This establishes containment of the exact family in the published theorem.
It does **not** establish that Nesterenko separately tabulated or optimized this
particular one-parameter choice. Those are distinct historical assertions.

## 3. The Fourier decomposition is also inside the primary theorem

The same paper, pp.545–547, §3.2, equation (29), Theorem1 and Lemma4,
represents derivative sums by finitely many integrals

\[
 I_r(u)=\frac1{2\pi i}\int R(s)
           \left(\frac\pi{\sin\pi s}\right)^r e^{\pi ius}\,ds.
\]

At z=1, take r=2m and u=0. Only the top-order part survives, giving
G₂ₘ(1)=Σⱼd₂ⱼ I₂ₘ(2j), with 1−m≤j≤m−1. Under the preceding integer shift,
sin^(2m)(πs)=sin^(2m)(πt) and exp(2πijs)=exp(2πijt), so

\[
 I_{2m}(2j)=T_{n,m,j},\qquad
 \boxed{d_{2j}=\frac{E(2m-1,j+m-1)}{(2m-1)!}.}
\]

The coefficient identification follows by differentiating cotangent or using
the Eulerian generating polynomial. Thus our named Eulerian formula is a
closed expression for the coefficients of the published finite decomposition;
the decomposition itself is not an unpreceded mechanism.

A second primary antecedent is Zudilin's 2001 preprint
[*Irrationality of values of zeta-function*, math/0104249](https://arxiv.org/abs/math/0104249),
§4, Lemmas2–3 and equations (11)–(18). It already uses cotangent derivatives,
reflection, a finite trigonometric/Fourier polynomial, and complex saddles.
Its chosen rational function in §3 has symmetric positive/negative zero
blocks and different pole multiplicities, so it should not be equated with
our parameter choice. The analytic mechanism is nonetheless directly shared.

## 4. The positive central period is another exact Nesterenko specialization

The multiple-integral identity in Nesterenko's Theorem2, pp.547–550, is also
stated cleanly as Proposition1, p.9, of Zudilin's
[*Arithmetic of linear forms involving odd zeta values*, math/0206176v2](https://arxiv.org/abs/math/0206176).
Use the latter notation with dimension M=2m+1, r=2m−1, a₀=a₁=…=a_M=n+1,
first 2m−1 parameters bᵢ=2n+2, and last two bᵢ=(m+1)n+2. Its gamma prefactor
is κ². Shifting its contour variable by n+1 gives Hₙ,ₘ and exactly Tₙ,ₘ,₀.

Equivalently the resulting positive unit-cube integral is

\[
 \int_{[0,1]^{2m+1}}
 \frac{\prod_{i=1}^{2m-1}v_i^n(1-v_i)^n\,
       x^n(1-x)^{mn}y^n(1-y)^{mn}}
 {\{\prod_i(1-v_i)+xy\prod_i v_i\}^{n+1}}
 \,d\mathbf v\,dx\,dy.
\]

All theorem inequalities hold: one may choose its vertical line between
−n−1 and 0. This also follows directly from our beta convolution by
uᵢ=vᵢ/(1−vᵢ). The dual Laurent coordinate description is another explicit
presentation of this established integral specialization.

## 5. Exact polynomial identification: confluent Jacobi–Piñeiro

The original Rodrigues source already checked in this project is Van Assche
and Coussement, *Some classical multiple orthogonal polynomials*, JCAM
**127** (2001), 317–347, [math/0103131](https://arxiv.org/abs/math/0103131),
§3.1, equation (3.24). It provides the commuting Rodrigues operators.

A still closer primary source is Christophe Smet and Walter Van Assche,
*Mellin transforms for multiple Jacobi–Piñeiro polynomials and a q-analogue*,
J. Approx. Theory **162** (2010), 782–806; the inspected
[author preprint](https://wis.kuleuven.be/analyse/members/walter/mellin.pdf)
is dated August6,2008 and has the shorter title without “and a q-analogue.”
Theorems2.2,2.4 and2.8 explicitly cover the Mellin numerator, coalescing
logarithmic weights, and Rodrigues formula.

For Theorem2.4 take r=m, α=0, balanced multiindex (n,…,n), and s=t−1.
Normalize the source's monic polynomial to constant term1, multiplying by

\[
 (-1)^{mn}\binom{(m+1)n}{n}^{m}.
\]

The result is exactly Jₙ,ₘ, and its Mellin transform is exactly Fₙ,ₘ.
Thus the repeated logarithmic orthogonality is explicitly classical, not
merely a conjectural limiting interpretation. Our elementary uniqueness and
Rolle proofs remain self-contained ways of explaining this specialization.

The author PDF was readable through the web tool's cached primary-source
extraction, but direct local download returned HTTP503 twice. Its formulas
and theorem numbers were inspected; no local PDF hash is claimed for it.

## 6. Relations that should not be conflated

**Rivoal.** The primary paper [math/0008051](https://arxiv.org/abs/math/0008051),
pp.1–3, uses a reflected pair of numerator blocks and pole order a, with a
free parameter r. Lemma1 uses symmetry to remove even zeta values. Our pure
form instead uses only double poles, a derivative of order2m−1, and vanishing
of the sum of simple-pole coefficients. Both are inside the general
rational-function framework; the explicit parameter families and arithmetic
objectives differ.

**Vasilyev-type integrals.** Zudilin's primary
[math/0206177](https://arxiv.org/abs/math/0206177), Theorem(4) and its remarks,
relates nested-denominator cube integrals to very-well-poised series. The
standard symmetric specialization has a rational linear combination of several
odd zeta values. Our positive central period has a product/sum denominator
and already at n=0 includes π²ζ(3) in the m=2 case. Do not identify the two
families merely because both generalize Beukers integrals. For example,
the source's theorem gives standard J₅,₀=2ζ(5), while our T₀,₂,₀ is
4ζ(5)+(4π²/3)ζ(3); even after the natural leading-coefficient normalization,
the extra positive term remains.

**Oscillation.** Fischler's primary paper
[*Nesterenko's criterion when the small linear forms oscillate*, 1201.2651](https://arxiv.org/abs/1201.2651),
Theorem1 and Proposition1, treats cosine-modulated exponential forms and
subsequences away from zeros. Our nonvanishing/sign analysis belongs in this
established setting. That criterion requires sufficiently small **integer**
linear forms, so it cannot turn our analytically small forms into an
irrationality theorem when denominator clearing overwhelms the decay.

## 7. What this bounded search did and did not establish

Exact published containment has been established for the polynomial/Mellin
family, pure derivative forms, finite Fourier/Barnes representation, and
positive multiple-integral specialization. The rational dictionaries above
are reproducible in `../computations/check_rodrigues_source_dictionary.py`.

The selected primary sections did not supply an explicit matched occurrence
of our Laurent critical polynomial, the product laws for its critical values,
or the sharp constants and certified dominant saddles for m=2,3,4,5. These
items remain **historically unattributed**, not established as new. A broader
investigation could examine general hypergeometric asymptotics, multiple
orthogonal-polynomial asymptotics, and parameter-specialization tables.

Searches covered the precise binomial pattern, the m=2 critical polynomial,
cotangent/Eulerian terminology, the named authors, and Jacobi–Piñeiro Mellin
transforms. Search-engine misses, especially for mathematical formulas, are
weak evidence and are not used for a novelty conclusion. Full downloaded
primary sources and SHA256 metadata are recorded in
`../sources/rodrigues-primary-context-provenance.json`.
