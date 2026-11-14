# Exact classical identifications: Koecher–Leshchiner and Askey–Wilson

Status: primary PDFs retrieved, relevant formulas read, and the Askey–Wilson
coefficient pages visually inspected. The parameterized generating function
and the entire denominator recurrence have exact published antecedents.
Our independent finite-balance derivation remains useful exposition, but
these two objects must not be described as historically new.

This note records a bounded comparison, not an exhaustive priority search.
In particular it does not claim that our associated-numerator construction
or its recurrence is new, nor that a search establishes absence of an
older continued-fraction formula.

## 1. The generating function is exactly Koecher–Leshchiner

Put sigma=x². The target already used in
[the positive Wilson deformation](logarithm-zeta3-wilson-deformation.md) is

\[
 Z_\sigma=\sum_{n=1}^\infty\frac1{n(n^2-\sigma)}
          =\sum_{r=0}^\infty\zeta(2r+3)\sigma^r,\qquad|\sigma|<1.
 \tag{1}
\]

The exact classical identity is

\[
 \boxed{
 Z_{x^2}=\frac12\sum_{n=1}^\infty
 \frac{(-1)^{n-1}}{n^3\binom{2n}{n}}
 \frac{5n^2-x^2}{n^2-x^2}
 \prod_{j=1}^{n-1}\left(1-\frac{x^2}{j^2}\right).
 }                                                               \tag{2}
\]

It is printed as:

- T. Rivoal, equation (1.2), printed p.2 of the archived author manuscript;
- Kh. and T. Hessami Pilehrood (2008), equation (3), printed p.115;
- R. Tauraso (2025), equation (2), printed p.1.

All three attribute it to Koecher and independently Leshchiner; Rivoal and
the 2008 paper specify that Leshchiner obtained an expanded form.
The original cited works are M. Koecher, Letter to the editor,
Mathematical Intelligencer 2 (1980), 62–64, and D. H. Leshchiner,
"Some new identities for zeta(k)," Journal of Number Theory 13 (1981),
355–362. Those original two papers were not retrieved in this bounded pass;
the attribution here is explicitly grounded in the three retrieved papers.

There is no extra change of variables between (1) and our target:
the earlier resolvent is S_tau(1/4)=2Z_(tau²). Thus its tau=1/2 endpoint is

\[
 Z_{1/4}=8\log2-4.
\]

Differentiating (2) at sigma=0 gives the established type of Apéry-like
formula

\[
 \zeta(5)=\sum_{n=1}^\infty
   \frac{(-1)^{n-1}}{n^3\binom{2n}{n}}
   \left(\frac2{n^2}-\frac52H_{n-1}^{(2)}\right).                \tag{3}
\]

This is an explicit example of how a whole odd-zeta family arises by
coefficient extraction, with harmonic sums entering immediately.
It does not by itself establish adequate denominator estimates for an
irrationality proof. Rivoal's introduction explicitly distinguishes the
fast Apéry series from the Diophantine approximations needed for that proof.

## 2. The denominator recurrence is an exact Askey–Wilson specialization

Richard Askey and J. A. Wilson, "A recurrence relation generalizing those of
Apéry," Journal of the Australian Mathematical Society, Series A,
36 (1984), 267–278,
[DOI 10.1017/S1446788700024691](https://doi.org/10.1017/S1446788700024691),
give, in equation (22), printed p.271,

\[
 F_n={}_4F_3\left(
 \begin{matrix}-n,-n-a,n+b+e+1,n+c+f+1\\d+1,e+1,f+1\end{matrix};1
 \right),\qquad a+d=b+c.
\]

Their equation (23) is

\[
 A_n^{AW}F_{n+1}+B_n^{AW}F_n+C_n^{AW}F_{n-1}=0,
\]

with explicit rational coefficients on printed pp.272–273.

Make the exact substitution

\[
 \boxed{a=d=0,\quad b=-\tau,\quad c=\tau,\quad
                  e=\tau,\quad f=-\tau,\qquad \sigma=\tau^2.}   \tag{4}
\]

The balance condition is satisfied, and their F_n becomes exactly

\[
 F_n={}_4F_3\left(
 \begin{matrix}-n,-n,n+1,n+1\\1,1+\tau,1-\tau\end{matrix};1
 \right)
 =\sum_{k=0}^n
  \binom nk^2\binom{n+k}k^2
       \frac{(k!)^2}{\prod_{j=1}^k(j^2-\sigma)}
 =A_n(\sigma).                                                \tag{5}
\]

There is no normalization factor in this F_n identification. If one starts
instead from their binomial g_n in equation (19), their equation (21)
adds the known prefactor binom(n,tau)binom(n,-tau); that different sequence
should not be confused with F_n.

### Checking every coefficient, not just the family name

Write h=n+1. Under (4), the two auxiliary denominators printed on p.273 are

\[
 D_1=4h^3(3h^2-\sigma),\qquad D_2=4n^3(3n^2-\sigma).
\]

Their forward and backward numerators become respectively
h^4(h²-sigma)² and n^4(n²-sigma)². Consequently

\[
 A_n^{AW}=\frac12\alpha_n,\qquad
 C_n^{AW}=\frac12\gamma_n,
\]

where our coefficients are

\[
 \alpha_n=\frac{(n+1)((n+1)^2-\sigma)^2}
                  {2(3(n+1)^2-\sigma)},\qquad
 \gamma_n=\frac{n(n^2-\sigma)^2}{2(3n^2-\sigma)}.
\]

Their middle coefficient simplifies to

\[
 B_n^{AW}=
 -\frac{h(h^2-\sigma)(5h^2-\sigma)}{4(3h^2-\sigma)}
 -\frac{n(n^2-\sigma)(5n^2-\sigma)}{4(3n^2-\sigma)}
 -(2n+1)n(n+1)
 =-\frac12\beta_n,
\]

\[
 \beta_n=(2n+1)[3n(n+1)+1-\sigma]-\alpha_n-\gamma_n.
\]

Thus multiplying the printed Askey–Wilson recurrence by two gives exactly
[our generic-sigma diagonal recurrence](generic-sigma-diagonal-recurrence.md),
coefficient by coefficient. The checker substitutes into the full printed
formulas with symbolic n and tau; it is not merely a comparison of finitely
many sequence values.

This is also conceptual motivation: the hypergeometric series is balanced
(the lower parameters sum to one more than the upper parameters), terminates,
and is evaluated at one. Askey and Wilson explain why balance-preserving
contiguous relations lead to a three-term relation. Our neighboring-grid
balance is an elementary reconstruction of that same phenomenon.

### What this match does not prove

The published F_n identity proves the denominator recurrence. It does not,
without further work, identify the particular associated numerator b_n
defined from our positive measure and its shifted rational approximation.
No statement about b_n has been inferred from the fact that A_n satisfies
the recurrence. A separate mixed-form or Cauchy-transform proof remains
necessary. Any later proof of that bridge should cite its own argument
rather than promote the denominator match to an unproved second-solution
identification.

## 3. Other retrieved sources: useful leads with different scope

**Rivoal's two-variable unification.** His Theorem 1, equation (1.4), unifies
(2) and the Almkvist–Granville generating function by replacing the left
side with sum n/(n^4-a²n²-b^4). The right product has numerator
(j²-a²)²+4b^4. The proof reduces to a finite rational/combinatorial identity,
then proves that identity. This is a primary proof of the simultaneous
series generation, not a proof of our diagonal Padé numerator.

**Markov–WZ acceleration.** Hessami Pilehrood and Hessami Pilehrood,
Theorem 2, equation (14), printed p.119, produce a simultaneous series with
geometric rate 2^(-10). Their two summation paths applied to one Markov–WZ
pair explain why the same generating function has several accelerated
representations. Their Theorem 1 contains a second-order recurrence for a
different auxiliary sequence L_n. In the specialization yielding Koecher's
formula (B_0=1,A_0=C_0=0), that L_n is identically zero. It is therefore
not our nontrivial A_n recurrence. This prevents a misleading identification
based only on the appearance of a second-order recurrence in the paper.

**Tauraso's q-analogue.** The 2025 note proves a q-version by telescoping
finite q-harmonic identities, then summing and taking a limit. It explicitly
prints the undeformed identity (2) and its harmonic coefficient formulas.
It is a useful source for finite telescopers and q-deformation ideas, but
does not supply an identified equality with our generic-sigma diagonal
Cauchy-transform numerator.

**Cohen's general acceleration construction (1981).** His Theorem 2.2
constructs an Apéry-type triangular acceleration for Laplace transforms
of entire functions g of exponential type less than one. Our elementary
summand fits that setup exactly via

\[
 g(t)=\frac{\cosh(\tau t)-1}{\tau^2},\qquad
 \int_0^\infty e^{-nt}g(t)\,dt
       =\frac1{n(n^2-\tau^2)},\qquad |\tau|<1,
\]

with g(t)=t²/2 at tau=0. Section 3 uses central-factorial polynomials;
the conclusion identifies the finite-difference operator 2sinh(D/2),
in contrast with the ordinary derivative in Euler–Maclaurin summation.
This is a concrete classical lead for the Newton/finite-difference viewpoint.
However, on printed p.280 Cohen explicitly says that this paper has not
generalized Apéry's barycentric coefficients or his recurrence relations.
It should not be cited as if it already identified our A_n,b_n pair.

## 4. Continued fractions and the boundary of the comparison

For the Koecher summand t_n(sigma) on the right of (2), direct cancellation
gives

\[
 \frac{t_{n+1}}{t_n}
 =-\frac{n(n^2-\sigma)^2(5(n+1)^2-\sigma)}
 {2(n+1)^2(2n+1)((n+1)^2-\sigma)(5n^2-\sigma)}.                 \tag{6}
\]

Euler's series-to-continued-fraction transformation can therefore be applied
directly. If q_n=t_(n+1)/t_n, its finite truncations are exactly those of

\[
 \frac{t_1}{
  1-\cfrac{q_1}{1+q_1-\cfrac{q_2}{1+q_2-\cfrac{q_3}{1+q_3-\cdots}}}}.
 \tag{7}
\]

The equality with the finite partial sums follows by induction, before any
limit is taken. This gives a continued fraction for Z_sigma, but does not
identify it with the diagonal Padé recurrence; their equivalence or a
specific Bauer–Muir acceleration between them needs an additional argument.

Henri Cohen's arXiv:2304.11727, section 14, derives the sigma=0 Euler
continued fraction from the Markov–Apéry series and gives a Bauer–Muir
family of continued fractions with the same zeta(3) limit. Section 16
separately records Apéry's fast continued fraction. These are concrete
primary references for acceleration, but their displayed parameters keep
the target zeta(3) fixed. Our sigma changes the target Z_sigma. The retrieved
formulas have not been asserted to identify a Bauer–Muir path between the
full parameterized Koecher series and our full parameterized Padé pair.

A next focused proof task is to apply the balanced contiguous relation to
the associated Cauchy transforms or the corresponding mixed finite forms.
That would prove the numerator recurrence directly. A second task is to
compare the resulting convergent fraction with (7) through an explicit,
checkable Bauer–Muir transformation. Neither question requires any claim of
historical novelty.

## 5. Reproduction and source provenance

Run

    python codex/computations/check_koecher_askey_wilson_dictionary.py

It verifies the five full symbolic coefficient identities, nine terminating
hypergeometric identifications, a full symbolic summand ratio, seven exact
finite Euler fractions, five numerical Koecher identities, 45 numerical
summand-ratio evaluations, the zeta(5) coefficient, and the log(2) endpoint.

All PDFs are archived under codex/sources; paired .txt files were generated
with pdftotext. Original coefficient pages were additionally rendered and
read as images because the text extraction omitted some formulas.

| Archive basename | Retrieved primary URL | SHA-256 |
|---|---|---|
| rivoal-aperylike.pdf | https://rivoal.perso.math.cnrs.fr/articles/aperylike.pdf | df3ab8875f3e4f0af061f26f3a9a296c96d55600553014b10f12cda623b91623 |
| markov-wz-dmtcs424.pdf | https://dmtcs.episciences.org/424/pdf | 69a819f8ac4a186becd2595dd4904a773b97ddd4896306736c73d384fc9d957a |
| tauraso-2511.18034.pdf | https://arxiv.org/pdf/2511.18034 | e3ffe982a6096d5a23ff64fbe016541e38e3280badea4c4c0ca9339cabdc3b39 |
| askey-wilson-apery-1984.pdf | https://doi.org/10.1017/S1446788700024691 | 26475a843753b6c0121c6f80338ff786e85394698d8a92b4aa60dcaa3f581ca4 |
| cohen-apery-generalisation-1981.pdf | https://www.numdam.org/article/BSMF_1981__109__269_0.pdf | cffe9ae687d12e584f4359188cff4d0576c8f65af81ea81c08d3bf27b12be85c |
| cohen-parametric-cf-2304.11727.pdf | https://arxiv.org/pdf/2304.11727 | b3ea1b2c159728c4ca3ab5fde041a8542ff828dd18aa6599f1e466c35c94d143 |

Rivoal's paper appeared in Experimental Mathematics 13(4) (2004), 503–508.
The Markov–WZ paper appeared in Discrete Mathematics and Theoretical Computer
Science 10:3 (2008), 115–124. Tauraso's retrieved manuscript is arXiv
2511.18034v1, dated 22 November 2025. Cohen's older paper is Bulletin de la
Société Mathématique de France 109 (1981), 269–281.

Priority status is now precise: exact published matches have been found for
Z_sigma and its accelerated Koecher–Leshchiner expansion, and for A_n with its
entire three-term recurrence. The finite-balance exposition is an independent
derivation. Other identifications require a specific proof or source, not an
inference from unsuccessful searching.

