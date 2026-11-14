# Primary-source context for the truncated rank-one zeta Hankel family

Status: bounded primary-source comparison, after the independent construction.
This identifies classical mechanisms and one exact Apéry/Wilson identity. It
does not establish that the precise truncated family is new or that its
conjectured adjacent-Apéry formula holds for all indices.

The construction being compared is
[rank-one-zeta-subspace-approximants.md](rank-one-zeta-subspace-approximants.md):
for s=2m+1, h=n+1 and K=2n+m+1, the matrix is
[μ_s(t^(i+j)/D_K)]_(0≤i,j<h), D_K=∏_(j=1)^K(t+j²).
Its single zeta-dependent entry is the bottom-right one. The rational root
q_(h,s) has error equal to the minimum monic degree-n squared norm for
μ_s/D_K. This is distinct from the paired Mellin–Rodrigues construction and
from the earlier fixed one-pole Wilson endpoint sequences.

## 1. Exact classical framework: rational Gaussian quadrature

Walter Gautschi, [*Gauss-type Quadrature Rules for Rational Functions*,
arXiv:math/9307223](https://arxiv.org/abs/math/9307223), Theorem1.1,
pp2–3 of the preprint, constructs rational Gaussian quadrature by dividing a
measure by its prescribed pole polynomial, applying ordinary Gaussian
quadrature, and multiplying the weights back. Its remark credits the case
of 2n distinct real simple poles to Van Assche–Vanherwegen,
[*Quadrature formulas based on rational interpolation*](https://arxiv.org/abs/math/9307221),
Math. Comp.61(1993),765–783. Gautschi's full theorem is the directly inspected
primary source here; the latter attribution is recorded as his attribution.

The exact parameter dictionary in our notation is

    dλ=dρ₀=dμ_s/D_(m+1),
    Ω_n=D_(2n+m+1)/D_(m+1),
    source ζ_ν=1/j²,  j=m+2,...,m+2n+1,
    source pole multiplicities all1, source m=2n.

The source ω_(2n)=∏(1+t/j²) differs from Ω_n only by a positive constant.
Thus its orthogonal polynomial is exactly our monic Schur minimizer p_n.
There is no negative-parameter analytic continuation in this identification.

For clarity, the associated Markov function in the project is
S(z)=∫dρ₀(t)/(z−t). All 2n samples S(−j²) are rational: partial fractions
of 1/[D_(m+1)(t)(t+j²)] have m+2 nodes, while the zeta coefficient is the
Lagrange sum of a degree-m polynomial. That sum is zero. The project's
multipoint-Padé proof makes this cancellation explicit. Gautschi supplies
the general quadrature mechanism, not a theorem about this specific zeta
measure or its rational samples.

The author preprint of Bultheel–González-Vera–Hendriksen–Njåstad,
[*Quadrature and orthogonal rational functions*](https://nalag.cs.kuleuven.be/papers/ade/qorf/preprint.pdf),
§8, Theorem8.1, also states this relation for varying measures, with the
unbounded-support case explicitly included. Published JCAM127(2001),67–91.
Its author PDF was inspected through the web tool; a local download returned
HTTP503, so no local file/hash is claimed for that paper.

## 2. Exact determinant umbrella, without a claimed closed evaluation

Christian Krattenthaler,
[*A determinant identity for moments of orthogonal polynomials that implies
Uvarov's formula for the orthogonal polynomials of rationally related
densities*](https://arxiv.org/abs/2103.03969), Theorem1, gives determinants of
moments after rational multiplication of a measure. Our substitution is:
source matrix size n=h; number of numerator factors m=0; number of denominator
factors k=K; poles y_j=−j²; base measure μ_s. Since h<K, it is explicitly
**equation(1.3)**, not equation(1.2), that applies. It expresses the Hankel
determinant using ordinary orthogonal polynomials and their Cauchy transforms.

This is a literal embedding of the determinant in a published general
identity. It does not by itself evaluate the rational root, reduce its
arithmetic cost, or prove an Apéry recurrence. The source's theorem is valid
for general linear functionals; our positive measure satisfies its ordinary
analytic hypotheses.

## 3. A close Hankel approximation precedent, but a different family

Timothy Ferguson,
[*Rational Approximations via Hankel Determinants*](https://arxiv.org/abs/2003.10616),
Theorem1 and its appendix, use positivity and Hilbert-space projections to
produce monotone determinant ratios approximating an unknown zeroth moment.
His zeta example (§2, pp3–4) uses the fixed [0,1] measure

    dρ(y)=(-log(1-y))^(s-1)/[(s-1)! y] dy.

Its positive-degree moments are finite rational sums. This is a clear
precedent for the Schur/projection mechanism, but not our changing measure
μ_s/D_K or bottom-right degree drop. The distinction can be checked already
at the first ζ(3) approximant: Ferguson's is 8/7, whereas ours is5/4.
No convergence-rate or denominator theorem is transferred between them.

## 4. Prévost proves the all-half Wilson–Apéry identity explicitly

Marc Prévost, [*A new proof of the irrationality of ζ(2) and ζ(3) using Padé
approximants*](https://www.math.ru.nl/~zudilin/zw/Prevost96.pdf), JCAM67(1996),
219–235, DOI10.1016/0377-0427(95)00019-4. Original scanned pp229–231 were
visually inspected, not inferred from an abstract or unreliable OCR formulas.

On p229, equations(53)–(54) give

    Π_(2n)(y) = 4F3(-n,n+1,y+1,-y;1,1,1;1).

The same page identifies these as Wilson polynomials with all four parameters
1/2, through y=it−1/2. On p230, equation(57) explicitly states
Π_(2n)(n)=A_n. Equation(58) gives
B_n=H_n^(3)Π_(2n)(n)+Θ_(2n)(n)/2, using his associated polynomial Θ.
His names for our A_n,B_n are b_n,a_n. This is exact prior attribution for
the all-half Wilson evaluation; it is not a conjectural numerical analogy.
The truncation's proposed formula involving 5A_h−A_(h−1) does not follow
from these cited equations without further work.

## 5. The negative-parameter Wilson expression is not an orthogonality theorem

For s=3, direct gamma recurrence gives

    w₃(y)/D_(2h)(y²)
      =(1/(4π)) |Γ(1+iy)^3 Γ(-2h+iy)/Γ(2iy)|².

Thus the density has the formal Wilson parameters(1,1,1,−2h). This is
outside the positive-parameter orthogonality theorem used elsewhere in the
project. Simply substituting these parameters in the terminating Wilson
polynomial gives the wrong Schur polynomial: for h=2 its degree-one monic
polynomial is t−11, while the actual minimizer is t−1/5. The exact checker
confirms both. Analytic continuation would need its contour/residue or
additional-mass analysis; this gamma identity alone cannot supply it.

## 6. What is established, and what remains an attribution question

- The varying-weight rational Gaussian/Padé framework is classical and has
  the exact source dictionary in §1.
- The rationally modified Hankel determinant is covered by Krattenthaler's
  explicit n<k identity, with the dictionary in §2.
- The all-half Wilson evaluation of Apéry's A_n and associated numerator is
  explicitly in Prévost1996.
- The exact all-h formula q_(h,3)=(5B_h−B_(h−1))/(5A_h−A_(h−1)) remains a
  conjecture in the current construction note, checked through h=12. The
  present source work neither proves it nor identifies an exact prior
  occurrence of this transformed approximant.

Searches of the first fractions did not identify a primary occurrence.
That search outcome provides no evidence of historical novelty. The useful
next step is an exact orthogonality/Padé proof of the adjacent-Apéry identity;
only then can its identity with a named historical family be assessed reliably.

## Reproducibility

Primary PDFs and extracted text are in `../sources/`; URLs, byte counts and
SHA256 hashes are recorded in `rank-one-primary-context-provenance.json`.
Prévost is a scanned PDF: plain pdftotext is empty, and any `.ocr.txt` is a
search aid only; equations(53)–(58) above were read from the page images.
The comparisons are checked by
`python codex/computations/check_rank_one_source_dictionary.py`.


## 7. The mixed numerator identity now has an elementary proof

Use Π_h(y) as shorthand for Prévost's Π_(2h)(y), and define his associated
polynomial directly, without OCR of a long binomial numerator:

\[
 \Theta_h(y)=\mathcal L_x\!\left(
          \frac{\Pi_h(x)-\Pi_h(y)}{x-y}\right),\qquad
 \mathcal L(1)=0,\quad\mathcal L(x^j)=jB_{j-1}\ (j\ge1),
\]

where B₁=−1/2. Put C_h(y)=H_y^(3)Π_h(y)+Θ_h(y)/2 at nonnegative integer y.
Prévost's exact formula gives C_h(h)=B_h. The all-index proof in
[apery-mixed-numerator-balance.md](apery-mixed-numerator-balance.md) now proves

\[
 C_h(h-1)=\frac{B_h+B_{h-1}}6.
\]

Together with the independently proved denominator identity
Π_h(h−1)=(A_h+A_(h−1))/6, this proves

\[
 \frac{C_h(h)-C_h(h-1)}{\Pi_h(h)-\Pi_h(h-1)}
 =\frac{5B_h-B_{h-1}}{5A_h-A_{h-1}}.
\]

Identifying this quotient with the truncated Schur root still requires a
separate orthogonality or interpolation argument. The mixed-numerator theorem
alone does not resolve that remaining conjecture. The source-dictionary checker
records its original finite checks; the linked proof has its own symbolic
rational certificates and independent pole-jet checks. Independent review of
the new proof is pending.
