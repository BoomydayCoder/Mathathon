# Independent review: Bernoulli repairs, Stieltjes positivity and Wilson polynomials

**Verdict: accepted.** The analytic, algebraic and classical-special-function
claims in [the construction](bernoulli-wilson-stieltjes.md) check out.
This review covers its completed Wilson section, not only the earlier
Bernoulli construction. No irrationality result or historical novelty is
inferred. The all-degree Wilson orthogonality and norm are explicitly
imported classical results; the finite Apéry dictionary is an elementary
identity proved in the note.

An independent Fraction-only checker is
[check_bernoulli_wilson_independent.py](../computations/check_bernoulli_wilson_independent.py).
It does not import the author's code or SymPy.

## 1. Algebra and the sign of the jump

For a monomial, direct binomial expansion gives
\[
 \mathcal D(t^k)=
 \sum_{\substack{1\le \ell<2k\\\ell\ {\rm odd}}}
 \binom{2k}{\ell}\frac{(-1)^{(\ell-1)/2}}{2^\ell}
 t^{\,k-(\ell+1)/2}.
\]
Its leading coefficient is \(k\). Thus the zero-constant inverse of
\(\mathcal D\) exists and is unique on every finite polynomial degree.
I constructed this inverse by triangular elimination, independently of
the Bernoulli formula, and checked the first fifteen resulting repairs.

In the Bernoulli formula, the two shifted arguments really are \(1-ix\)
and \(-ix\), where \(t=x^2\). The Bernoulli translation identity therefore
gives
\[
 R_k((x+i/2)^2)-R_k((x-i/2)^2)=2ikx^{2k-1}.
\]
The sign \((-1)^k\), the factor \(2i\), and the denominator \(x\) all
agree. Reflection removes odd powers of \(x\), and expansion at \(1/2\)
gives the stated coefficients and examples.

There is no need to assume Euler's evaluations of even zeta values in
order to obtain their signs: first define \(R_k\) by the Bernoulli
formula without any assertion about signs. The positive integral below
then implies \(R_k'(0)>0\), hence
\((-1)^jB_{2j}(1/2)>0\). The author has added this independent sign
argument.

For \(0<v<1\), the upper boundary value of \(F\) has negative imaginary
part:
\[
 F_+(-v)=(\sqrt{G(v)}-i/2)^2.
\]
Consequently its polynomial composition has jump
\[
 P(F_+(-v))-P(F_-(-v))
 =-2i\sqrt{G(v)}\,(\mathcal DP)(G(v)).
\]
A transform \(\int f(v)/(z+v)\,dv\) has upper-minus-lower jump
\(-2\pi i f(v)\). This independently fixes the coefficient \(1/\pi\)
in the claimed representation, and \(k/\pi\) for \(R_k\). In particular
the sign would be wrong if the upper and lower boundary values were
interchanged.

## 2. Analytic continuation and the entire positivity cone

The real integral for \(F\), obtained by integration by parts and
differentiation, supplies its single-valued continuation to the slit
plane. It uses neither a zeta evaluation nor Wilson orthogonality.
For a fixed polynomial \(P\) with \(P(0)=0\), \(P(F(z))=O(z^{-1})\)
at infinity. Near zero the continued logarithm gives a bound by a fixed
power of \(1+|\log|z||\); near \(-1\) the function is bounded with a
square-root correction. These estimates are uniform on the small
endpoint arcs in the slit plane. Their lengths are \(O(\epsilon)\),
so the endpoint contour contributions vanish.

The resulting density is integrable: it has at most a power of a
logarithm at zero and is \(O(\sqrt{1-v})\) at one. There are no omitted
endpoint atoms or residual poles. The contour representation is valid
for every polynomial in the claimed class.

Uniqueness of the finite signed measure from its Cauchy transform gives
both directions of
\[
 P(F)\text{ positive Stieltjes on }[0,1]
 \quad\Longleftrightarrow\quad \mathcal DP\ge0\text{ on }[0,\infty).
\]
The fact that \(G:(0,1)\to(0,\infty)\) is onto is essential for the
forward implication. Continuity handles \(t=0\).

The factorization proof that every real nonnegative polynomial on the
half-line is \(A^2+tB^2\) is correct, including positive real roots,
negative real roots, conjugate pairs and zero roots. The zero polynomial
is handled by \(A=B=0\). The displayed multiplicative identity preserves
this form. Over rational coefficients the note claims only a sufficient
construction, which is appropriate: it does not infer a rational
two-square representation from a real one.

For a nonzero \(Q\ge0\), its density is positive except at finitely many
interior points. This gives strictly positive Gram matrices and the usual
\([N-1/N]\) Padé approximants at infinity, with denominator roots on
\((-1,0)\). I suggested specifying this Padé index explicitly to avoid
an unnecessarily broad reading.

## 3. Moments and the noncircular recovery of Euler's formula

In the substitution \(v=\operatorname{sech}^2r\), the density for
\(R_k(F)\) becomes
\[
 \frac{2k}{\pi^{2k}}r^{2k-1}\operatorname{sech}^2r\tanh r\,dr.
\]
Multiplication by \(v^m\), followed by one integration by parts, gives
exactly the coefficient
\(k(2k-1)/((m+1)\pi^{2k})\) in equation (13). The \(r=0\) boundary
vanishes also for \(k=1\).

The base moments \(a_m\) follow from a beta integral. Composition of
its convergent expansion with the finite polynomial \(R_k\) gives the
finite convolution formula (16), with sign \((-1)^{j-1}\).
The measures on \([0,1]\) generally have moments in
\(\mathbb Q[\pi^{-2}]\), not rational moments. For example the
normalized first moment at \(k=2\) is \(1/3-2/\pi^2\).

The different measure in the design variable is
\[
 d\nu(t)=\pi^2\operatorname{sech}^2(\pi\sqrt t)
                         \tanh(\pi\sqrt t)\,dt.
\]
Taking the \(1/z\) coefficient gives
\[
 P'(0)=\pi\int_0^1\sqrt{G(v)}Q(G(v))\,dv
       =\int_0^\infty Q(t)\,d\nu(t).
\]
Thus this is a probability measure with the claimed rational moments
\((2j+1)|B_{2j}(1/2)|\). The two moment measures and their arithmetic
properties are correctly distinguished.

For the quadratic repair, integration by parts in the second-derivative
identity for \(\operatorname{sech}^{q}\) gives
\[
 \frac{V_{m+1}}{U_{m+1}}-\frac{V_m}{U_m}
             =-\frac1{2(m+1)^2}.
\]
The initial value \(V_0=\pi^2/12\) follows from the proved Stieltjes
representation of the explicit polynomial \(R_2=t^2+t/2\); no even-zeta
formula is required. This proves both the Euler-tail formula and its
rational convolution counterpart.

The proposed concentration argument is valid. For instance
\(\log\cosh r\le r^2/2\) on the real axis yields
\(U_m\ge e^{-1}/\sqrt{m+1}\) by integrating up to
\(1/\sqrt{m+1}\). The fixed-\(\delta\) tail bound then forces
\(V_m/U_m\to0\). Hence the recovery of \(\zeta(2)=\pi^2/6\) is
noncircular. The author's numerical moment checker uses known even-zeta
values for some auxiliary initial data; that is corroboration only, not
the logical basis of this noncircular proof.

## 4. Exact Wilson measure, norms and the Apéry evaluation

Gamma reflection gives
\[
 \frac1{2\pi}
 \left|\frac{\Gamma(1/2+iy)^4}{\Gamma(2iy)}\right|^2
 =2\pi^2y\,\operatorname{sech}^2(\pi y)\tanh(\pi y).
\]
The right side is precisely \(d\nu(y^2)/dy\). The factor \(1/(2\pi)\)
is necessary.

I independently checked DLMF §18.25 through the web tool
(retrieval reference **turn700929view0**, 2026-11-14), including the
parameter domain, the variable \(t=y^2\), equations 18.25.2–18.25.5,
and the references to Wilson (1980). Direct urllib requests for the
individual TeX endpoints returned HTTP 403; the subsequent web retrieval
did succeed. The exact formulas are also transcribed in
[the saved retrieval note](../sources/dlmf-wilson-norm-retrieval.md).

For four half-parameters, the six pair sums are all \(1\), their total
sum is \(2\), and the classical squared norm is
\(2\pi(n!)^6/(2n+1)\) in the unnormalized measure. Dividing the measure
by \(2\pi\) gives equation (26). The leading coefficient of \(W_n\)
is \((-1)^n(2n)!/n!\), so its monic norm and the Hankel determinant
product (27) follow with exactly the stated factorial powers.

The finite evaluation needs no orthogonality theorem. At
\(t=-(m+1/2)^2\), each product in the terminating Wilson sum becomes
\((-m)_k(m+1)_k\). Combining both pairs of Pochhammer symbols gives
the product of four binomial coefficients defining \(K(n,m)\).
The \(n=m\) specialization is exactly \(A_n\).

The evaluation point moves with the degree and is negative, whereas the
orthogonality support is positive. The note correctly warns that an
ordinary fixed-variable three-term recurrence does not automatically
give the diagonal Apéry recurrence.

The classical polynomial family and norm are not novel. The exact
Wilson bibliographic attribution now uses DOI **10.1137/0511064**.
The original 1980 paper itself was not retrieved in this audit. The
locally retrieved Cheng–Chiang paper explicitly supplies the operator,
hypergeometric definition, weight and attribution; its PDF hash matches
the saved provenance:
eb32ab545612543d38f98afb2683f42b9b56ddbb684c198ed0f63fba7d00171c.

## 5. Reproducible checks and limits

The independent checker passes:

- 15 triangular Wilson inversions, compared with independently generated
  Bernoulli numbers;
- 16 rational square-cone inverses and 15 rational design moments;
- 66 exact Wilson Gram entries, including all norms through degree 10;
- 143 finite-grid evaluations and 11 diagonal Apéry evaluations;
- 10 direct Hankel determinants and 150 rational Euler convolutions.

The author's original checker also passed, including its direct complex
Stieltjes quadratures. These checks support the written proofs and catch
normalization errors; they do not establish analytic continuation,
measure uniqueness, all-degree classical orthogonality, or arithmetic
approximation bounds by finite extrapolation.

No substantive mathematical gap was found. The established result is a
precise positive-function construction and a classical special-function
dictionary. Denominator growth, useful zeta linear forms, and any claim
of historical novelty remain separate questions.
