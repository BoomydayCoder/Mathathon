# Independent audit: Newton multiplication with a Hilbert or beta Gram matrix

**Status: accepted.** This reviews
[newton-hilbert-selfadjoint-classification.md](newton-hilbert-selfadjoint-classification.md),
its checker, and the cited Ismail–Koelink primary source. The all-index
classification and positive-functional converse are correct. The exact
classical operator match is also correct, including its sign, scale, spectral
shift, and both positive Wilson specializations.

Algebraic symmetry of multiplication, positivity of a moment functional, and
essential selfadjointness of a closed operator are distinct assertions.
The note correctly proves the first two directly and imports the last,
where used, from the stated classical spectral theorem. No historical
novelty conclusion is part of this audit.

## 1. Necessity uses two anchor columns; sufficiency checks every pair

For a Newton basis satisfying
\[
 tQ_j=a_jQ_{j+1}-\lambda_jQ_j,\qquad Q_0=1,\quad a_j\ne0,
\]
and beta moments
\[
 h_r=\kappa\mathrm B(r+A,B-A),\qquad A>0,\quad B>A,
\]
the ratio \(h_{r+1}/h_r=(r+A)/(r+B)\) is exact. Thus symmetry is equivalent to
\[
 (a_j-a_k)(j+k+A)=(\lambda_j-\lambda_k)(j+k+B)
\]
for every \(j,k\ge0\). All moments used here are strictly positive.

Write \(d=a_j-a_0\) and \(c=(a_1-a_0)/(B+1)\).
The \(k=0\) equation gives
\(\lambda_j-\lambda_0=d(j+A)/(j+B)\), while the case \(j=1\) gives
\(\lambda_1-\lambda_0=c(A+1)\). Substitution into the \(k=1\) equation
reduces it exactly to
\[
 \frac{B-A}{j+B}\,[d-cj(j+B)]=0.
\]
Since \(B>A>0\), neither the prefactor nor its denominator vanishes.
There is no division by \(j\) or \(j-1\); the indices zero and one are
covered by their definitions. Therefore, for all indices,
\[
 a_j=a_0+cj(j+B),\qquad
 \lambda_j=\lambda_0+cj(j+A).
\]

Conversely the two differences factor as
\[
 a_j-a_k=c(j-k)(j+k+B),\qquad
 \lambda_j-\lambda_k=c(j-k)(j+k+A),
\]
which proves the original symmetry condition for every pair, including the
diagonal. Necessity from the anchor equations alone would not suffice
without this converse; the note supplies both.
The uniform Hilbert matrix is \(A=1,B=2,\kappa=1\), giving the claimed
\(j(j+2)\) and \(j(j+1)\) formulas.

## 2. The converse really is a positive moment functional

For any allowed real constants, the recursive polynomials have successive
degrees and nonzero leading coefficients. Thus \(UQ_j=x^j\) defines a
linear bijection of real polynomial spaces. Transporting the positive beta
inner product through \(U\) gives a strictly positive definite form.
Its multiplication-by-\(t\) operator is symmetric by the all-pair calculation.

Repeatedly transferring a factor of \(t\) between the arguments shows
\[
 \langle t^r,t^s\rangle=\langle t^{r+s},1\rangle.
\]
With \(L(p)=\langle p,1\rangle\), bilinearity now gives
\(L(pq)=\langle p,q\rangle\) for arbitrary polynomials. This step is
essential: an arbitrary positive transported inner product need not have
the moment-product property. Here symmetry proves it.
Consequently \(L(p^2)>0\) for every nonzero polynomial \(p\), and the
Hamburger moment theorem applies.

No assumption \(c>0\) is needed for this argument. At \(c=0\), the explicit
measure is the pushforward of the beta measure under
\(t=a_0x-\lambda_0\). In the uniform case its interval endpoints are
\(-\lambda_0\) and \(a_0-\lambda_0\), even when \(a_0<0\).
It has nonnegative support exactly when
\(\lambda_0\le\min(0,a_0)\).
If some \(a_j=0\), the full polynomial-basis construction fails; no such
degenerate case is included in the theorem.

The generic proof makes no assertion of nonnegative support or uniqueness
of the representing measure. Those require additional input.

## 3. Differential operator, tridiagonal coefficients, and endpoints

The weighted differential coefficient is
\[
 \frac{[x^2(1-x)w(x)]'}{w(x)}
 =(A+1)x-(B+1)x^2,\qquad
 w=x^{A-1}(1-x)^{B-A-1}.
\]
Acting on \(x^j\) therefore gives exactly
\(a_jx^{j+1}-\lambda_jx^j\), with the claimed signs.
The integration-by-parts boundary factor is
\[
 x^2(1-x)w=x^{A+1}(1-x)^{B-A},
\]
which vanishes at both endpoints for polynomials under the stated strict
inequalities. The derivative integrals are also finite.

For the Jacobi operator \(J=-w^{-1}D[x(1-x)wD]\), direct product
differentiation gives
\[
 T=\frac c2(xJ+Jx)+(a_0-cB/2)x+cA/2-\lambda_0.
\]
The weighted square-completion formula is correct. Its cross-term boundary
factor is \(x(1-x)w=x^A(1-x)^{B-A}\), which also vanishes.
The apparent \(f/x\) singularity in the square is integrable because its
squared contribution is bounded by a constant times \(x^{A-1}\) near zero.
No unstated condition \(A\ge1\) is needed.

In the uniform case, multiplication by \(x\) in the basis
\(L_n=P_n(1-2x)\) has coefficients
\[
 -\frac{n+1}{2(2n+1)},\quad \frac12,\quad
 -\frac{n}{2(2n+1)}.
\]
The anticommutator multiplies the upper coefficient by
\(a_0+cn(n+2)=a_n\), and the lower by \(a_{n-1}\).
Its diagonal is \(cn(n+1)/2+a_0/2-\lambda_0\).
This proves the full tridiagonal formula, including \(n=0\), where the
backward term is absent.

Requiring the additional binomial evaluation at \(k=n=1\) gives
\(-2c/a_0=-2\), hence \(a_0=c\ne0\). The factored product then proves every
evaluation, also \(k>n\), when both sides vanish. It leaves \(\lambda_0\)
free. Choosing \(c=1,\lambda_0=1/4\) gives the half-square products,
and the linear Legendre coefficients times \(a_n=(n+1)^2\) give the cubes.

Evaluation at a negative half-square produces the mixed binomial kernel
and its stated recurrence. The adjacent balance is an additional proved
identity; only after using it do the coefficients \(17\) and \(5\) follow.
They are not forced by Gram symmetry alone.

## 4. Wilson series and the beta(2,2) normalization

For \(c\ne0\), choose \(\eta\) with
\(\eta(B-\eta)=a_0/c\). Factoring the Newton products gives
\[
 Q_j(t)=
 \frac{(A/2+iY)_j(A/2-iY)_j}{(\eta)_j(B-\eta)_j},
 \qquad
 Y^2=\frac{t+\lambda_0-cA^2/4}{c}.
\]
Both numerator and denominator reproduce the required quadratic factors.
Nonvanishing of every \(a_j\) ensures the denominator factors in this
terminating formula never vanish.
Transporting
\({}_2F_1(-n,n+B-1;A;x)\) therefore gives exactly the displayed terminating
\({}_4F_3\), with Wilson parameters
\[
 A/2,\quad A/2,\quad \eta-A/2,\quad B-\eta-A/2.
\]
Their sum is \(B\), as required by the upper parameter \(n+B-1\).
This algebraic identification does not assume a positive real Wilson
integral for arbitrary parameter choices.

At \(A=2,B=4,a_0=4,c=1,\lambda_0=1\), take \(\eta=2\).
All four parameters are one, and \(Q_k=D_k/((k+1)!)^2\).
The conventional Wilson measure with factor \(1/(2\pi)\) is twice
\(t\,d\mu_3\), so the note correctly uses \(\kappa=1/2\).
Its mass is \(1/12\), matching
\(\tfrac12\mathrm B(2,2)\). The spectral identification and this mass
fix the functional; the claimed all-index Gram equality is not being
deduced from the finite checker alone.

## 5. Direct comparison with the primary operator and spectral edge

I read Ismail–Koelink, arXiv:1108.5716, §3, equation (3.1), the Wilson
parameter display immediately following it, and Theorem 3.1 in the archived
primary text.
Under
\[
 \xi=1-2x,\quad
 \alpha=A-1,\quad\beta=B-A-1,\quad\delta=\eta-A,
\]
one has \(D_\xi=-D_x/2\) and \(D_\xi^2=D_x^2/4\).
Substitution into the actual source operator gives
\[
 T_{\rm paper}
 =2\{x^2(1-x)D_x^2+
 [(A+1)x-(B+1)x^2]D_x-\eta(B-\eta)x\}.
\]
Therefore
\[
 T_{\rm ours}=-\frac c2T_{\rm paper}-\lambda_0
\]
with exactly the stated negative scale.
The source Wilson parameters become the same multiset as above.
Its eigenvalue relation
\(\lambda_{\rm paper}=-A^2/2-2Y^2\) becomes
\[
 t=cY^2+cA^2/4-\lambda_0.
\]

In both selected cases, \(\alpha,\beta>-1\) and the source's
\(\gamma=-\eta(B-\eta)\) is real, so Theorem 3.1 applies.
The uniform case has \(\alpha=\beta=\delta=0\); the companion has
\(\alpha=\beta=1,\delta=0\).
All their Wilson parameters are positive, so no discrete extra masses occur.
The gamma identities in the all-half case give the probability density
\[
 2\pi^2y\operatorname{sech}^2(\pi y)\tanh(\pi y)\,dy.
\]
For \(a_0=c>0\), its affine spectral image has lower edge
\(c/4-\lambda_0\). Essential selfadjointness in this selected classical
family supplies the corresponding unique spectral realization.
Thus \(\lambda_0=c/4\) is indeed the maximal shift compatible with
nonnegative support in that family. The elementary symmetry equations
alone do not select this shift.

## 6. Independent checks and scope

The author's checker passed in full. A second
[independent checker](../computations/check_newton_hilbert_independent.py)
solves the finite symmetry systems from scratch for three beta choices,
recovering precisely the three-dimensional solution space.
For five parameter choices, including \(c=0\), a negative constant scale,
and \(c<0\), it reconstructs moments only from \(L(Q_j)=h_j\), then
checks the full product Gram identities and positive Hankel determinants.
It also starts from the untransformed primary-source differential
coefficients and verifies the coordinate and spectral maps. All checks passed.

These finite tests support the all-index proofs above. No general support
classification, arbitrary Wilson positivity claim, new spectral-family
claim, or identification of the separate \(D_{2h}\) optimized norm follows.
The accepted result is the exact classification within the specified
Newton/Gram design, together with its correctly attributed classical
spectral interpretation.
