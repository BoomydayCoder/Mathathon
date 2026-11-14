# Independent pencil review and stronger enclosures from positivity on a half-line

**Status.** The proof in `real-rooted-zeta-pencils.md` has been independently reviewed from its matrix identities, not merely its finite tests; no gap was found. This note develops further consequences of standard positive-matrix and Stieltjes-moment theory. Historical novelty is not claimed. There are no claims about favorable denominators, irrationality, or an effective convergence rate.

## 1. Audit of the original proof

Let \(s\ge3\) be odd, \(K>N\ge0\), \(r\ge2\) even, \(h=K-N\), and \(\xi=\zeta(s)\). The original matrix is
\[
 G(X)_{ab}=\mu_{s,X}\left(\frac{D_N(t)^rt^{a+b}}{D_K(t)}\right),
 \quad 0\le a,b<h.
\]
The positive moment representation makes \(G(\xi)\) strictly positive definite: a nonzero polynomial has only finitely many zeros, while the measure has positive density on the whole positive half-line. All moments converge because the density decays exponentially in \(\sqrt t\).

Writing \(C=[X]G\), the congruence
\[
 G(\xi)^{-1/2}G(X)G(\xi)^{-1/2}=I+(X-\xi)H,
 \quad H=G(\xi)^{-1/2}CG(\xi)^{-1/2}
\]
is valid. The matrix \(H\) is real symmetric. Its nonzero eigenvalues give all finite roots, with multiplicities.

The residue decomposition
\[
 C=V\operatorname{diag}(c_j)V^T,
 \quad V_{aj}=(-j^2)^a,
 \quad c_j=\frac{j^{s-1}D_N(-j^2)^{r-1}}
 {(D_K/D_N)'(-j^2)}
\]
has invertible square \(V\), nonzero \(c_j\), and
\(\operatorname{sign}c_j=(-1)^{j-1}\). The denominator has exactly \(j-N-1\) negative factors and the numerator has sign \((-1)^{N(r-1)}\); since \(r-1\) is odd, their quotient has the stated sign. Therefore the degree is exactly \(h\), every root is real, and the inertia gives the claimed numbers on either side of \(\xi\). For the preprint parameters, the count \(\lfloor37n/2\rfloor\) below is correct for both parities of \(n\).

The positive-definite region is precisely the interval between the nearest roots on either side, not an arbitrary interval on which the determinant is positive. Repeated roots do not invalidate the statement. When only one parity occurs (necessarily \(h=1\)), the region is a half-line instead. For a principal submatrix with singular \(C\), the same argument gives real finite roots and degree equal to \(\operatorname{rank}C\); no assertion of full dimension is warranted there.

The displayed two-by-two determinants and algebraic endpoints in the original note check out exactly.

## 2. Add the localizing matrix for \(t\ge0\)

The integral support is the positive half-line. Consequently we have not just positivity of \(q(t)^2\), but also positivity of \(tq(t)^2\). Define
\[
 G^{(1)}(X)_{ab}=\mu_{s,X}\left(\frac{tD_N(t)^rt^{a+b}}{D_K(t)}\right).
\]
It too is strictly positive definite at \(\xi\). Its coefficient of \(X\) is
\[
 C^{(1)}=V\operatorname{diag}(-j^2c_j)V^T.
\]
All residue signs are reversed. Thus its determinant has exactly \(h\) real roots, with the number of roots below \(\xi\) equal to the number of **even** integers in \(\{N+1,\ldots,K\}\). The number above equals the number of odd integers.

Intersecting the positive-definite intervals of \(G\) and \(G^{(1)}\) combines both constraints. It improves either side when the other matrix supplies a nearer endpoint.

### Exact two-by-two improvements

For \(N=0,K=2\), the additional determinant for \(s=3\) is
\[
 \det G^{(1)}(X)=-\frac{768X^2-1904X+1179}{48}.
\]
The new lower endpoint and the original upper endpoint give
\[
 \boxed{\frac{119-\sqrt{13}}{96}<\zeta(3)<\frac{17+\sqrt5}{16}.}
\]
The interval is approximately \((1.2020255075,1.2022542486)\).

For \(s=5\),
\[
 \det G^{(1)}(X)=-\frac{18432X^2-42396X+24143}{288},
\]
so
\[
 \boxed{\frac{3533-\sqrt{120873}}{3072}<\zeta(5)
 <\frac{317+\sqrt{6601}}{384}.}
\]
The interval is approximately \((1.0368919450,1.0371003604)\).

### Three-by-three exact algebraic enclosures

For \(N=0,K=3\), the middle root of \(\det G\) is below the target and the middle root of \(\det G^{(1)}\) is above it. Exact rational root isolation gives the following weaker, simpler rational consequences:
\[
 \boxed{1.202056695706<\zeta(3)<1.202056932326,}
\]
\[
 \boxed{1.036927616686<\zeta(5)<1.036927779168.}
\]
These are certified by exact polynomial arithmetic, not by comparison against numerical values of zeta. The script accompanying this note constructs the rational matrices, computes the determinant polynomials, and isolates their roots with rational intervals.

## 3. Joint intervals are nested even when \(K\) changes

Fix \(N,r,s\), let \(K\) increase, and retain matrix size \(h=K-N\). Denote the common positive-definite interval by \(I_K\). Then
\[
 \boxed{I_{K+1}\subset I_K.}
\]
This is stronger than the nesting of principal submatrices, because the weight changes with \(K\).

To prove it, write \(a=K+1\). For any nonzero \(q\) of degree less than \(h\),
\[
 \frac{q^2}{D_K}=\frac{a^2q^2+tq^2}{D_{K+1}},
\]
\[
 \frac{tq^2}{D_K}=\frac{(tq)^2+a^2tq^2}{D_{K+1}}.
\]
Multiply these identities by the common numerator \(D_N^r\) and apply \(\mu_{s,X}\). Positivity of both new matrices implies positivity of both old quadratic forms. The polynomial \(tq\) has degree at most \(h\), so it fits in the new dimension \(h+1\). This proves the inclusion. The same argument proves inclusion for positive-semidefinite regions.

This section alone does not prove that these particular intervals shrink to a point. The following fixed-denominator version does establish convergence.

## 4. Fixed denominator, increasing dimension: convergence theorem

Fix \(s,N,K,r\), but now let the matrix size \(m\) increase independently of \(h=K-N\). Define
\[
 G_m^{(\ell)}(X)_{ab}=\mu_{s,X}\left(
 \frac{t^{a+b+\ell}D_N(t)^r}{D_K(t)}\right),
 \quad0\le a,b<m,\quad\ell=0,1.
\]
Let
\[
 J_m=\{X\in\mathbb R:G_m^{(0)}(X)\succeq0,
                         G_m^{(1)}(X)\succeq0\}.
\]

**Theorem.** For \(m\ge h\), \(J_m=[L_m,U_m]\) is a compact interval with algebraic endpoints and
\[
 L_m<\xi<U_m,\qquad
 L_m\nearrow\xi,\qquad U_m\searrow\xi.
\]
Monotonicity need not be strict at each step. Both endpoints are roots of one of the two rational determinant polynomials. Each determinant has exactly \(h\) finite real roots, including multiplicities, even when \(m>h\).

**Proof of finite-dimensional statements.** The residue Vandermonde now has \(m\) rows and full column rank \(h\). Its coefficient matrix has exactly the nonzero inertia of \(\operatorname{diag}(c_j)\), plus \(m-h\) zero eigenvalues. The positive-metric argument therefore gives degree exactly \(h\), all roots real, and the same sign counts as before; the shifted pencil reverses them. Each PSD region is a closed interval or half-line. Together the two inertia patterns supply a finite endpoint on each side. Their intersection contains \(\xi\) in its interior. Principal-submatrix inclusion makes \(J_{m+1}\subset J_m\).

**Proof that the intersection is a singleton.** Put
\[
 D(t)=D_K(t)/D_N(t),\qquad W(t)=D_N(t)^{r-1},
 \qquad\mathcal L_X(q)=\mu_{s,X}(Wq/D).
\]
If \(X\) belongs to every \(J_m\), then
\(\mathcal L_X(q^2)\ge0\) and \(\mathcal L_X(tq^2)\ge0\) for every real polynomial \(q\). The Stieltjes moment-existence theorem gives a positive measure \(\nu\) on \([0,\infty)\) representing \(\mathcal L_X\). The exact theorem used is: a real sequence is the moment sequence of such a measure if and only if both all its Hankel matrices and all its once-shifted Hankel matrices are positive semidefinite.

Multiplication by \(D\), which is strictly positive on \([0,\infty)\), gives
\[
 \int t^kD(t)\,d\nu(t)=\mathcal L_X(Dt^k)=\mu_s(Wt^k).
\]
The right side is independent of \(X\). Let \(\eta\) be the positive measure whose moments are these numbers, obtained by pushing forward
\(W(y^2)w_s(y)\,dy\) under \(t=y^2\).

These moments determine \(\eta\) uniquely on \([0,\infty)\). Here is a direct determinacy argument, avoiding an unexplained application of Carleman's criterion. The density satisfies
\[
 w_s(y)=O(y^se^{-2\pi y})\quad(y\to\infty)
\]
and is bounded near zero. After the fixed polynomial factor \(W(y^2)\), the symmetric square-root lift of \(\eta\) to the real line still has a finite exponential moment \(\int e^{c|y|}\) for every \(0<c<2\pi\). Any other positive measure with the same moments has the same finite integral of \(\cosh(c\sqrt t)\), by its nonnegative power series and monotone convergence. Its symmetric square-root lift therefore also has a finite exponential moment. The two lifted measures have identical moments; their analytic Laplace transforms coincide near zero and hence throughout their common strip. Uniqueness of the Fourier transform gives equality of the lifted measures and then of the original measures.

Thus \(D\nu=\eta\), and division by the positive function \(D\) forces \(\nu\) to equal the original positive measure representing \(\mathcal L_\xi\). Therefore \(\mathcal L_X=\mathcal L_\xi\). Their difference on a polynomial is
\[
 (X-\xi)\sum_{N<j\le K}c_jq(-j^2).
\]
A Lagrange polynomial supported at any one of these distinct nodes makes the sum nonzero. Hence \(X=\xi\).

The nested compact intervals have intersection \(\{\xi\}\), so their endpoints converge to \(\xi\). This completes the proof.

The existence theorem is a standard result of the Stieltjes moment problem; see N. I. Akhiezer, *The Classical Moment Problem* (English translation, 1965), or K. Schmüdgen, *The Moment Problem* (2017). The exact theorem invoked has been stated above; no claim to have newly proved the general representation theorem is made.

## 5. A one-pole specialization gives rational, rather than algebraic, endpoints

Take \(N=0,K=1\), so \(D=t+1\) and \(h=1\). The coefficient of \(X\) in each pencil has rank one. Consequently each determinant is linear in \(X\), and both interval endpoints are rational.

There is an especially transparent Schur-complement formula. Write
\[
 M_e=\mu_s(t^e)
 =\frac{(-1)^eB_{2e+2}(2e+s)!}{(s-1)!(2e+2)!},
 \quad a_s=\frac12+\frac1{s-1},\quad M_0=\frac{s}{12}.
\]
Use the polynomial basis
\(1,(t+1),t(t+1),\ldots,t^{m-2}(t+1)\).
For \(0\le i,j<m-1\), put
\[
 A_{ij}=M_{i+j+1}+M_{i+j},\quad b_i=M_i,
\]
\[
 B_{ij}=M_{i+j+2}+M_{i+j+1},\quad c_i=M_{i+1}.
\]
Both \(A\) and \(B\) are positive definite Gram matrices, independent of \(X\). The two transformed pencils are
\[
 \begin{pmatrix}X-a_s&b^T\\b&A\end{pmatrix},\qquad
 \begin{pmatrix}a_s+M_0-X&c^T\\c&B\end{pmatrix}.
\]
Thus
\[
 \boxed{L_m=a_s+b^TA^{-1}b<\zeta(s)
 <a_s+M_0-c^TB^{-1}c=U_m.}
\]
For \(m=1\), the empty quadratic forms are zero. These rational endpoints are nested and converge by §4.

This also identifies the errors as **best weighted polynomial approximation errors**. If \(d\rho\) is the positive measure representing \(\mu_s\), then
\[
 \zeta(s)-L_m=
 \min_{\deg p<m-1}\int_0^\infty
 (t+1)\left(p(t)-\frac1{t+1}\right)^2d\rho(t),
\]
\[
 U_m-\zeta(s)=
 \min_{\deg p<m-1}\int_0^\infty
 t(t+1)\left(p(t)-\frac1{t+1}\right)^2d\rho(t).
\]
For \(m=1\), read the minimization as the zero polynomial alone. Expanding the squares and minimizing the resulting positive quadratic forms proves the identities. The two Gram matrices approximate the same reciprocal function in complementary weights. This explains both the enclosure and why adding polynomial directions improves it.

For example:

| \(s\) | \(m\) | \(L_m\) | \(U_m\) |
|---|---:|---:|---:|
| 3 | 1 | \(1\) | \(5/4\) |
| 3 | 2 | \(19/16\) | \(29/24\) |
| 3 | 3 | \(259/216\) | \(65/54\) |
| 5 | 1 | \(3/4\) | \(7/6\) |
| 5 | 2 | \(203/204\) | \(161/152\) |
| 5 | 3 | \(3315/3232\) | \(8951/8576\) |

**Subsequently proved Wilson-polynomial identity.** For \(s=3\), the first five computed pairs match
\[
 L_m=H_m^{(3)}+\frac{m-1}{2m^3},\qquad
 U_m=H_m^{(3)}+\frac1{2m(m+1)}.
\]
The identification is now proved for every matrix size in `one-pole-schur-structure.md`: the residual minimizers are Wilson polynomials with parameters (1,1,1,0) and (1,1,1,1). Their factorial norms give telescoping endpoint increments. The general convergence theorem above does not depend on that special evaluation.

## 6. Scope and reproducibility

The matrix construction gives rational or algebraic enclosures and an intuitive least-squares interpretation. It does not make the denominators of the rational approximants small relative to their errors, and does not prove irrationality. The fixed-denominator convergence theorem has no effective rate here. The increasing-\(K\) nesting theorem is separate from it.

Run `python codex/computations/check_localizing_zeta_enclosures.py` for exact determinant, root-isolation, and Schur-complement checks. The general convergence and inertia assertions are proved above, not inferred from those finite computations.
