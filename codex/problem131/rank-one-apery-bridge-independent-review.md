# Independent review: the moving-pole Wilson–Apéry bridge

**Status: accepted for the proved moving-pole family and general functional
identity. The original \(D_{2h}\) Schur-root identification remains open.**
This reviews [the bridge note](rank-one-apery-bridge-draft.md), together with
the already reviewed [mixed-numerator theorem](apery-mixed-numerator-independent-review.md).
Classical Wilson orthogonality is an explicit input. No historical novelty
or new higher-zeta irrationality statement is asserted.

## 1. Measures and the exact variational problem

Use the pushforward under \(t=y^2\) of

\[
 d\mu_3(t)=w_3(y)\,dy,\qquad
 w_3(y)=\pi^2y^3\coth(\pi y)\operatorname{csch}^2(\pi y),
 \qquad d\nu=t\,d\mu_3.
\]

The first notation denotes a pushforward, so no additional Jacobian belongs
in an integral already written in \(y\). The polynomial moments and integer
pole values are

\[
 \mu_3(t^j)=\frac{(-1)^j(2j+3)B_{2j+2}}2,\qquad
 \mu_3((t+a^2)^{-1})
 =a^2[\zeta(3)-H_a^{(3)}]-\frac12+\frac1{2a}.
\]

These normalizations agree with the reviewed
[positive moment construction](odd-moment-generalization.md) and
[one-pole Wilson review](one-pole-independent-review.md). In particular,
\(\mu_3(1)=1/4\) and \(\nu(1)=1/12\). Polynomial division gives

\[
 \int\frac{d\nu(t)}{t+a^2}
 =\frac14-a^4\zeta(3)+a^4H_a^{(3)}+\frac{a^2}2-\frac a2.
\]

All moments are finite; the density is positive on the positive half-line
and decays exponentially in \(y\).

Let \(P_n\) be a nonzero real scalar multiple of the degree-\(n\) orthogonal
polynomial for \(\nu\), and let \(a\) be a positive integer. Its zeros are
positive, so \(P_n(-a^2)\ne0\). The precise variational statement is

\[
 \boxed{
 U_n(a)-\zeta(3)
 =\frac1{a^4}
 \min_{\substack{\deg q\le n\\q(-a^2)=1}}
 \int\frac{q(t)^2}{t+a^2}\,d\nu(t)
 =\frac{\displaystyle\int P_n(t)^2(t+a^2)^{-1}\,d\nu(t)}
 {a^4P_n(-a^2)^2}>0.
 }
\]

Indeed the allowed variations are \((t+a^2)r(t)\), \(\deg r<n\).
For \(q_*=P_n/P_n(-a^2)\), the cross term vanishes by orthogonality for
the **original measure \(\nu\)**. The remaining quadratic term is strictly
positive unless the variation is zero. For \(n=0\), the sole residual is
\(q=1\); there is no omitted minimization variable. This verifies the
normalization, degree constraint, uniqueness, and strict sign.

Equivalently,
\[
 \int\frac{P_n^2}{t+a^2}\,d\nu
 =P_n(-a^2)\int\frac{P_n}{t+a^2}\,d\nu,
\]
because \((P_n-P_n(-a^2))/(t+a^2)\) has degree \(n-1\).
This identity and the displayed pole moment give exactly the note's formula
for \(U_n(a)\). Its rationality follows from the polynomial integrand in
\(R_n(a)\), the rational moments, and the rational Wilson coefficients.

## 2. The all-degree contiguous identity

For \(h\ge1\), set

\[
 \Pi_h(y)=\sum_{k=0}^h
 \frac{(-h)_k(h+1)_k(-y)_k(y+1)_k}{(k!)^4},
 \qquad
 d_h(y)=\frac{\Pi_h(y)-\Pi_h(y-1)}{2y}.
\]

Reflection \(\Pi_h(y)=\Pi_h(-1-y)\) makes the numerator odd and the quotient
an even polynomial, including its removable value at zero. Direct division
of the \(k\)-th term gives
\[
 \frac{\phi_k(y)-\phi_k(y-1)}{2y}
 =\frac{k\phi_k(y)}{y(y+k)}.
\]
At \(y=ix\), cancellation of the factors with indices zero and \(k\)
gives \(-kD_{k-1}(x^2)/(k!)^2\), with exactly this minus sign.
Reindexing proves
\[
 d_h(ix)=h(h+1)
 {}_4F_3\left(\begin{matrix}1-h,h+2,1+ix,1-ix\\2,2,2\end{matrix};1\right)
 =\frac{h(h+1)}{(h!)^3}W_{h-1}(x^2;1,1,1,1).
\]
The conventional Wilson prefactor is
\((2)_{h-1}^3=(h!)^3\), so no scale factor is missing.
The measure \(\nu\) is one half of the usual positive Wilson weight with
these four parameters. Multiplying a measure by a positive constant does
not change its orthogonal polynomials. Thus the variational argument applies
for every \(h\ge1\), with no negative-parameter continuation.

## 3. The functional correction is exact, not an orthogonality guess

Write \(\mathcal L(f)=\mathcal B(f')\), with
\(\mathcal B(x^j)=B_j\) and \(B_1=-1/2\). Its exponential generating series
is \(z^2/(e^z-1)\). It gives
\[
 \mathcal L(f(x+1)-f(x))=f''(0),\qquad
 \mathcal L(f(-1-x))=-\mathcal L(f(x)),\qquad
 \mu_3(Q)=\tfrac12\mathcal L[x^3Q(-x^2)].
\]
The last formula was checked on each monomial, including its sign.

For any reflection-invariant polynomial \(\Pi\), define \(d,P,\Theta,C\)
as in the bridge note and suppose \(d(a)\ne0\). The general identity is

\[
 \frac{C(a)-C(a-1)}{\Pi(a)-\Pi(a-1)}-U(P,a)
 =\frac{\mathcal L(x\Pi(x))}
 {a^3[\Pi(a)-\Pi(a-1)]}.
\]

The following independent algebra checks the sensitive intermediate steps.
For \(d=d(a)\) and \(H(x)=[d(x)-d]/(x^2-a^2)\), one has
\[
 R(P,a)=-\frac1{2d}\mathcal L(x^5H).
\]
Let \(g(x)=[\Pi(x-1)-\Pi(a-1)]/(x-a)\). Taking the shift difference of
\(g\), rather than silently shifting the functional, gives
\[
 \Theta(a)-\Theta(a-1)
 =2\mathcal L(x^2H+axH)-g''(0).
\]
The endpoint derivative is
\[
 g''(0)=-\frac{\Pi''(0)}a
 +\frac{2\Pi'(0)}{a^2}
 -\frac{2[\Pi(0)-\Pi(a-1)]}{a^3}.
\]
Since \(H\) is even and \(\mathcal L(x^2)=-1\),
\(\mathcal L(x^2H)=[\Pi'(0)-d]/a^2\). These formulas give precisely the
intermediate identity stated in the note.

Finally,
\[
 \mathcal L(xd(x))=\frac{\Pi''(0)}2,\qquad
 \mathcal L(x^3d(x))=\Pi(0)-\mathcal L(x\Pi(x)).
\]
For the second equality, apply the shift relation to
\(f(x)=(x+1)^2\Pi(x)\), evaluating \(f''(-1)=2\Pi(0)\); this verifies the
otherwise easy-to-miss correction term. Substitution proves the general
identity without requiring \(\Pi\) to belong to an orthogonal family.
In particular an added constant in \(\Pi\) changes the correction, and
cannot simply be discarded.

For \(\Pi=\Pi_h\), the needed cancellation
\(\mathcal L(x\Pi_h)=0\), \(h\ge1\), already has the elementary finite
Legendre-sum proof in the
[mixed-numerator review](apery-mixed-numerator-independent-review.md).
The bridge note's separate Wilson-integral proof is also correct:
\[
 \mu_3(D_k)=\frac{(k+1)(k!)^2}{2(k+2)},\qquad
 {}_3F_2(1-h,h+2,1;2,3;1)=\frac{2}{h(h+1)}.
\]
The first follows from the half-normalized Wilson integral at
\((k+1,1,1,0)\), with its already justified integrable zero-parameter limit.
For the second, Saalschütz with \(n=h-1\) gives
\[
 \frac{(-h)_{h-1}(1)_{h-1}}{(2)_{h-1}(-h-1)_{h-1}}
 =\frac2{h(h+1)}.
\]
Thus \(\mu_3(d_h(i\sqrt t))=1/2\) as claimed. No finite-sample inference is
needed for either cancellation.

## 4. Constructive corollary: the moving-pole family has full Apéry arithmetic

Combining the established functional identity with the independently
reviewed mixed-numerator balance yields, for every integer \(h\ge1\),

\[
 \boxed{
 U_{h-1}(h)=\frac{5B_h-B_{h-1}}{5A_h-A_{h-1}}>\zeta(3).
 }
\]

Here \(A_0=1,A_1=5,B_0=0,B_1=6\), with the usual Apéry recurrence.
The denominator is positive, since \(A_h\ge5A_{h-1}>0\).
At \(h=1\) the formula gives \(30/24=5/4\), matching the constant-residual
integral bound. There is no \(h=0\) claim.

Let \(\ell_h=\operatorname{lcm}(1,\ldots,h)\), and write
\[
 \widehat A_h=5A_h-A_{h-1},\quad
 \widehat B_h=5B_h-B_{h-1},\quad
 f_h=A_h\zeta(3)-B_h.
\]
The classical denominator property immediately gives
\[
 \widehat A_h\in\mathbb Z,\qquad
 \ell_h^3\widehat B_h\in\mathbb Z.
\]
The positive integral supplies the strict lower bound in
\[
 0<\widehat B_h-\widehat A_h\zeta(3)
 =f_{h-1}-5f_h<f_{h-1}.
\]
The [elementary recurrence proof](discovery-guide.md) gives
\[
 0<f_{h-1}\le\frac5{4h^3A_{h-1}},\qquad
 A_{h-1}\ge C(135/4)^{h-1},\qquad
 \log\ell_h\le ch+O(\log^2(h+1)),
\]
where \(c=\tfrac12\log2+\tfrac34\log3\) and \(e^{3c}<135/4\).
Consequently
\[
 0<\ell_h^3(\widehat B_h-\widehat A_h\zeta(3))\longrightarrow0.
\]
Under rationality \(\zeta(3)=p/q\), multiplying by \(q\) produces a
positive integer tending to zero. Thus this **moving-pole family already
has a complete irrationality argument through its proved Apéry
identification**. This is a reformulation of existing Apéry arithmetic,
not an independent new proof of that arithmetic or a higher-zeta result.

## 5. The remaining family distinction and reproducibility

The original determinant family minimizes a *monic* degree-\(h-1\) norm
for \(\mu_3/D_{2h}\). The proved moving-pole family minimizes a residual
normalized at \(-h^2\) for \(\nu/(t+h^2)\). These are different optimization
problems with different minimizers. Neither Wilson orthogonality nor the
mixed-numerator identity proves equality of their minimum values.
That norm-duality equality, and hence the claimed formula for the original
\(D_{2h}\) determinant root, remains conjectural.

The original checker passed 24 functional identities and twelve finite
Wilson/moving-pole/determinant matches. A second
[independent checker](../computations/check_rank_one_bridge_independent.py)
constructs eight monic minimizers directly from Bernoulli moment Gram
matrices, then checks their equality with the contiguous polynomial and
their bounds against independently generated Apéry recurrence values.
It also checks eighteen functional corrections for mixed polynomial bases
with nonzero constant terms. All checks passed. This second checker does
not read or use the conjectural determinant-root data.

The general proofs, rather than these finite computations, establish the
moving-pole theorem. No unresolved gap was found in that theorem or the
general correction identity; the separate determinant-root equality is
explicitly outside the accepted result.
