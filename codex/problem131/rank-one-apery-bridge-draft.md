# The remaining Apéry bridge: a proved moving-pole identity and one explicit gap

**Later resolution (15:20 UTC):** the original Schur–Apéry identity is now
proved and accepted in two independent audits; see
[the Heun/RII proof](rank-one-original-apery-heun-rii-proof.md).
Unresolved labels below describe this earlier investigation only.


**Status: partial proof, with the unresolved step isolated.** The terminating Wilson contiguous identity and the general associated-polynomial identity below are proved. They identify a classical moving-pole Schur upper bound with the Prévost finite-difference quotient. The equality of that moving-pole bound with the newly constructed D_(2h) rank-one determinant root is still unproved, despite exact agreement through h=12. No sample-based theorem is claimed.

This is a checkpoint of the bounded proof investigation requested by the parent. Primary context and notation are in `rank-one-primary-source-context.md`; the independent construction is in `rank-one-zeta-subspace-approximants.md`.

## 1. A terminating identity identifies the relevant positive Wilson polynomial

Write Prévost's polynomial as
\[
 \Pi_h(y)=\sum_{k=0}^h
 \frac{(-h)_k(h+1)_k(-y)_k(y+1)_k}{(k!)^4}.
\]
It satisfies Pi_h(y)=Pi_h(-1-y). Put
\[
 d_h(y)=\frac{\Pi_h(y)-\Pi_h(y-1)}{2y}.
\]
The quotient is an even polynomial. For
\(\phi_k(y)=(-y)_k(y+1)_k/(k!)^2\), direct cancellation gives
\[
 \phi_k(y-1)/\phi_k(y)=(y-k)/(y+k).
\]
Consequently, after y=ix and t=x^2,
\[
 \frac{\phi_k(y)-\phi_k(y-1)}{2y}
 =-\frac{k}{(k!)^2}D_{k-1}(t).
\]
Reindexing the finite sum proves
\[
 \boxed{d_h(ix)=h(h+1)
 {}_4F_3\!\left(\begin{matrix}1-h,h+2,1+ix,1-ix\\2,2,2\end{matrix};1\right)
 =\frac{h(h+1)}{(h!)^3}W_{h-1}(t;1,1,1,1).}       \tag{1}
\]
These are ordinary positive-parameter Wilson polynomials for the measure dnu=t dmu_3. No negative-parameter continuation is being used.

This also explains why the naive odd-part polynomial was not the D_(2h) minimizer: it belongs to a different, fixed positive measure. For example its monic form at h=2 is t-1, whereas the D_4 minimizer is t-1/5.

## 2. A moving-pole Schur bound with an explicit associated numerator

Let P_n be any nonzero scalar multiple of the degree-n orthogonal polynomial for dnu=t dmu_3, and a a positive integer. Define
\[
 R_n(a)=\int\frac{1-P_n(t)/P_n(-a^2)}{t+a^2}\,d\nu(t).
\]
The integrand is a polynomial, so this is rational for the Wilson polynomial. Set
\[
 U_n(a)=H_a^{(3)}+\frac1{2a^2}-\frac1{2a^3}+\frac1{4a^4}
                         -\frac{R_n(a)}{a^4}.       \tag{2}
\]
Then
\[
 \boxed{U_n(a)-\zeta(3)=
 \frac1{a^4P_n(-a^2)^2}
 \int\frac{P_n(t)^2}{t+a^2}\,d\nu(t)>0.}           \tag{3}
\]
To prove the constants, use
\[
 \int\frac{d\nu(t)}{t+a^2}
 =\frac14-a^4\zeta(3)+a^4H_a^{(3)}+\frac{a^2}2-\frac a2.
\]
Subtract R_n(a). Orthogonality and the fact that
\((P_n(t)-P_n(-a^2))/(t+a^2)\) has degree n-1 give
\[
 \int\frac{P_n(t)^2}{t+a^2}\,d\nu
 =P_n(-a^2)\int\frac{P_n(t)}{t+a^2}\,d\nu.
\]
This proves (3). It is the optimized one-pole Schur upper bound, now with a varying pole a, rather than the fixed pole a=1 studied earlier.

## 3. A general polynomial identity connects this bound to Prévost's numerator

Define the polynomial functional
\[
 \mathcal L(1)=0,\qquad\mathcal L(x^j)=jB_{j-1}\quad(j\ge1),
 \qquad B_1=-1/2.
\]
Its generating series is z^2/(e^z-1), so
\[
 \mathcal L[f(x+1)-f(x)]=f''(0),\qquad
 \mathcal L[f(-1-x)]=-\mathcal L[f(x)].              \tag{4}
\]
Also the polynomial moment formula gives
\[
 \mu_3(Q(t))=\tfrac12\mathcal L[x^3Q(-x^2)].         \tag{5}
\]

Let Pi be *any* polynomial invariant under x->-1-x, and put
\[
 d(x)=\frac{\Pi(x)-\Pi(x-1)}{2x},\quad
 P(t)=d(i\sqrt t),\quad
 \Theta(y)=\mathcal L_x\left[\frac{\Pi(x)-\Pi(y)}{x-y}\right].
\]
At a positive integer a with d(a)!=0, set
\(C(a)=H_a^{(3)}\Pi(a)+\Theta(a)/2\), and define U(P,a) by (2), even if P is not orthogonal. Then the exact identity is
\[
 \boxed{\frac{C(a)-C(a-1)}{\Pi(a)-\Pi(a-1)}-U(P,a)
 =\frac{\mathcal L[x\Pi(x)]}
 {a^3[\Pi(a)-\Pi(a-1)]}.}                          \tag{6}
\]
Thus the one scalar orthogonality condition L[x Pi]=0 is enough to make the two expressions coincide.

Here are algebraic details to expose the shifts and signs. Put d=d(a) and
\[
 H(x)=\frac{d(x)-d}{x^2-a^2}.
\]
This is an even polynomial. Formula (5) gives
\[
 R(P,a)=-\frac1{2d}\mathcal L[x^5H(x)].             \tag{7}
\]
For g(x)=(Pi(x-1)-Pi(a-1))/(x-a), (4) gives
\[
 \Theta(a)-\Theta(a-1)
 =2\mathcal L[x^2H+a xH]-g''(0).
\]
Even powers above two have L-value zero, while L[x^2]=-1. Reflection symmetry gives Pi(-1)=Pi(0), Pi'(-1)=-Pi'(0), Pi''(-1)=Pi''(0), and d(0)=Pi'(0). Expanding g at zero therefore yields
\[
 \Pi(a-1)+\frac{a^3}{2}[\Theta(a)-\Theta(a-1)]
 =\Pi(0)-a d+\frac{a^2}{2}\Pi''(0)+a^4\mathcal L[xH].
\]
Finally (4), applied also to (x+1)^2 Pi(x), gives
\[
 \mathcal L[xd(x)]=\Pi''(0)/2,\qquad
 \mathcal L[x^3d(x)]=\Pi(0)-\mathcal L[x\Pi(x)].
\]
Substitution of these identities and (7) proves (6). The calculation does not assume a conjectural orthogonality relation for the negative Wilson parameter.

For Pi=Pi_h, the condition L[x Pi_h]=0 is the first nontrivial Pade orthogonality condition in Prévost's construction. Here is an independent all-h verification. The positive Wilson integral, with parameters (k+1,1,1,0) and its integrable zero-parameter limit, gives
\[
 \mu_3(D_k)=\frac{(k+1)(k!)^2}{2(k+2)}.
\]
Substituting this in the terminating expansion (1) yields
\[
 \mu_3(d_h(i\sqrt t))
 =\frac{h(h+1)}4
 {}_3F_2(1-h,h+2,1;2,3;1)=\frac12.
\]
The last evaluation is the terminating Saalschutz identity, whose factorial quotient is 2/[h(h+1)]. Formula (5) and the already proved relation L[x^3 d_h]=Pi_h(0)-L[x Pi_h], with Pi_h(0)=1, now give L[x Pi_h]=0. Thus this step is not left as a finite check. Combining this condition with (1) and (6) proves
\[
 \boxed{U_{h-1}(a)=
 \frac{C_h(a)-C_h(a-1)}{\Pi_h(a)-\Pi_h(a-1)}.}        \tag{8}
\]
The separately developed all-h mixed numerator identity then gives the adjacent-Apéry formula for this *moving-pole* family at a=h. That separate identity does not by itself identify the D_(2h) family.

## 4. The remaining statement is one precise optimized-norm identity

The desired Schur identification is equivalent to
\[
 \boxed{
 \min_{\substack{p\ \mathrm{monic}\\\deg p=h-1}}
 \int\frac{p(t)^2}{D_{2h}(t)}\,d\mu_3(t)
 \stackrel{?}=
 \frac1{h^4P_{h-1}(-h^2)^2}
 \int\frac{P_{h-1}(t)^2}{t+h^2}\,t\,d\mu_3(t),}
\]
where P_(h-1) is the ordinary Wilson polynomial for t mu_3. Both sides are positive and optimized Schur errors for their respective spaces. Exact agreement was checked for h=1,...,12. Their minimizers are different, so identifying the two polynomials directly is false.

For reference, the first D_(2h) minimizers are
\[
 1,\quad t-1/5,\quad t^2-(191/149)t+20/149,
\]
\[
 t^3-(5038/1363)t^2+(3427/1363)t-252/1363.
\]
A possible proof could use rational Gaussian quadrature duality, the rational-modification determinant identity, or a closed Christoffel/Toda relation. The elementary relation
\[
 \mu_3(t^e/D_{K-1})=
 \mu_3(t^{e+1}/D_K)+K^2\mu_3(t^e/D_K)
\]
is exact, but by itself it does not close a recurrence along K=2h. No such closure or norm-duality theorem has been proved in this checkpoint.

## 5. Reproducibility and scope

The exact computations in `rank_one_subspace_approximants.json` establish the finite D_(2h) roots. The separate checker accompanying this note verifies the contiguous identity, general functional identity (6), moving-pole equality through h=12, and the relevant finite Padé orthogonality conditions. These checks support the proved identities; they do not convert the unresolved equality in Section 4 into a theorem.

Independent audit: [rank-one-apery-bridge-independent-review.md](rank-one-apery-bridge-independent-review.md) accepts the proved moving-pole identities and explicitly retains the unresolved D_(2h) norm equality. Its separate moment-Gram checker passes 8 minimizers/adjacent-Apéry bounds and 18 general functional correction cases. The variational characterization on the moving-pole side minimizes over degree-at-most h−1 residuals q with q(−a²)=1, against t dμ₃(t)/(t+a²), and then divides by a⁴; its minimizer is P_(h−1)(t)/P_(h−1)(−a²).
