# Repeated Rodrigues, simple roots, and confluent Jacobi–Piñeiro orthogonality

**Status.** The root, Mellin, confluence, and root-statistic statements below have complete proofs and exact computational checks. They have not yet received a separate agent's review. The Jacobi–Piñeiro identification is a comparison with a primary paper, not a historical novelty claim. Root location does **not** prove positivity of the higher-zeta forms; Section 6 gives an exactly certified counterexample to that inference.

## 1. A root theorem needing only Rolle's theorem

Put \(N=mn\), where \(m,n\ge1\), and
\[
 R_nf=\frac1{n!}D^n(x^nf),\qquad
 J_{n,m}=R_n^m(1-x)^N
 =\sum_{k=0}^N(-1)^k\binom Nk\binom{n+k}k^m x^k.
\]
**Theorem.** Every one of the \(N\) roots of \(J_{n,m}\) is simple and belongs to \((0,1)\).

Here is a stronger induction that makes the mechanism visible. After \(r\) applications of \(R_n\), the polynomial has degree \(N\), exactly \(rn\) simple roots in \((0,1)\), a root of multiplicity \(N-rn\) at 1, and no other roots. Its value at 0 is 1. This holds for \(r=0\).

Suppose it holds for \(r<m\), and multiply that polynomial by \(x^n\). Before differentiating, there are \(rn\) simple interior roots, an endpoint root of multiplicity \(n\) at 0, and one of multiplicity \(N-rn\ge n\) at 1. After \(j\) derivatives, \(0\le j\le n\), Rolle's theorem gives at least \(rn+j\) distinct interior roots; the endpoint multiplicities are at least \(n-j\) and \(N-rn-j\). The total counted multiplicity is
\[
 rn+j+(n-j)+(N-rn-j)=N+n-j,
\]
exactly the degree of that derivative. Thus every asserted multiplicity is exact, every interior root is simple, and there are no uncounted roots. At \(j=n\) this proves the next induction step. Division by \(n!\) has no effect on roots. The constant term stays 1 because \(R_n1=1\).

The picture is elementary: each differentiation moves one zero from the available endpoint multiplicities into every intervening Rolle interval; degree counting prevents any additional multiplicity. Repeating the block \(m\) times consumes the entire reservoir of zeros at 1.

## 2. The same reservoir creates repeated Mellin zeros

If \(f\) has at least \(n\) zeros at 1, integration by parts \(n\) times gives
\[
 \int_0^1x^{t-1}R_nf(x)\,dx
 =\frac{(-1)^n}{n!}\prod_{j=1}^n(t-j)
   \int_0^1x^{t-1}f(x)\,dx,\qquad\Re t>0.
\]
The boundary terms at 1 vanish by the zero multiplicity. At 0 they are \(O(x^{\Re t})\), so vanish as well. Every intermediate Rodrigues polynomial retains the needed \(n\) zeros at 1 until the final block.

Starting from
\(\int_0^1x^{t-1}(1-x)^Ndx=N!/\prod_{k=0}^N(t+k)\), we obtain
\[
 \boxed{\int_0^1x^{t-1}J_{n,m}(x)\,dx
 =\frac{(-1)^NN!}{(n!)^m}
 \frac{\prod_{j=1}^n(t-j)^m}{\prod_{k=0}^N(t+k)}.}
\]
Consequently
\[
 \int_0^1 x^k(\log x)^rJ_{n,m}(x)\,dx=0,
 \quad0\le k<n,\quad0\le r<m.
\]
The same construction therefore explains two apparently different features: consuming endpoint zeros gives interior polynomial zeros, while integration by parts converts that reservoir into multiple Mellin zeros. No residue computation is needed for this derivation.

## 3. A precise primary-source identification

Walter Van Assche and Els Coussement, *Some classical multiple orthogonal polynomials*, Journal of Computational and Applied Mathematics **127** (2001), 317–347, DOI [10.1016/S0377-0427(00)00503-3](https://doi.org/10.1016/S0377-0427(00)00503-3), [author preprint math/0103131](https://arxiv.org/abs/math/0103131), Section 3.1, equations (3.20)–(3.24), defines Jacobi–Piñeiro polynomials for weights
\(x^{\alpha_j}(1-x)^{\alpha_0}\), with \(\alpha_j>-1\) and pairwise nonintegral differences. It gives their commuting Rodrigues operators explicitly. The downloaded primary text is `codex/sources/van-assche-coussement-classical-multiple-0103131.pdf` and its `.txt` extraction.

Set \(\alpha_0=0\), use the balanced multiindex \((n,\ldots,n)\), and normalize the constant term to 1 instead of making the polynomial monic. Define
\[
 S_{n,\alpha}f=\frac{x^{-\alpha}}{(\alpha+1)_n}
 D^n[x^{n+\alpha}f].
\]
It maps \(x^k\) to
\((k+\alpha+1)_n/(\alpha+1)_n\,x^k\). Thus these operators commute: they are polynomials in \(xD\), divided by constant normalizers. Applying them to \((1-x)^N\) gives
\[
 J_{\boldsymbol\alpha}(x)=\sum_{k=0}^N(-1)^k\binom Nk
 \prod_{j=1}^m\frac{(k+\alpha_j+1)_n}{(\alpha_j+1)_n}x^k.
\]
Its Mellin transform is
\[
 \frac{(-1)^NN!}{\prod_j(\alpha_j+1)_n}
 \frac{\prod_{j=1}^m\prod_{r=1}^n(t-\alpha_j-r)}
      {\prod_{k=0}^N(t+k)}.
\tag{*}
\]
One may prove (*) by the same integration by parts, initially for sufficiently large \(\Re t\), and then rational continuation. Each operator consumes at most \(n\) endpoint zeros. The result is orthogonal to \(x^{\alpha_j+k}\), \(0\le k<n\), exactly the defining conditions in the primary source. When the parameter differences are nonintegral, the \(N\) Mellin zeros are distinct. A polynomial of degree at most \(N\) has a Mellin numerator of degree at most \(N\), so these zeros and the constant normalization give uniqueness directly.

Let the distinct \(\alpha_j\) approach zero. The displayed finite coefficient formula converges to \(J_{n,m}\); the separate Mellin zeros coalesce into order-\(m\) zeros. This proves that our polynomial is a **confluent limit of Jacobi–Piñeiro polynomials**. It is not legitimate simply to substitute identical parameters into the source's orthogonality definition and claim that its generic uniqueness theorem applies: the ordinary weights would then coincide. The confluent conditions retain the derivatives of those weights, namely \(x^k(\log x)^r\), and do give uniqueness.

Equivalently, take divided differences of the vanishing moments at the distinct \(\alpha_j\), then let the parameters merge. Their limits are logarithmic moments divided by \(r!\). Interchanging the limit and integral is justified by an integrable bound \(C x^{-\delta}(1+|\log x|^{m-1})\), with \(0<\delta<1\).

## 4. Exact root statistics without solving for any roots

Write the roots as \(\xi_1,\ldots,\xi_N\). The leading coefficient is
\((-1)^N\binom{(m+1)n}{n}^m\). The endpoint value is
\[
 J_{n,m}(1)=(-1)^N\frac{(mn)!}{(n!)^m}.
\]
For example this follows by multiplying its Mellin transform by \(t\) and letting real \(t\to\infty\): for any continuous polynomial \(P\), \(t\int_0^1x^{t-1}P(x)dx\to P(1)\).

Vieta's formulas now give
\[
 \prod_i\xi_i=\binom{(m+1)n}{n}^{-m},\qquad
 \prod_i(1-\xi_i)=\frac{(mn)!}{(n!)^m}\binom{(m+1)n}{n}^{-m},
\]
\[
 \boxed{\sum_i\xi_i=mn\left(\frac m{m+1}\right)^m.}
\]
The last equality follows by taking the ratio of the last two coefficients. Thus the average root is \((m/(m+1))^m\), tending to \(e^{-1}\) as \(m\to\infty\). These identities are consequences of the explicit classical family, not a novelty assertion.

## 5. An unbalanced extension

The same proof works with positive integers \(n_1,\ldots,n_m\) and \(N=\sum n_j\):
\[
 J_{\boldsymbol n}(x)=\prod_{j=1}^mR_{n_j}(1-x)^N
 =\sum_{k=0}^N(-1)^k\binom Nk\prod_j\binom{n_j+k}k x^k.
\]
All its \(N\) roots are simple in \((0,1)\), and
\[
 \int_0^1x^{t-1}J_{\boldsymbol n}(x)dx
 =\frac{(-1)^NN!}{\prod_jn_j!}
 \frac{\prod_j\prod_{r=1}^{n_j}(t-r)}{\prod_{k=0}^N(t+k)}.
\]
The zero at \(t=r\) has multiplicity \(\#\{j:n_j\ge r\}\). This makes it possible to design a nonuniform staircase of logarithmic orthogonality conditions. Whether such a design improves a useful irrationality form remains open here; the root theorem by itself supplies no denominator or decay improvement.

## 6. Why this does not repair positivity of the zeta(5) form

The pure form from `multiplicity-lifted-zeta-family.md` is
\[
 \mathscr L_{n,2}=\frac1{3!}\iint_{(0,1)^2}
 J_{n,2}(x)J_{n,2}(y)\frac{[-\log(xy)]^3}{1-xy}\,dx\,dy.
\]
Although the kernel is pointwise positive, this is not a pointwise positive integrand: the polynomial factors change sign. Nor is this kernel positive semidefinite on the entire polynomial space in the required sense. The exact same family supplies opposite signs:
\[
 \mathscr L_{1,2}>0,\qquad\mathscr L_{3,2}<0.
\]
Their approximate values are printed by the verification script, but the sign certificates use rational inequalities only. Write
\[
 \zeta(5)=\frac{16}{15}\sum_{k\ge0}2^{-k-1}
 \sum_{j=0}^k(-1)^j\binom kj\frac1{(j+1)^5}.
\]
Each inner difference equals
\(\int_0^1(1-x)^k(-\log x)^4dx/24\), hence lies in \([0,1]\). Truncating at \(k=M-1\) gives a rational lower bound and an upper bound only \(16/(15\cdot2^M)\) larger. With \(M=80\), substitution into the exact rational forms from formula (7) of the original note proves the two signs, including the small negative value at \(n=3\).

Thus simple interior zeros, simultaneous logarithmic orthogonality, and a positive scalar kernel do not imply positivity of this quadratic integral. The source of classical Beukers positivity is an additional integral transformation; it is not a consequence of root location alone. The positive mixed-zeta period in `positive-period-amplitude-mixing.md` makes the complementary tradeoff explicit: positivity survives there, but additional zeta values appear.

## 7. Reproduction and useful next directions

Run `python codex/computations/check_confluent_roots.py`. It checks 16 polynomial cases by exact root counts and gcds, the endpoint and root-sum formulas, shifted Mellin identities, and the rational sign enclosures. The theorem proofs above do not depend on those finite checks.

The useful constructive conclusion is an exact design language: Rodrigues blocks correspond to selected Mellin zeros; collisions of ordinary Jacobi weights produce logarithmic orthogonality; endpoint multiplicity guarantees real simple roots. Further progress would need an additional ingredient linking this design to a positive transformed integral or sharper primitive denominator estimates. Those ingredients have not been proved by the root theorem.
