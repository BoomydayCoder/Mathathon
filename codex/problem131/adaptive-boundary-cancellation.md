# Cancelling the total boundary contribution with one adaptive parameter

**Status.** The exact cancellation formula, eventual existence of the tuning parameter, and the sharp analytic asymptotic below are proved. The asymptotic uses the previously proved and independently reviewed derivative-tail asymptotic in `truncated-derivative.md`. Exact rational computations at eleven indices through 40 measure the arithmetic cost; logarithmic error sizes in the table are high-precision numerical evaluations, not interval certificates. This is an analytically successful construction with a serious observed arithmetic cost, not an irrationality proof or a general impossibility theorem. Historical novelty is unresolved.

## 1. Cancel a sum instead of imposing repeated zeros separately

Write
\[
 P_n(t)=\frac{\prod_{j=1}^n(t-j)}{\prod_{k=0}^n(t+k)},
 \quad R_n=P_n^2,\quad b=n+1,
 \quad w_{n,a}(t)=\frac{t-a}{t+b}.
\]
The proposed form is
\[
 \mathcal L_{n,a}=-\frac16\sum_{t\ge1}(R_nw_{n,a})'''(t).
\]
Unlike the higher-multiplicity construction, this leaves the double zeros of \(R_n\) in place. We choose \(a\) to make the **sum** of their third derivatives vanish.

For \(1\le j\le n\), define
\[
 u_j=P_n'(j)=(-1)^{n-j}\frac{(j-1)!^2(n-j)!}{(n+j)!},
 \qquad\delta_j=H_{n-j}+H_{n+j}-2H_{j-1}.
\]
Then \(R_n''(j)=2u_j^2\) and \(R_n'''(j)=-12u_j^2\delta_j\). Set
\[
 D=\sum_j u_j^2\delta_j,\quad
 E=\sum_j\frac{u_j^2\delta_j}{j+b},\quad
 G=\sum_j\frac{u_j^2}{(j+b)^2},\quad W=2E+G.
\]
Because \(w=1-(a+b)/(t+b)\), direct differentiation gives
\[
 \sum_{j=1}^n(R_nw)'''(j)=-12D+6(a+b)W.
\]
Thus, whenever \(W\ne0\), the unique tuning is
\[
 \boxed{a_n=\frac{2D}{W}-(n+1).}
\tag{1}
\]
The resulting full form is exactly its tail:
\[
 \boxed{\mathcal L_n=-\frac16\sum_{t>n}(R_nw_{n,a_n})'''(t).}
\tag{2}
\]
There is no subtraction of a rational boundary correction after forming the integral.

## 2. The partner is only a combination of two adjacent Legendre polynomials

Let \(L_n(x)=D^n[x^n(1-x)^n]/n!\), so its Mellin transform is \((-1)^nP_n(t)\). Since
\[
 P_{n+1}=P_n\frac{t-b}{t+b},
\]
putting \(d_n=(a_n+b)/(2b)=D/(bW)\) gives
\[
 P_nw=(1-d_n)P_n+d_nP_{n+1}.
\]
The polynomial partner is therefore
\[
 \boxed{Q_n=(1-d_n)L_n-d_nL_{n+1}.}
\tag{3}
\]
The minus sign in the polynomial identity is essential: adjacent Legendre Mellin transforms have opposite prefactors. Its degree is \(n+1\) whenever \(d_n\ne0\).

The ordinary monomial-kernel identity now gives
\[
 \mathcal L_n=\frac16\iint\frac{L_n(x)Q_n(y)[-\log(xy)]^3}{1-xy}\,dx\,dy
 =A_n^*\zeta(5)-B_n^*,
\]
with rational coefficients. Thus purity is automatic, even though the pair is asymmetric. Once \(Q_n\) is multiplied by an integer to have integer coefficients, the usual sufficient denominator bound is \(D_{n+1}^5\) for its rational part. The cost of making \(Q_n\) integral must be counted separately.

The first tunings are
\[
 a_1=\frac7{10},\qquad a_2=\frac{2942}{3561},\qquad
 a_3=\frac{563739}{651481}.
\]
For example, \(d_1=27/40\) and \(20Q_1=-7+68x-81x^2\), giving
\[
 20\mathcal L_1=-572\zeta(5)+\frac{18979}{32}<0.
\]
This isolated negative sign is compatible with the eventual positivity proved below.

## 3. The tuning exists and stays bounded for every sufficiently large index

Let \(\delta_1=H_{n-1}+H_{n+1}\). The estimates in `log-kernel-extension.md` imply
\[
 \sum_{j=2}^n u_j^2|\delta_j|=O(n^{-8}\log n),\qquad
 \sum_{j=2}^n u_j^2=O(n^{-8}).
\]
For completeness, write \(|u_j|=1/[j^2\binom nj\binom{n+j}j]\). For \(2\le j\le n\), the product of binomials is bounded below by a fixed positive multiple of \(n^4\); its ratio crosses 1 at most once, so its minimum is at an endpoint. The \(j=2\) endpoint is polynomial of degree four, and the central-binomial endpoint is larger for \(n\ge5\). Summing \(j^{-4}\) proves the estimates, since \(|\delta_j|=O(\log n)\).

As \(u_1^2=1/[n^2(n+1)^2]\), we obtain
\[
 D=u_1^2\delta_1+O(n^{-8}\log n),
\]
\[
 W=u_1^2\left(\frac{2\delta_1}{n+2}+\frac1{(n+2)^2}\right)
       +O(n^{-9}\log n).
\]
Both \(D\) and \(W\) are positive eventually. Their ratio yields
\[
 \boxed{a_n=1-\frac{n+2}{2\delta_1(n+2)+1}+O(n^{-3})
 =1-\frac1{4\log n}+O((\log n)^{-2}).}
\tag{4}
\]
In particular \(0<a_n<1\) for all sufficiently large \(n\), and \(d_n\to1/2\). This proves eventual nondegeneracy without relying on the finite samples. It does not assert existence or the inequalities at every small index; those were checked exactly at the listed sample indices.

## 4. A contour comparison proves the actual tail asymptotic

This argument handles any bounded real sequence \(a_n\), so it separates the analytic estimate from the complicated rational tuning.

Define
\[
 U_n(z)=\frac{\Gamma(z)^2\Gamma(n+1-z)}{\Gamma(n+1+z)}.
\]
Reflection gives \(R_n(z)=\sin^2(\pi z)U_n(z)^2/\pi^2\). Take the upward line \(c=n+1/2\). A clockwise cotangent contour, followed by three integrations by parts, gives
\[
 T_n(w):=-\frac16\sum_{t>n}(R_nw)'''(t)
 =-\frac\pi{12}\int_{\mathbb R}
 [4-6\operatorname{sech}^2(\pi y)]U_n(c+iy)^2w(c+iy)\,dy.
\tag{5}
\]
Indeed \(K(z)=\pi\cot\pi z\) satisfies
\(K'''=4\pi^4\csc^2\pi z-6\pi^4\csc^4\pi z\). This checks the factor and sign. The line meets no integer, the added pole \(-b\) is outside the right closure, and all integration-by-parts boundary terms vanish by rational decay before applying reflection. The original contour can be closed on half-integer right edges; thus the calculation does not require an unproved representation of the tuned linear form.

Factor
\[
 U_n(z)=B(z,n+1-z)B(z,n+1).
\]
The beta integral gives \(|B(c+iy,n+1)|\le B(c,n+1)\). Also
\(|\Gamma(c+iy)|\le\Gamma(c)\) by its integral and
\(|\Gamma(1/2-iy)|^2=\pi/\cosh\pi y\). Consequently
\[
 \boxed{|U_n(c+iy)|^2\le U_n(c)^2\operatorname{sech}(\pi y),
 \qquad U_n(c)^2\sim\frac{\pi^2}{2n^2}16^{-n}.}
\tag{6}
\]
The last relation is Stirling's formula. Since \(|4-6\operatorname{sech}^2|\le4\), (5) immediately gives \(T_n(w)=O(n^{-2}16^{-n})\) for bounded \(a_n\).

More precisely,
\[
 w(c+iy)-\frac12=\frac{-1/2-2a_n+iy}{2(c+b+iy)},
\]
so its absolute value is \(O((1+|y|)/n)\). The majorant in (6) has finite first moment, hence
\[
 \boxed{T_n(w)=\frac12T_n(1)+O(n^{-3}16^{-n}).}
\tag{7}
\]
The previously established asymptotic in `truncated-derivative.md` is
\[
 E_n:=-\frac1{24}\sum_{t>n}R_n'''(t)
 \sim\frac\pi2\frac{\log n}{n^3}16^{-n}.
\]
Here \(T_n(1)=4E_n\), so (2), (4), and (7) prove
\[
 \boxed{\mathcal L_n=2E_n+O(n^{-3}16^{-n})
 \sim\pi\frac{\log n}{n^3}16^{-n}.}
\tag{8}
\]
In particular the tuned forms are strictly positive for all sufficiently large \(n\) and have exact exponential rate \(1/16\). The construction is therefore not a vanishing identity or a merely formal cancellation.

## 5. Exact arithmetic cost: small examples already contain large new primes

Let \(s_n>0\) be the rational scale making \(s_nQ_n\) an integer polynomial with coprime coefficients. This scale is computed exactly from all its coefficients, not estimated from \(d_n\)'s denominator alone. The first values are
\[
 s_1=20,\quad s_2=10683=9\cdot1187,\quad
 s_3=2605924=4\cdot651481.
\]
Thus the adaptive parameter has introduced primes much larger than \(n+1\), even though the partner degree has increased by only one. These primes are **not** covered by the usual \(D_{n+1}^5\) coefficient bound. This is a finite, exact obstruction to ignoring the cost of tuning; it is not an infinite-prime or asymptotic lower-bound theorem.

For a completely general sufficient bound, the earlier factorial-valuation estimate gives
\(D_nD_{2n}u_j\in\mathbb Z\). Therefore
\[
 K_n=D_n^2D_{2n+1}^4
\]
clears both \(D\) and \(W\). If \(V_n=K_nW\in\mathbb Z\), then \(d_n=(K_nD)/[(n+1)V_n]\), so \((n+1)|V_n|\) clears the partner coefficients. This sufficient bound is large and no claim is made that it captures their primitive size sharply.

The table uses exact rational coefficients and high-precision evaluations of the errors. In the final column, the *linear form itself* has been made primitive: if \(s_n\mathcal L_n=A\zeta(5)-B\), with \(A\in\mathbb Z\) and reduced \(B=p/q\), multiply by \(q/\gcd(A,p)\). This accounts for possible cancellation beyond primitive polynomial normalization.

| n | a_n (approx.) | log(s_n)/n | log|L_n|/n | log|primitive integer form|/n |
|---:|---:|---:|---:|---:|
| 2 | .826172 | 4.63820 | −4.47528 | 2.90945 |
| 5 | .891960 | 5.99061 | −3.52163 | 4.90452 |
| 10 | .915167 | 6.06405 | −3.28501 | 4.92293 |
| 20 | .930712 | 7.14210 | −3.11529 | 6.31807 |
| 30 | .937536 | 7.33367 | −3.03645 | 7.05239 |
| 40 | .941646 | 7.29300 | −2.98968 | 6.79326 |

These bounded observations show no approach to a small primitive integer form. They do not prove an asymptotic arithmetic obstruction, and the numerical signs in this table are distinct from the eventual-sign theorem in Section 4.

## 6. Reproduction and interpretation

Run `python codex/computations/check_adaptive_boundary_cancellation.py` for independent exact derivative cancellations, Mellin partner identities, and cotangent normalization at five indices. Run `python codex/computations/explore_adaptive_boundary_cancellation.py`; exact fractions and the eleven sample rows are preserved in `adaptive_boundary_cancellation_samples.json`. The script constructs the partner directly from Legendre coefficients, evaluates the rational monomial integral coefficient by coefficient, and verifies the standard denominator bound after making the partner integral. Its high-precision evaluations are explicitly exploratory.

This route answers the proposed question positively on the analytic side: cancelling the total early-root contribution works, keeps degrees \(n,n+1\), isolates a single zeta value, and gives an eventually positive exponential tail. Its unresolved problem is the arithmetic size of the adaptive coefficient. The fixed harmonic-polynomial no-go results do not imply failure here; any stronger negative conclusion would require new control of the reduced rational quantities in (1), or of their cancellation in the final primitive linear form.
