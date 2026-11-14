# A table-free upper bound for Fauzan's arithmetic normalization

Status: independently derived elementary estimate, with exact rational computation.
It replaces the integration work in source §5.3 and Appendix B, **conditional
on the source's local arithmetic envelope and Proposition 5.2**. It does not
independently establish the local valuation lemmas, integrality Proposition 5.1,
or all of the claimed irrationality proof. Source: the downloaded preprint
[`../sources/fauzan-zeta5.txt`](../sources/fauzan-zeta5.txt), equations
(5.4)–(5.14). The source's own 143-piece inner integration table is unnecessary
for this weaker, sufficient arithmetic estimate.

## 1. Result and why the weaker estimate is enough

Put α=3/40, λ=37/40 and use the source's normalizer m_{K,M}, with M fixed
and K tending to infinity through multiples of 40. Then Proposition 5.2 implies

\[
 \limsup K^{-2}\log m_{K,1000}
 \le \frac{439979}{320000}=1.374934375<\frac{11}{8}.
\]

Together with the separately certified analytic upper bound −69/50=−1.38,
this leaves a strict negative margin greater than 1/200 in the normalized
determinant exponent. The analytic estimate and its transfer to discrete
determinants are separate dependencies; see
[`fauzan-equilibrium-certificate.md`](fauzan-equilibrium-certificate.md).
The source's condition K≥200M² now reads K≥200,000,000. This large threshold
does not affect an asymptotic irrationality criterion, but makes this estimate
unsuitable for practical small approximants.

The motivation is simple: retain the cancellation of the main quadratic terms
in the prime-density envelope, but bound its remaining variance term uniformly.
The surviving fractional-part averages can be integrated exactly using Euler's
constant. A crude elementary lower bound for that constant is already enough.

## 2. Algebraic cancellation before integration

Write H=1+2α, f={x}, g={αx}, τ={2Hx}, σ={2x}, η={2λx}. Set

\[
 \ell(x,z)=\lfloor x-z\rfloor+\lfloor x+z\rfloor+1,\quad
 A(z)=\ell(x,z)-2x,\quad B(z)=\ell(\alpha x,z)-2\alpha x.
\]

Integrals in this section are over 0<z<1/2. The function ℓ(x,z) takes only
the values q=⌊2x⌋ and q+1; the latter occurs on a set of length σ/2.
Consequently ∫A=0 and ∫A²=σ(1−σ)/2≤1/8. The same statements hold for B.

Expanding the source definitions of Γ and N directly gives

\[
\begin{aligned}
\Gamma(x)&=10\alpha\lambda x^2-5\lambda x
 +\frac{-\tau^2+\tau\sigma+(\tau-\sigma)_+}{2}
 +3\int AB-9\int B^2,\\
N(x)&=-10\alpha\lambda x^2+\lambda x-2\lambda xf
 +12\lambda xg-\frac{\eta(1-\eta)}2.
\end{aligned}
\]

Thus R=−Γ−N=xF+Q, where

\[
 F=4\lambda+2\lambda\{x\}-12\lambda\{\alpha x\}\ge-8\lambda,
\]
\[
 Q=\frac{\tau(\tau-\sigma)-(\tau-\sigma)_+
              +\eta(1-\eta)}2+9\int B^2-3\int AB\le\frac{13}{8}.
\]

For the last inequality, τ(τ−σ)−(τ−σ)₊ lies in [−1/4,0],
η(1−η) lies in [0,1/4], and Cauchy–Schwarz gives |∫AB|≤1/8.
This verifies the source's decomposition and bound without its periodic
antiderivatives or inner integration table.

## 3. Integrate the fractional parts globally

For a positive integer n, summing the integrals over unit intervals gives

\[
 \int_n^\infty\frac{\{x\}}{x^2}\,dx=H_n-\log n-\gamma.
\]

For 0<y<1, splitting at 1 gives instead 1−log y−γ. Since 3α<1,

\[
\begin{aligned}
J:=\int_3^\infty F(x)x^{-2}\,dx
={}&\frac{4\lambda}{3}
 +2\lambda(H_3-\log3-\gamma)\\
 &-12\lambda\alpha(1-\log(3\alpha)-\gamma).
\end{aligned}
\]

The elementary bound γ>H_n−log n−1/(2n), at n=3, gives
γ>5/3−log3. For completeness the sequence H_n−log n−1/(2n)
increases to γ: its increment is the trapezoid approximation to
∫_n^{n+1}dt/t minus that integral, strictly positive by convexity.
The coefficient of γ in J is λ(12α−2)<0, so substitution yields

\[
 J<\lambda\left(\frac53+8\alpha+12\alpha\log\alpha\right).
\]

For fixed M≥40, F≥−8λ and Q≤13/8 imply

\[
 \int_3^M\frac{R(x)}{x^3}\,dx
 \le J+\frac{8\lambda}{M}+\frac{13}{144}.
\]

Here ∫_3^∞x⁻³dx=1/18. Applying source Proposition 5.2 gives

\[
 \limsup K^{-2}\log m_{K,M}
 \le I_{\rm out}+\frac{13}{144}
 +\lambda\left(\frac53+8\alpha+12\alpha\log\alpha\right)
 +\frac{14\lambda}{M}. \tag{*}
\]

The order of limits matters: M is chosen first and kept fixed. No uniform
prime-number-theorem estimate in growing M is used.

## 4. The outer integral without its 11-piece table

Use the source's R₀(y) and d(y), equations (5.8)–(5.9).
For 0<α<1/8 the identity min(α,v)=v−(v−α)₊ gives directly

\[
\begin{aligned}
\int_{1/3}^{1/2}R_0
 &=\frac{41}{72}-\frac{4\alpha}{3}
       +\frac54(\tfrac13-\alpha)^2-\frac{5\alpha^2}{6},\\
\int_{1/2}^{1}R_0
 &=\frac18+3(\tfrac12-\alpha)^2+\frac{5\alpha^2}{2}.
\end{aligned}
\]

Every term follows from ∫(c−dy)₊dy, evaluated at its single zero.
Moreover, on (1/3,1/2),

\[
 d(y)=(1+4\alpha-3y)_+-(1+\alpha-3y)_+,
 \qquad\int d=\frac{5\alpha^2}{2}.
\]

The floor contribution is −5λ/3, since ⌊1/y⌋ is 2 on (1/3,1/2),
1 on (1/2,1), and 0 thereafter. All five hinge zeros 2λ/j lie above 1/3,
so their integral up to 2λ is

\[
 \sum_{j=1}^5\frac{(2\lambda-j/3)^2}{2j}
 =2\lambda^2H_5-\frac{10\lambda}{3}+\frac56.
\]

Adding these expressions proves

\[
 I_{\rm out}=\frac{29}{12}-\frac{31\alpha}{6}
  +\frac{41\alpha^2}{12}-5\lambda+2\lambda^2H_5
  =\frac{127751}{96000}.
\]

## 5. A five-term rational logarithm certificate

Since α=(3/5)/8,

\[
 \log\alpha=-2\operatorname{atanh}(1/4)-6\operatorname{atanh}(1/3).
\]

The power series atanh z=∑_{k≥0}z^{2k+1}/(2k+1) has positive terms.
Keeping two terms in the first series and three in the second therefore gives

\[
 \log\alpha< -2\left(\frac14+\frac{(1/4)^3}{3}\right)
 -6\left(\frac13+\frac{(1/3)^3}{3}+\frac{(1/3)^5}{5}\right)
 =-\frac{33559}{12960}.
\]

Substitute this and M=1000 into (*). Exact fraction arithmetic gives

\[
 \frac{127751}{96000}+\frac{13}{144}
 +\frac{37}{40}\left(\frac53+\frac35
           -\frac9{10}\frac{33559}{12960}\right)
 +\frac{14(37/40)}{1000}
 =\frac{439979}{320000}=\frac{11}{8}-\frac{21}{320000}.
\]

Reproduce with `python codex/computations/check_table_free_arithmetic.py`.
No floating-point approximation enters this certificate.

## 6. What remains outside this simplification

The replacement needs the correctness of source (5.7), the outer local
valuation envelope, and the prime-sum conversion Proposition 5.2. Integrality
also needs its own Proposition 5.1 and underlying local lemmas. Our fresh
calculation verifies the subsequent envelope integration, not those premises.
The analytic certificate is stronger than needed for this arithmetic bound;
one should keep its equilibrium and discrete-energy justification separately
reviewable. No historical novelty is asserted for fractional-part integration,
the logarithm bound, or this way of combining them.

## 7. Optional second simplification: discard the refined outer rank saving

The separate local arithmetic audit by the `simplify` agent observes that one
can use the plain rank bound instead of the refined `min(r,z)` penalty. At the
envelope level, this means dropping −d(y) from the outer integrand. Its exact
cost is ∫d=9/640, already calculated above. This weakens the arithmetic estimate
but removes a further case distinction in the outer local argument.

Keep the same elementary calculation and take M=10000. The resulting bound is

\[
 A_{\rm plain}
 =\frac{439979}{320000}-\frac{14\lambda}{1000}
       +\frac{14\lambda}{10000}+\frac9{640}
 =\frac{2203747}{1600000}=1.377341875<\frac{69}{50}.
\]

The gap is 4253/1600000>1/400, so, subject to the local integrality and analytic
transfer hypotheses, the normalized determinant eventually satisfies
0<Q_{40n,10000}(ζ(5))≤exp(−4n²). This has ample decay relative to degree 37n.
The source's sufficient threshold becomes K≥20,000,000,000. Again this is an
asymptotic simplification, not an efficient finite computation. The script
checks the combined bound exactly. The justification for replacing the local
rank estimate is a separate lemma; adding ∫d does not by itself prove that lemma.
