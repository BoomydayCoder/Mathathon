# A direct variational improvement: rank-one error rate below 1/500

Status: proved elementary trial bounds and exact trial-norm rate; independent
adversarial review requested. The universal comparison theorem has passed
the separate review rank-one-universal-rate-independent-review.md.
No adjacent-Apéry identity or general orthogonal-polynomial asymptotic
theorem is used. Historical novelty is not assessed.

The parent requested an independent variational improvement of the previous
\(1/16\) rate. The monomial trial below was derived from the integral.
The parent independently suggested exploring a repeated positive root;
the fixed rational choice below yields a short exact certificate.

## 1. Results

For the rank-one optimized errors
\[
E_{h,s}=\min_{\substack{p\ {\rm monic}\\\deg p=h-1}}
\int_0^\infty\frac{p(y^2)^2}{D_{2h+(s-3)/2}(y^2)}w_s(y)\,dy,
\]
every fixed odd \(s\ge3\) satisfies
\[
\boxed{\limsup_{h\to\infty}E_{h,s}^{1/h}
<\frac{25}{12544}<\frac1{500}.}
\tag{1}
\]

This improves the previous \(1/16\) bound by an explicit elementary trial.
It still does not identify the exact optimized rate or prove that this
optimized rate is a limit.

The simpler monomial trial \(p(t)=t^{h-1}\) has an exact norm rate
\[
\boxed{\rho_0=(1+v^2)^{-2}
=0.002226938498\ldots,\qquad v-\arctan v=\pi.}
\tag{2}
\]
The equation has a unique positive solution \(v=4.4934094579\ldots\).
The decimal displays are diagnostic; all inequalities used in (1) are
certified rationally below.

## 2. Scaling turns the trial into a one-variable potential

By the universal comparison it is enough to bound the \(s=3\) problem.
Write \(A=2\pi\) and, for a fixed shift \(a\ge0\), use
\[
p_{h,a}(t)=(t-a^2h^2)^{h-1},\qquad
T_h(a)=\int_0^\infty
\frac{(y^2-a^2h^2)^{2h-2}}{D_{2h}(y^2)}w_3(y)\,dy.
\]
These are admissible monic trial polynomials, so \(E_{h,3}\le T_h(a)\).
When \(a\) is rational their coefficients are rational too.

Put \(y=hx\), and define
\[
I(x)=\int_0^2\log(x^2+u^2)\,du
=2\log(x^2+4)-4+2x\arctan(2/x).
\]
Since \(\log(x^2+u^2)\) increases with \(u\), right Riemann sums give
the useful two-sided inequality, for every \(x>0\),
\[
\boxed{
e^{hI(x)}
\le\prod_{j=1}^{2h}\left(x^2+(j/h)^2\right)
\le\frac{x^2+4}{x^2}e^{hI(x)}.}
\tag{3}
\]
The upper error bound follows by telescoping the endpoint variation on
the \(2h\) cells; it is not an asymptotic assertion.

Define the trial potential
\[
\phi_a(x)=2\log|x^2-a^2|-I(x)-Ax,\qquad
M(a)=\sup_{x\ge0}e^{\phi_a(x)}.
\tag{4}
\]
At \(x=a\), interpret the exponential as zero. If \(a>0\), its value
at zero is its continuous limit; for \(a=0\) it is zero there.
At infinity it tends to zero, so \(0<M(a)<\infty\).

The terms in (4) have a direct meaning: a repeated polynomial root,
the uniform distribution of negative-square poles, and the exponential
tail of the positive density. No equilibrium theorem is required.

## 3. Uniform domination and the exact trial rate

The density has the exact form
\[
w_3(y)=A^2y^3e^{-Ay}\frac{1+e^{-Ay}}{(1-e^{-Ay})^3}.
\]
Using \(y/(1-e^{-Ay})\le y+A^{-1}\), one gets
\[
A^2y^3e^{-Ay}\le w_3(y)
\le2A^2(y+A^{-1})^3e^{-Ay}.
\tag{5}
\]

Combining the lower denominator bound in (3) with the upper density
bound gives, for \(h\ge2\),
\[
\begin{aligned}
T_h(a)
&\le2A^2\int_0^\infty
\frac{(x+(Ah)^{-1})^3}{(x^2-a^2)^2}
e^{h\phi_a(x)}\,dx\\
&\le C(a)M(a)^{h-2},
\end{aligned}
\tag{6}
\]
where one may take the finite explicit constant
\[
C(a)=\frac{A^2e^8}{128}
\int_0^\infty (x+A^{-1})^3(x^2+a^2)^2e^{-2Ax}\,dx.
\tag{7}
\]

To justify the apparent singularity in (6), after taking out
\(M(a)^{h-2}\) the remaining squared factor is
\((x^2-a^2)^2\), so it vanishes instead of diverging at \(x=a\).
Also \(I(x)\ge I(0)=4\log2-4\), giving the factor \(e^8/256\);
replacing \((x^2-a^2)^2\) by \((x^2+a^2)^2\) gives (7).
Thus the domination is uniform on the entire half-line.

The constant is a finite elementary sum:
\[
C(a)=\frac{A^2e^8}{128}
\sum_{i=0}^3\sum_{j=0}^2
\binom3i\binom2j A^{i-3}a^{4-2j}
\frac{(i+2j)!}{(2A)^{i+2j+1}}.
\]
For \(a=0\), it simplifies to \(687e^8/(2048A^6)\).

Conversely, the upper denominator bound and lower density bound imply
\[
T_h(a)\ge A^2\int_0^\infty
\frac{x^5}{(x^2+4)(x^2-a^2)^2}e^{h\phi_a(x)}\,dx.
\tag{8}
\]
Restrict to a compact interval of positive \(x\) away from \(a\), on
which \(\phi_a\) is within any prescribed error of its supremum.
Such an interval exists even if the maximum is at zero, by continuity
from the right. The prefactor has a positive lower bound there.
Equations (6) and (8) therefore prove
\[
\boxed{\lim_{h\to\infty}T_h(a)^{1/h}=M(a).}
\tag{9}
\]

This is an exact rate for a specified trial norm, not for the optimized
error.

## 4. Solve the monomial potential exactly

At \(a=0\),
\[
\phi_0'(x)=4/x-2\arctan(2/x)-2\pi.
\]
Set \(v=2/x\). The function \(v-\arctan v\) increases strictly from
zero to infinity. Hence the unique maximizer is \(x_*=2/v\), where
\(v-\arctan v=\pi\). Substituting the stationary equation cancels the
linear terms and gives
\[
e^{\phi_0(x_*)}
=\left(\frac{x_*^2}{x_*^2+4}\right)^2
=(1+v^2)^{-2}.
\]
This proves (2) by (9).

A simple strict bound is
\[
v>\sqrt{20},\qquad \rho_0<1/441.
\tag{10}
\]
Indeed \(\pi>157/50\) and
\[
\sqrt{20}+\arctan(1/\sqrt{20})
<21/\sqrt{20}<3\pi/2.
\]
The first inequality uses \(\arctan u<u\), and the second is verified
by squaring rational quantities after using the displayed lower bound
for \(\pi\). This says
\(\sqrt{20}-\arctan\sqrt{20}<\pi\), so strict monotonicity yields (10).

## 5. A repeated positive root gives a better elementary bound

Choose
\[
a=3/20,\qquad b=3/5.
\]
We bound \(e^{\phi_a}\) on three intervals.

**For \(0\le x<a\):** all three terms in the derivative have negative
sign, so the potential decreases. Its largest value is at zero:
\[
e^{\phi_a(0)}=\frac{a^4e^4}{16}
<\frac{a^4(11/4)^4}{16}<\frac{25}{12544}.
\tag{11}
\]

**For \(a<x\le b\):**
\[
e^{\phi_a(x)}
=e^{\phi_0(x)}(1-a^2/x^2)^2
\le\rho_0(1-a^2/b^2)^2
<\frac1{441}\left(\frac{15}{16}\right)^2
=\frac{25}{12544}.
\tag{12}
\]

**For \(x\ge b\):** the monomial potential is decreasing there, since
\(x_*=2/v<2/\sqrt{20}<b\). Thus
\[
e^{\phi_a(x)}\le e^{\phi_0(b)}
<\frac{81}{11881}\exp(-323/250)
<\frac{25}{12544}.
\tag{13}
\]
For the middle inequality, rewrite
\(\arctan(10/3)=\pi/2-\arctan(3/10)\), use
\(\arctan(3/10)<3/10\) and \(\pi>157/50\).
For the last inequality, the cubic partial sum of the positive
exponential series is already a sufficient lower bound on
\(\exp(323/250)\). Every remaining comparison is rational.

These three bounds show \(M(3/20)<25/12544<1/500\).
Equations (6) and (9) prove the claimed optimized-error bound for
\(s=3\). The universal comparison transfers it to every fixed odd
\(s\). More explicitly, for \(h\ge2\),
\[
E_{h,s}\le
(1+4\pi^2)^{(s-3)/2}\,C(3/20)
\left(\frac{25}{12544}\right)^{h-2}.
\]

## 6. Exact certification and limits

The standard-library checker check_rank_one_elementary_trial.py proves
all rational comparisons in (10)--(13). It obtains \(\pi>157/50\)
from Machin's identity and alternating-series bounds. It obtains
\(e<11/4\) from the ordinary exponential series and a geometric bound
on its tail. No floating-point value affects acceptance.

It separately encloses the monomial saddle by
\[
4.493409<v<4.493410,
\]
again using exact rational arctangent bounds. The rate interval printed
in decimals is diagnostic only.

The improvement concerns real decay. It provides no new theorem on
reduced rational denominators and no irrationality proof. The optimized
rate may be substantially smaller than the trial rate. No adjacent-Apéry
formula or claimed exact value for the optimized limit has been used.
