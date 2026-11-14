# Independent review of the short elementary hybrid proof

**Verdict:** [the short proof](short-elementary-hybrid-proof.md) is a complete
standalone proof of \(\zeta(3)\)'s irrationality. Its arithmetic cancellation,
integral transformations, positivity, exact decay comparison, and final integer
contradiction are valid. A short absolute-convergence argument between the
two integrations by parts has been inserted directly into the proof.

No recurrence, prime number theorem, complex contour, or geometric divergence
certificate is needed. This assembles classical Apéry/Beukers/Zudilin
ingredients; no historical novelty is claimed.

## 1. Mellin normalization and the rational arithmetic

For
\[
L_n(x)=\frac1{n!}\frac{d^n}{dx^n}[x^n(1-x)^n],
\]
integration by parts gives
\(\int_0^1x^{t-1}L_n(x)\,dx=(-1)^nP_n(t)\).
The sign is important but disappears when two transforms are multiplied.
The displayed coefficient expansion makes the left side a rational function
of \(t\), so extension from \(t>n\) to every real \(t>0\) is legitimate.
There are no poles in that interval.

The residue at \(-k\) is
\[
c_k=(-1)^{n-k}\binom nk\binom{n+k}k.
\]
Squaring the partial fractions gives
\[
a_k=c_k^2,\qquad
\frac{b_k}{2}=c_k\sum_{j\ne k}\frac{c_j}{j-k}.
\]
Thus \(a_k\in\mathbb Z\) and \(D_nb_k/2\in\mathbb Z\). The latter bound
uses only the nonzero differences \(j-k\), of absolute value at most \(n\).
Since \(R_n=P_n^2=O(t^{-2})\), its simple-pole coefficients sum to zero.

Differentiation and summation give
\[
-\sum_{t\ge1}R_n'(t)
=2\sum_ka_k\zeta(3)+\sum_kb_k\zeta(2)
-\sum_k\left(2a_kH_k^{(3)}+b_kH_k^{(2)}\right).
\]
The \(\zeta(2)\) coefficient vanishes exactly. Every resulting rational
term is cleared by \(D_n^3\): use \(D_n^3H_k^{(3)}\in\mathbb Z\) and
\((D_nb_k/2)(D_n^2H_k^{(2)})\in\mathbb Z\).
The series is absolutely convergent because \(R_n'(t)=O(t^{-3})\).
At \(n=0\), \(I_0=2\zeta(3)\); at \(n=1\), the arithmetic gives
\(I_1=10\zeta(3)-12\), confirming the normalization.

## 2. The logarithm produces exactly the derivative

The Stieltjes identity
\[
\int_0^\infty\frac{dt}{(1+t)(v+t)}
=-\frac{\log v}{1-v}
\]
follows by an elementary partial-fraction antiderivative; at \(v=1\)
both sides equal \(1\).

Expanding the kernel at \(v=xy\) gives
\[
\sum_{t\ge1}
\iint L_n(x)L_n(y)(-\log(xy))(xy)^{t-1}\,dx\,dy.
\]
Each inner integral is the negative derivative of the product of Mellin
transforms and hence equals \(-R_n'(t)\). If
\(M_n=\max_{[0,1]}|L_n|\), the sum of absolute integrals is bounded by
\[
M_n^2\sum_{t\ge1}\iint(-\log(xy))(xy)^{t-1}\,dx\,dy
=2M_n^2\sum_{t\ge1}t^{-3}<\infty.
\]
This justifies the series, derivative, and integral interchanges without
assuming positivity of the Legendre factors.

## 3. The two integrations by parts and their improper limits

For fixed \(y,t>0\), the first integration gives
\[
\int_0^1\frac{L_n(x)}{xy+t}\,dx
=\int_0^1\frac{x^n(1-x)^ny^n}{(xy+t)^{n+1}}\,dx
=J_n(y,t)\ge0.
\]
The two factors \((-1)^n\), one from integration by parts and one from
differentiating the reciprocal kernel, cancel. Every endpoint term vanishes
because the first \(n-1\) derivatives of \(x^n(1-x)^n\) vanish at both endpoints.

One must justify reordering before the second integration, since \(L_n(y)\)
is still signed. The useful estimate is
\[
J_n(y,t)\le M_n\int_0^1\frac{dx}{xy+t}.
\]
Its integral against \(|L_n(y)|/(1+t)\) is at most
\(2M_n^2\zeta(3)\). Hence the first transformed triple integral is absolutely
convergent.

Now \(t=yu\) has Jacobian \(y\), cancelling the entire factor
\(y^n/y^{n+1}\). The remaining \(y\)-kernel is \((1+yu)^{-1}\).
The second Rodrigues integration gives the positive factor
\(y^n(1-y)^nu^n/(1+yu)^{n+1}\).
Fubini is justified by the preceding absolute bound. Its inner positive
integral equals the previous signed inner integral and is therefore bounded
by the corresponding absolute integral. Tonelli then proves finiteness of
the resulting nonnegative triple integral.

Finally \(u=(1-z)/(yz)\), with \(|du/dz|=1/(yz^2)\), gives exactly
\[
I_n=\int_{(0,1)^3}
\frac{[x(1-x)y(1-y)z(1-z)]^n}
     {(1-z+xyz)^{n+1}}\,dx\,dy\,dz.
\]
No factorial, sign, or power of \(y\) is missing. The integrand is strictly
positive in the open cube, so \(I_n>0\). When \(n=0\), no integrations by
parts are required and the same convergence argument applies.

## 4. The odds substitution and the exact minimum

With \(x=X/(1+X)\), and similarly for \(y,z\), direct cancellation gives
\[
\frac{x(1-x)y(1-y)z(1-z)}{1-z+xyz}=\Lambda^{-1},
\qquad
\frac{dx\,dy\,dz}{1-z+xyz}
=\Lambda^{-1}\frac{dX\,dY\,dZ}{XYZ}.
\]
These two identities prove the reciprocal-period normalization for every
\(n\), not only for its exponential part.

Put \(r=\sqrt{XY}\), \(A=\sqrt{(1+X)(1+Y)}\).
First minimizing in \(Z\), then using \(A\ge1+r\), gives
\[
\Lambda\ge(A+A^2/r)^2
\ge(3+2r+1/r)^2\ge(3+2\sqrt2)^2.
\]
The function \(A+A^2/r\) is increasing in \(A>0\), so the middle
substitution preserves the inequality. Equality is simultaneous at
\(X=Y=1/\sqrt2\), \(Z=1+\sqrt2\).

The already finite integral \(I_0\) therefore supplies the domination
\[
0<I_n\le I_0(17+12\sqrt2)^{-n}.
\]
The additional constant-term identity mentioned in the proof is not a
logical input to this estimate or to the final contradiction.

## 5. The lcm rate is elementary and strictly sufficient

For real \(x\ge0\), factorial prime-power counting gives the displayed
four-factorial sum
\[
T(x)=\sum_{p,a\ge1}h(x/p^a)\log p.
\]
Floors inside factorials cause no ambiguity:
\(\lfloor\lfloor x/m\rfloor/p^a\rfloor=\lfloor x/(mp^a)\rfloor\).
The function \(h\) has period \(6\), with values \(0,1,1,0,0,1\) on
the consecutive unit intervals. Hence it is nonnegative and equals \(1\)
on \([1,3)\), including the required endpoint \(1\). This proves
\(\psi(x)-\psi(x/3)\le T(x)\).

The estimate
\(\log(m!)=m\log m-m+O(\log(m+2))\) follows from upper and lower integrals
of the increasing function \(\log x\); no full Stirling expansion is needed.
In the factorial combination the \(x\log x\) and \(x\) terms cancel,
leaving
\[
T(x)=\left(\frac{\log2}{3}+\frac{\log3}{2}\right)x
+O(\log(x+2)).
\]
Telescoping at \(x,x/3,x/9,\ldots\) gives
\[
\log D_n\le
\left(\frac{\log2}{2}+\frac{3\log3}{4}\right)n
+O(\log^2(n+2)).
\]
The comparison with the analytic rate is exact:
\(\sqrt2>7/5\) and
\[
64\cdot3^9\cdot5^4=787320000
<815730721=169^4.
\]
Thus \(\exp(3c)<169/5<17+12\sqrt2\).
The strict linear exponential margin absorbs the \(O(\log^2n)\) term,
proving \(D_n^3I_n\to0\).

Under \(\zeta(3)=a/b\), the quantity
\(bD_n^3I_n=2(aD_n^3A_n-bD_n^3B_n)\)
is a positive integer. This closes the contradiction; positivity was proved
independently of the rationality assumption.

Exact rational checks for \(0\le n\le12\) additionally confirm the pole
coefficients, cancellation, \(D_n^3B_n\) integrality and initial normalization.
The period-six table, strict integer comparison, and odds Jacobian were also
checked exactly. These finite checks supplement the all-\(n\) arguments above.
