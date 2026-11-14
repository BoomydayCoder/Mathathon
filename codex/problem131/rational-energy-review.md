# Independent review of the standard-library energy certificate

Status: independently audited `certify_equilibrium_energy_rational.py`, tested
its new primitives against separate exact rational calculations, and reran the
complete certificate. No mathematical or rounding error was found. The result
removes the Arb dependency from this scalar energy certificate; it does not
replace the separate equilibrium existence, global inequality, or arithmetic
arguments in the full construction.

The certified objective enclosure is
\[
[-1.385510140104,\,-1.385509421992],
\]
where the printed decimal endpoints are rounded outward. In particular the
objective is strictly less than
\[
-\frac{2771}{2000}=-1.3855<-\frac{69}{50}.
\]
The accompanying JSON stores the exact interval endpoints as integers divided
by \(2^{120}\). Neither the comparison nor the certificate relies on decimal
floating-point values.

## 1. Arithmetic foundation and scope of the endpoint box

The script imports the previously independently audited dyadic interval and
Taylor-jet primitives from `certify_equilibrium_endpoints.py`. Integer addition,
min/max products, floor/ceiling division and integer square roots all round
outward. It does not alter that source file. The only new primitives are the
logarithm and its jet recurrence.

The parameter intervals are exactly
\[
a\in[0.00007119508,0.00007119517],\qquad
b\in[0.75394957,0.75394967],
\]
with \(\alpha=3/40\) and \(\lambda=37/40\). All subsequent estimates are
uniform over this box. The endpoint certificate separately proves that the
actual solution lies in it. The energy formulas involving the endpoint equations
are needed only at that solution: evaluating their algebraic expressions over
the entire box encloses the true value, even though arbitrary other points in
the box need not define an equilibrium measure.

## 2. Logarithm: range reduction, monotonicity and remainder

For \(1\le x\le2\), put \(z=(x-1)/(x+1)\), so \(0\le z\le1/3\). Then
\[
\log x=2\sum_{k=0}^{L-1}\frac{z^{2k+1}}{2k+1}+R_L,
\qquad
0\le R_L\le\frac{2z^{2L+1}}{(2L+1)(1-z^2)}.
\]
The script uses \(L=42\). Its loop leaves `term` equal to the interval
containing \(z^{85}\), and its denominator is \(85(1-z^2)\), so the remainder
index is correct. Positivity allows addition of the one-sided error interval
from zero to the computed upper bound.

For a positive dyadic endpoint \(x=\texttt{raw}/2^{120}\), the bit-length
formula chooses an integer \(e\) so that \(x/2^e\in[1,2)\). Outward rounding
of the reduced interval preserves inclusion in \([1,2]\). The formula
\(\log x=\log(x/2^e)+e\log2\) is valid also when \(e<0\); ordinary interval
multiplication reverses the endpoints automatically in that case.

Finally, logarithm is increasing. Evaluating its lower input endpoint and its
upper input endpoint separately, then keeping their respective outer bounds,
therefore gives a valid enclosure of the full input interval. Every call has
strictly positive lower endpoint. The elapsed-time counter and its printed
rounding have no role in any bound.

## 3. Logarithmic jets and Simpson's constant

The jet coefficients represent \(a_k=a^{(k)}/k!\). If
\(c=\log a\), the identity \(a'=ac'\) gives
\[
c_k=\frac{k a_k-\sum_{j=1}^{k-1}j c_j a_{k-j}}{k a_0},
\qquad k\ge1.
\]
This is exactly the implemented recurrence. With a whole cell supplied as the
zeroth interval coefficient, interval operations enclose the derivative values
at every point of that cell, uniformly in \(a,b\). The parameter intervals
are constant jets, so no spurious differentiation with respect to parameters
is introduced.

For a cell of length \(h\), Simpson's error is at most
\(h^5\sup|f^{(4)}|/2880\). Since the fourth jet coefficient is
\(f^{(4)}/24\), the code's divisor is \(2880/24=120\), as required.
The three point evaluations are rounded outwards, and the symmetric derivative
error is added on every complete cell. The panels cover exactly
\([0,\alpha]\), \([\alpha,1]\), and the transformed tail \([0,1]\).

## 4. The removable tail and its uniform remainder

Write \(\delta=b-a\),
\[
z_0(v)=\frac{\delta}
{(\sqrt{1+av^2}+\sqrt{1+bv^2})^2},\qquad
\tau_s=\frac{\delta}{(\sqrt{s+a}+\sqrt{s+b})^2}.
\]
For \(0\le v\le1\), one has \(0<z_0\le\delta/4\) and
\(0<\tau_s<1\). Under \(u=1/v\), a single logarithmic tail term is
\[
\frac{\log(1-\tau_s z_0(v)v^2)}{v^2}
=-\tau_s z_0(v)H(\tau_s z_0(v)v^2),
\quad
H(z)=\sum_{k\ge0}\frac{z^k}{k+1}.
\]
This is analytic at zero with the stated limiting value.

Despite the constant's name `TAIL_TERMS=14`, Horner evaluation includes powers
\(z^0\) through \(z^{14}\), namely 15 coefficients, ending in \(z^{14}/15\).
Its remainder is therefore
\[
0\le H(z)-H_{14}(z)\le\frac{z^{15}}{16(1-z)}.
\]
Put \(B=\delta/4\). Multiplication by \(\tau_s z_0\le B\) gives the
uniform bound \(B^{16}/[16(1-B)]\) for a single tail term. The absolute values
of the two coefficients sum to \(1+6\alpha\). The code's `TAIL_ERROR` is
exactly
\[
(1+6\alpha)\frac{B^{16}}{16(1-B)}.
\]
The transformed integration interval has length one, so this also bounds the
integral error. Simpson is applied only to the finite polynomial expression;
no derivative bound for the omitted infinite tail is needed. Its separate
uniform value error is added afterwards. The final factor 2 for the charge
on the tail is applied to the whole enclosed integral, including this error.

## 5. Energy assembly and normalization

The code is the same potential/virial expression independently derived in
`equilibrium-independent-review.md`, with the elliptic and hypergeometric
numerical calls replaced by rational bounds. Define
\[
D_s=\frac{(\sqrt{s+a}+\sqrt{s+b})^2}{4},\qquad
f(u)=\log(1-\tau_1\tau_{u^2})
      -6\alpha\log(1-\tau_{\alpha^2}\tau_{u^2}).
\]
The charge integral is
\[
\int q(u)f(u)du
=3\int_0^1 f(u)du-6\int_0^\alpha f(u)du
  +2\int_1^\infty f(u)du.
\]
Thus the combination `3*medium-6*small+2*large` has the correct signs and
normalization. Reuse of the `small` interval can only widen the resulting bound
through lost correlation; it cannot invalidate inclusion.

The virial identity and off-support potential formulas then give
\[
I_V=2\lambda^2+\lambda(\log D_1-6\alpha\log D_{\alpha^2})
 +\lambda(-2+12\alpha)+\int qf.
\]
The plateau expression is
\[
M=2\lambda\log((b-a)/4)-V_0+kappa(1)
 -6\alpha\kappa(\alpha)-4\lambda,
\quad V_0=-2+12\alpha(1-\log\alpha),
\]
where \(\kappa(u)=\log(4u^2/(\sqrt{u^2+a}+\sqrt{u^2+b})^2)\).
Finally the code evaluates
\[
\mathcal E=\frac{\lambda M-I_V}{2}+C_*,\quad
C_*=-2\lambda+12\alpha\lambda(1-\log\alpha)
       +3\lambda^2-2\lambda^2\log(2\lambda).
\]
All factors agree with both previously audited energy routes.

## 6. Independent verification and cross-check

Run:

```
python codex/computations/check_rational_energy_primitives.py
python codex/computations/certify_equilibrium_energy_rational.py
```

The first script uses only the standard library. It compares logarithms with
separate 100-term exact rational series bounds at 69 endpoints across many
scales and at 20 intervals. It checks 80 logarithm jets through the independent
reciprocal-series formula for \(a'/a\), and 21 tail-series remainders. All pass.

The independently rerun full certificate returns:

| Quantity | Outward decimal enclosure |
|---|---|
| Regular integral on \([0,\alpha]\) | [0.039641077630, 0.039641142678] |
| Regular integral on \([\alpha,1]\) | [0.118993681963, 0.118993813006] |
| Transformed tail | [0.041802114458, 0.041802148657] |
| \(I_V\) | [1.915868444229, 1.915869572215] |
| \(\mathcal E\) | [−1.385510140104, −1.385509421992] |

These intervals lie inside the displayed intervals from the previous Arb
potential/virial certificate and agree with the Carlson certificate. Agreement
is a cross-check, not the reason the rational certificate is rigorous: the
outward arithmetic and proved remainder bounds give that conclusion directly.

The final integer comparison uses the lower dyadic enclosure of the negative
rational threshold. Thus `objective.hi < I(-2771/2000).lo` is slightly stronger
than the desired strict rational inequality, never weaker.

Conditionally combining this separately certified analytic constant with the
parent's arithmetic upper bound \(2203747/1600000\) leaves the exact margin
\[
\frac{2771}{2000}-\frac{2203747}{1600000}
=\frac{13053}{1600000}>0.
\]
At \(K=40n\) this is \(13.053n^2\) before lower-order errors, so a bound by
\(\exp(-13n^2)\) holds eventually if those separate arithmetic and asymptotic
hypotheses are supplied. This final comparison does not independently audit
those hypotheses.
