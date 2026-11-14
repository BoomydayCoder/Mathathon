# The hyperbolic system is a confluent Wilson path with an explicit Bernoulli transform

Status: direct identities and endpoint analysis proved below; independent
review pending. This verifies the root agent's Wilson-path proposal and adds
an exact polynomial transform. The factorial formula for p(0) and the ratio
p(−1)/p(0) remain conjectural for m>1. No irrationality or novelty claim.

Fix m≥1. For z>0 define

\[
 \Omega_{m,z}(y)=\frac{y^{2m+1}\coth(\pi y)}{\sinh^2(\pi y)+z},
 \qquad \mathcal T_m[p](z)=\int_0^\infty p(y^2)\Omega_{m,z}(y)\,dy.
\]

The components in the multiple-orthogonality problem are
ωⱼ(y)=y^(2m+1)coth(πy)csch(πy)^(2j), 1≤j≤m.

## 1. Exact Wilson normalization

For 0<a<1 put z=sin²(πa). Gamma reflection and recurrence give

\[
 \boxed{\left|
 \frac{\Gamma(a+iy)\Gamma(1-a+iy)\Gamma(iy)\Gamma(1+iy)}
      {\Gamma(2iy)}\right|^2
 =\frac{4\pi^3 y\coth(\pi y)}{\sinh^2(\pi y)+\sin^2(\pi a)}.} \tag{1}
\]

Indeed the first two gamma factors have squared modulus
π²/[sinh²(πy)+sin²(πa)]; the quotient involving the other factors has
squared modulus4πy coth(πy). Thus Ω_(m,z) is tᵐ/(4π³) times the ordinary
Wilson weight with parameters(a,1−a,0,1), t=y². Equivalently, multiplying
once by t changes the zero gamma parameter to1, so it is t^(m−1)/(4π³)
times the Wilson weight(a,1−a,1,1). The latter description has no zero
parameter for0<a<1 and uses only a polynomial change of the measure.

The parameter sum is fixed (2 in the first description,3 in the second).
Both a and1−a move. This is not the single-variable parameter family in
Beckermann–Coussement–Van Assche, arXiv:math/0303093, Theorem3.3. An appeal
to that theorem alone therefore supplies no formula for the present
confluent family.

## 2. The confluence is rigorous but is not an analytic Taylor series

For every polynomial p and0≤r≤m−1, differentiation under the integral,
including a one-sided limit at z=0, gives

\[
 \boxed{\mathcal T_m[p]^{(r)}(0)
 =(-1)^r r!\int_0^\infty p(y^2)\omega_{r+1}(y)\,dy.} \tag{2}
\]

Near zero the absolute integrand after r derivatives is bounded by a
constant times y^(2m)/(y²+z)^(r+1), hence by y^(2m−2r−2). This is
integrable for r≤m−1. At infinity there is a common exponential majorant
for a polynomial p. Dominated convergence proves(2).

Consequently the balanced(n,...,n) conditions are exactly

\[
 \mathcal T_m[t^kp]^{(r)}(0)=0,
 \quad0\le k<n,\quad0\le r<m. \tag{3}
\]

These are m confluent parameter conditions at each of n polynomial powers.
The separate triangular change of weights in
[hyperbolic-multiple-orthogonality.md](hyperbolic-multiple-orthogonality.md)
turns their moment system into rational coefficients.

The next derivative need not exist at the endpoint. In particular the
unmodified weight with two zero gamma parameters at a=0 behaves as y⁻².
The factor tᵐ repairs the required integrals but does not justify arbitrarily
many parameter derivatives or an unrestricted Wilson continuation.

## 3. Every monomial moment has a closed Bernoulli-polynomial formula

For every integer ℓ≥0 and0<a<1,

\[
 \boxed{\int_0^\infty
 \frac{y^{2\ell+1}\coth(\pi y)}{\sinh^2(\pi y)+\sin^2(\pi a)}\,dy
 =\frac{(-1)^{\ell+1}[B_{2\ell+2}(a)-B_{2\ell+2}]}
             {2(\ell+1)\sin^2(\pi a)}.} \tag{4}
\]

Here Bᵣ(a) is the Bernoulli polynomial with B₁=−1/2. The proof does not
require Wilson orthogonality. First use the elementary hyperbolic identity

\[
 \frac{z\coth(\pi y)}{\sinh^2(\pi y)+z}
 =\coth(\pi y)-\frac{\coth(\pi(y+ia))+\coth(\pi(y-ia))}{2}
 =2\sum_{r\ge1}e^{-2\pi r y}[1-\cos(2\pi r a)].
\]

The final summands are nonnegative, so Tonelli applies after multiplication
by y^(2ℓ+1). Integrate each exponential and use the classical even Bernoulli
Fourier series. Its coefficient follows by repeated integration by parts
from B₁(x)=x−1/2 on(0,1), and it converges absolutely here since its order
is2ℓ+2≥2. The resulting factorial ratio is exactly(4).

For example ℓ=0 gives a(1−a)/(2sin²πa). This also checks the endpoint
singularity of the base weight. The symmetry B_(2r)(1−a)=B_(2r)(a)
ensures that both solutions a,1−a of z=sin²πa give the same integral.

## 4. The antidifference Q is precisely the Wilson parameter transform

For an arbitrary polynomial p, let Q be the normalized antidifference

\[
 Q(x)-Q(x-1)=x^{2m+1}p(-x^2),\qquad Q(0)=0.
\]

The earlier construction proves Q(x)=Q(−x−1). Equivalently, if p(t)=Σpₑtᵉ,

\[
 Q(x)=\sum_e\frac{(-1)^e p_e}{2m+2e+2}
       [B_{2m+2e+2}(x+1)-B_{2m+2e+2}].
\]

Applying(4) term by term proves the exact spectral-parameter identity

\[
 \boxed{\mathcal T_m[p](\sin^2\pi a)
             =\frac{(-1)^{m+1}Q(a-1)}{\sin^2\pi a}.} \tag{5}
\]

Thus the Bernoulli-antidifference construction is the polynomial numerator
of a concrete Wilson-weight integral. It is not merely a device for solving
abstract rational moment equations.

If p satisfies the first m conditions (k=0), the independently proved
parameterization is Q(x)=[x(x+1)]^(2m+1)V(x(x+1)). With v=a(1−a), (5)
becomes

\[
 \boxed{\mathcal T_m[p](\sin^2\pi a)
       =\frac{(-1)^m v^{2m+1}V(-v)}{\sin^2\pi a}.} \tag{6}
\]

For general k, apply the same statement to tᵏp. A full balanced polynomial
is characterized by the divisibility of all n corresponding Bernoulli
numerators, not just the one at k=0.

## 5. The first fractional-power term remembers p(0)

For arbitrary p, the Taylor expansion in a of Q(−a) has no odd powers
below2m+1, since Q(x)−Q(−x)=x^(2m+1)p(−x²). Its coefficient at
 a^(2m+1) is−p(0)/2. Set a=arcsin(√z)/π for0<z<1. Its odd power
series and(5) give

\[
 \mathcal T_m[p](z)=
 \sum_{r=0}^{m-1}\frac{\mathcal T_m[p]^{(r)}(0)}{r!}z^r
 +\frac{(-1)^m p(0)}{2\pi^{2m+1}}z^{m-1/2}+O(z^m). \tag{7}
\]

The even powers of a contribute an analytic series in z; the displayed odd
power is the first possible fractional term. For p(0)≠0 this proves that
 the endpoint function is C^(m−1) but not Cᵐ. For a balanced polynomial,
all m displayed Taylor coefficients vanish, leaving the explicit branch
amplitude

\[
 \mathcal T_m[p](z)\sim
        \frac{(-1)^m p(0)}{2\pi^{2m+1}}z^{m-1/2}.
\]

For m=1 and p=t−1/3, the amplitude is1/(6π³), a direct normalization check.
The conjectural factorial expression for p_(m,n)(0) is therefore also a
conjecture about an explicit endpoint amplitude of this Wilson path. That
interpretation does not prove the factorial expression.

## 6. What remains open

The path and moment transform are exact for all m and polynomial degrees.
They do not yet give a closed balanced Rodrigues formula, an all-m proof of
p_(m,n)(0), or of p_(m,n)(−1)/p_(m,n)(0)=(n+1)^(2m). The last ratio was
independently observed by the valuations agent through m≤3,n≤6.

The complete balanced normality theorem is proved separately in
[hyperbolic-balanced-normality.md](hyperbolic-balanced-normality.md). It
ensures that the exact rational moment construction gives a growing family,
but a coefficient formula or arithmetic cancellation theorem is additional
work. The finite rational small-form experiments with these polynomials have
not established an irrationality argument.

Reproduce the gamma normalization, direct quadrature, exact Bernoulli
transform and fractional-power checks with
`python codex/computations/check_confluent_wilson_path.py`.
