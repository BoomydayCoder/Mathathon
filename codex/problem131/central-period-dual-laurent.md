# Why the positive higher-zeta period has a different integer sequence

**Status:** exact identities proved below. This builds on the independent
beta-convolution observation in
[eulerian-periods-independent-review.md](eulerian-periods-independent-review.md).
The new coordinate transformation, dual constant-term sum and initial-value
calculation make the obstruction to simply retaining the positive component
concrete. No historical novelty claim is made.

**Primary-source context.** The underlying rational forms, finite Fourier/Barnes
decomposition and positive central integral are exact specializations of
Nesterenko’s published general constructions; the logarithmic polynomials
are confluent Jacobi–Piñeiro polynomials. See the explicit parameter dictionaries
in [rodrigues-primary-source-comparison.md](rodrigues-primary-source-comparison.md).
The specific geometric analysis and rate certificates here were independently
developed, but their historical priority remains unresolved.

## 1. Transform the positive central convolution into a second Laurent model

The central component from the Eulerian decomposition has the positive integral

\[
 T_{n,m,0}=\int_{(0,\infty)^{2m-1}}\int_0^1\int_0^1
 [(1-x)(1-y)]^{mn}
 \left(\frac{Uxy}{1+Uxy}\right)^{n+1}
 \prod_{i=1}^{2m-1}(1+u_i)^{-n-1}
 \frac{dx}{x}\frac{dy}{y}\prod_i\frac{du_i}{u_i},
 \quad U=\prod_i u_i.
\]

Use the odds substitutions x=a/(1+a), y=b/(1+b). Direct cancellation gives

\[
 \boxed{\widetilde\Lambda_m=
 (1+a)^m(1+b)^m\prod_{i=1}^{2m-1}(1+u_i)
 \frac{(1+a)(1+b)+Uab}{Uab},} \tag{1}
\]

\[
 \boxed{T_{n,m,0}=\int_{(0,\infty)^{2m+1}}
 \widetilde\Lambda_m^{-n-1}(1+a)^{m-1}(1+b)^{m-1}
 \frac{da}{a}\frac{db}{b}\prod_i\frac{du_i}{u_i}.} \tag{2}
\]

In particular the measure has a nonconstant positive weight when m>1.
The original beta-convolution proof supplies convergence. Formula (2) is not
asserting equality with the unweighted inverse period of the first Λₘ.

Expanding the final parenthesis in (1) and extracting every constant term gives

\[
 \boxed{b_{n,m}:=\operatorname{CT}\widetilde\Lambda_m^n
 =\sum_{k=0}^n\binom nk^{2m}\binom{mn+k}{k}^{2}.} \tag{3}
\]

Indeed its kth term has outer coefficient binom(n,k); each of the 2m−1
u-variables contributes another binom(n,k), while a,b each contribute
binom(mn+k,k). Compare this with the original integer coefficient

\[
 a_{n,m}=\sum_{k=0}^{mn}\binom{mn}k^2\binom{n+k}{k}^{2m}.
\]

The roles of the powers 2 and 2m and of the ranges n and mn have exchanged.
For m=1 the two sums and the two Laurent polynomials coincide. For m=2 the
new sequence begins b₀=1,b₁=10,b₂=626, whereas a₁=146. Thus retaining the
positive contour component changes its natural integer period as well.

## 2. Its positive minimum is exactly the central Barnes saddle

At a positive critical point let Q=Uab/[(1+a)(1+b)]. The logarithmic
uᵢ equations force every uᵢ=s=1/Q. The a and b equations then give

\[
 a=b=\frac{s}{m(1+s)},\qquad s^{m+1}=(m+1)s+m. \tag{4}
\]

There is one positive root s>1 of the last equation. Equivalently, put
t₀=s/(1+s); it is the unique root in (0,1) of

\[
 t_0^{m+1}=(1-t_0)^m(t_0+m),
\]

the central Barnes saddle equation. Substitution in (1) gives

\[
 \boxed{\min\widetilde\Lambda_m=
 \left(\frac{m+t_0}{m(1-t_0)}\right)^{2m}=\eta_0^{-1}.} \tag{5}
\]

To justify the word minimum globally, log of a positive Laurent polynomial
is convex in log coordinates. Here its exponent support contains 0 and every
positive coordinate unit vector, so its affine span is full and its Hessian
is positive definite everywhere. The displayed positive stationary point is
therefore the unique global minimum. Existence of a zero exponent mean with
strictly positive coefficients and full span also puts 0 inside the Newton
polytope, proving coercivity. Thus the usual torus Gaussian argument gives
lim bₙ,ₘ^(1/n)=η₀^(−1).

For m=2, t₀=2/3, s=2, a=b=1/3, and the minimum is 256. This explains the
real central scale 256⁻ⁿ by a positive Laurent model. It is different from
the first model's positive minimum 7953.13744… and from the pure form's
complex-pair modulus 0.01121323… . These are three distinct geometric
quantities, not conflicting estimates for one integral.

## 3. The positive initial value retains lower odd-zeta periods

At n=0 the central Barnes integral is well defined on 0<Re z<1 and simplifies
to

\[
 T_{0,m,0}=\frac1{2\pi i}\int
      \frac{\pi^{2m}}{z^2\sin^{2m}(\pi z)}\,dz. \tag{6}
\]

Close the upward line to the right using half-integer vertical edges. The
orientation is clockwise; the right edge tends to zero like its squared
real coordinate's reciprocal, and the horizontal edges decay exponentially.
Thus (6) is minus the sum of the residues at positive integers.

Write the finite relevant part of the even Taylor series as

\[
 \left(\frac w{\sin w}\right)^{2m}
 =\sum_{\ell\ge0}c_{m,\ell}w^{2\ell},\qquad c_{m,\ell}\in\mathbb Q,
 \quad c_{m,0}=1.
\]

At z=k, the residue equals

\[
 -\sum_{\ell=0}^{m-1}(2m-2\ell)c_{m,\ell}\pi^{2\ell}
                   k^{-(2m-2\ell+1)}.
\]

Summing proves the exact identity

\[
 \boxed{T_{0,m,0}=\sum_{\ell=0}^{m-1}
 (2m-2\ell)c_{m,\ell}\pi^{2\ell}\zeta(2m-2\ell+1).} \tag{7}
\]

For m=1 this is 2ζ(3), consistent with the original Beukers period. For m=2,
c₂,₁=2/3, so

\[
 \boxed{T_{0,2,0}=4\zeta(5)+\frac{4\pi^2}{3}\zeta(3).} \tag{8}
\]

For m=3, c₃,₁=1 and c₃,₂=8/15, giving

\[
 T_{0,3,0}=6\zeta(7)+4\pi^2\zeta(5)+\frac{16\pi^4}{15}\zeta(3).
\]

The Eulerian Fourier combination instead has initial value 2mζ(2m+1).
Already at n=0, its noncentral components must cancel these lower periods.
This is an explicit reason why positivity alone does not provide a pure
higher-odd-zeta irrationality form. The positive integral is valuable, but
its arithmetic target is different. A successful simplification would need
an additional cancellation mechanism that preserves both positivity and
sufficiently cheap rational denominators; neither has been supplied here.
