# Independent review: the positive logarithm-to-zeta deformation

**Status: accepted, including the fully real finite-Chebyshev resolvent
proof.** This reviews
[logarithm-zeta3-wilson-deformation.md](logarithm-zeta3-wilson-deformation.md)
and its latest checker. No missing endpoint mass, interchange hypothesis,
normalization factor, or arithmetic overclaim was found.
The underlying Wilson family is classical. Historical priority of the
specific reconstruction is not established here.

## 1. Density, endpoints, and mass

For real \(|\tau|\le1/2\), the prefactor
\((\sin(\pi\tau)/(\pi\tau))^2\) lies in \((0,1]\), interpreted continuously
at zero. Also
\(\cosh^2(\pi y)-\sin^2(\pi\tau)\ge\sinh^2(\pi y)\).
Thus the common bound
\[
 0\le w_\tau(y)\le\frac{4\pi^2y}{\sinh(2\pi y)}
\]
is correct. It is bounded near zero and integrable with every polynomial
weight at infinity. The boundary density is exactly
\[
 w_{\pm1/2}(y)=\frac{16y}{\sinh(2\pi y)}.
\]
Dominated convergence holds for all polynomial moments and for every fixed
resolvent with \(a>0\). In particular the endpoint limit occurs in total
variation for these absolutely continuous measures, so no atom or escaping
mass can appear. The independent moment computation below gives mass one,
including at both endpoints and at \(\tau=0\).

## 2. Finite Chebyshev products and the real resolvent proof

For \(x>0\), the numbers
\(\cos(2\pi(k+1/2+\tau)/N)\), \(0\le k<N\), are the roots of
\(T_N(z)+\cos(2\pi\tau)\). They are distinct at generic \(\tau\).
Dividing its monic factorization by \(2^N\) gives the product identity
(7a), with factor \(2^{1-2N}\). Continuity supplies the coincident-root
cases \(\tau=0,\pm1/2\). Logarithmic differentiation is legitimate for
\(x>0\), where every factor is positive.

For the limiting partial fractions, pair \(k\) with \(N-1-k\).
On the paired half of the index set, both angles corresponding to
\(b\pm\tau\), \(b=k+1/2\), lie in \((0,\pi/2]\) for
\(|\tau|<1/2\). The inequalities
\(\sin\theta\ge2\theta/\pi\) and
\(N\sinh(2x/N)\le\sinh(2x)\) give exactly
\[
 \frac{\sinh(2x/N)}{2N}
 \frac1{\sinh^2(x/N)+\sin^2(\pi(b\pm\tau)/N)}
 \le\frac{\sinh(2x)}{8(b\pm\tau)^2}.
\]
The right sides are summable for each fixed interior parameter.
For odd \(N\), the unpaired central sine is
\(\cos(\pi\tau/N)\), bounded away from zero, so its term is \(O_x(N^{-2})\).
Dominated convergence therefore proves the infinite real partial-fraction
formula. No uniform bound as \(\tau\to1/2\) is being silently used at this step.

The algebra converting its central difference into \(w_\tau/y\) is correct.
For the subsequent signed integral, direct differentiation gives
\[
 \partial_c^2\frac{y}{y^2+c^2}
 =\frac{2y(3c^2-y^2)}{(y^2+c^2)^3},\qquad
 \left|\partial_c^2h_c(y)\right|
 \le\frac{6y}{(y^2+c^2)^2}.
\]
Since \(h_a\le1/(2a)\), the integral of this upper bound against \(h_a\)
is at most \(3/(2ac^2)\). The triangular central-difference formula and
\(c\ge(1-2|\tau|)b\) now prove absolute summability of the integrated
series. Hence exchanging its sum and integral is justified.

The real integral
\(\int_0^\infty h_a(y)h_c(y)\,dy=\pi/[2(a+c)]\) includes \(a=c\) by
continuity. Substitution yields
\[
 S_\tau(a^2)=2\sum_{\ell\ge0}
 \frac1{(\ell+a+1/2)[(\ell+a+1/2)^2-\tau^2]}.
\]
The density majorant and a common summable bound for this series extend
the identity to zero and the boundary. For every fixed \(a>0\), its first
denominator remains positive even at \(|\tau|=1/2\).
The digamma and positive Laplace versions have the same signs and factors.
Thus the coalescence \(S_0=2\zeta(3,a+1/2)\) is proved by real arguments.

## 3. Logarithmic values, moments, and the Hilbert isometry

The changes of variables \(x=e^{-u}\), then \(x=v^q\), give the stated
factor \(q^3\). Differentiating the three proposed elementary primitives
verifies the exact logarithmic members:
\[
 S_{1/2}(1/4)=16\log2-8,\quad
 S_{1/3}(1/4)=27(\log3-1),\quad
 S_{1/4}(1/4)=96\log2-64.
\]
No special digamma evaluation is needed.

For nonzero interior \(\tau\), the integrated absolute Fourier series
is bounded by a constant times \(\sum r^{-2j-2}\), including \(j=0\).
The even Bernoulli Fourier series therefore gives
\[
 \nu_\tau(t^j)=
 \frac{(-1)^j}{(j+1)\tau^2}
 [B_{2j+2}(1/2+\tau)-B_{2j+2}(1/2)].
\]
This has the claimed removable limit and is a rational polynomial in
\(\sigma=\tau^2\) of degree \(j\). Its zeroth value is one.
The endpoint formulas follow from the common density majorant, rather
than an unjustified limiting interchange in a nonuniform Fourier bound.

The primitive-functional identity (13) follows on each monomial by
integrating \((x+1/2)g(-(x+1/2)^2)\) and using symmetry of the even
Bernoulli polynomial about \(1/2\). The earlier polynomial difference
certificate gives the proposed primitive involving \(Z'=E\), with the
factor \(1/[2(j+k+1)]\). Applying Bernoulli translation yields
\[
 \nu_\tau(D_jD_k)
 =\frac{d_j(\tau)d_k(\tau)}{j+k+1}.
\]
All \(d_k\) are positive on the entire parameter interval. Dividing by them
proves the Hilbert Gram identity and therefore the transported Legendre
orthogonality and norm. This is a polynomial proof, with no hidden Wilson
integral import.

## 4. Positive Padé error and its exact arithmetic scope

At the negative half-square, the signs in the Legendre coefficient and
Newton product cancel. Thus \(A_n(\sigma)\) has the stated positive
sum, and polynomial division makes \(N_n^\tau\) rational whenever
\(\sigma\in\mathbb Q\). The finite tail \(T_n(\tau)\) is rational under
the same hypothesis; \(\tau\) need not itself be rational.

The squared-error identity follows from orthogonality against
\((\Phi_n^\tau-A_n)/(t+(n+1/2)^2)\). Its factors \(1/2\) and \(A_n^{-2}\)
are correct. The exact norm gives the claimed upper bound, including \(n=0\).
For \(0\le\sigma\le1/4\),
\[
 0\le\log\prod_{r=1}^k(1-\sigma/r^2)^{-1}
 \le\frac{\sigma}{1-\sigma}\sum_{r\ge1}r^{-2}\le\frac23,
\]
which proves the uniform comparison with the Apéry denominator.
This transfers the stated sufficient real decay, not an integer
normalization theorem.

The unscaled integrality assertion is exactly as narrow as claimed:
\(A_1\in[5,19/3]\) leaves only the integer values \(5,6\), corresponding
to \(\sigma=0,1/5\). At \(1/5\), \(A_2=1774/19\) is not an integer.
Therefore simultaneous unscaled integrality already at indices one and
two forces \(\sigma=0\). This neither excludes other scalings nor bounds
their denominator growth.

## 5. Higher derivatives and the classical source

On a sufficiently small complex disk in \(\sigma\), the entire function
\(\sin^2(\pi\sqrt{\sigma})\) has modulus below \(1/2\).
The density consequently has a common integrable majorant with every
polynomial weight, permitting all parameter derivatives.
The resolvent series yields \(2\zeta(2m+3,a+1/2)\) at derivative order \(m\).
The moment polynomial has leading coefficient \((-1)^j/(j+1)\), so the
derived measure at order \(m\ge1\) has mass zero and a nonzero \(m\)-th
moment. It is a finite signed measure and must have both signs.
The displayed first derivative density agrees with direct differentiation.
This is an obstruction to retaining this positive ordinary Stieltjes
measure, not to other kernels or constructions.

Gamma reflection gives the Wilson weight with parameters
\((1/2,1/2,1/2+\tau,1/2-\tau)\).
Its conventional mass is
\(\Gamma(1+\tau)^2\Gamma(1-\tau)^2
=(\pi\tau/\sin\pi\tau)^2\), yielding exactly the normalization in (1).
The previously audited Ismail–Koelink dictionary gives
\(A=1,B=2,c=1,a_0=1-\tau^2,\lambda_0=1/4,\eta=1\pm\tau\),
hence \(\alpha=\beta=0,\delta=\pm\tau\). The source match and the
elementary quadratic-form identity are correct.
The direct endpoint proof avoids relying on a zero-parameter continuation
of the classical Wilson integral.

## 6. Reproducible checks

The author's checker passed its generic Gram identities, exact polynomial
families, logarithmic antiderivatives, finite products, resolvent and moment
quadratures, and positive error checks.
The [independent checker](../computations/check_logarithm_deformation_independent.py)
uses exact algebraic resultants for nine finite-product identities:
\[
 \operatorname{Res}_w(w^N+r,w^2-2zw+1)
 =2rT_N(z)+r^2+1.
\]
This is the Chebyshev factorization without numerical trigonometry.
It separately reconstructs Bernoulli moments using rational arithmetic and
verifies 144 Gram identities, including coalescence, the endpoint, and
\(\sigma=1/5\), plus the integrality distinction. All checks passed.
The domination and interchange proofs above, rather than finite checks,
establish the full analytic statements.
