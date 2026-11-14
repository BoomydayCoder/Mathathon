# Eulerian numbers split higher Rodrigues zeta forms into geometric periods

**Status:** exact all-m decomposition, saddle/critical-value identity, and
conditional component asymptotic proved below. Exact rational certificates
also establish the sharp oscillating asymptotic for m=2,3,4,5. The central
component has an unconditional asymptotic for every m. The entire higher-odd-zeta form is not
claimed to decay for every m. In fact the simplest vertical-contour argument
encounters a concrete obstruction already at m=6. This extends
[higher-rodrigues-laurent-geometry.md](higher-rodrigues-laurent-geometry.md).
No novelty claim is made.

## 1. The exact decomposition

Retain the minimal Mellin-zero polynomials Jₙ,ₘ and their pure forms

\[
 \mathscr L_{n,m}=2m\left(\sum_{k=0}^{mn}\binom{mn}{k}^2
                      \binom{n+k}{k}^{2m}\right)\zeta(2m+1)-B_{n,m}.
\]

Write κₙ,ₘ=(mn)!/(n!)^m. Their squared Mellin transform is

\[
 R_{n,m}(z)=\kappa_{n,m}^2
 \frac{\Gamma(z)^{2m+2}}{\Gamma(z-n)^{2m}\Gamma(z+mn+1)^2}.
\]

Reflection converts it to

\[
 R_{n,m}(z)=\kappa_{n,m}^2\frac{\sin^{2m}(\pi z)}{\pi^{2m}}H_{n,m}(z),
 \quad H_{n,m}(z)=
 \frac{\Gamma(z)^{2m+2}\Gamma(n+1-z)^{2m}}{\Gamma(z+mn+1)^2}. \tag{1}
\]

Let E(r,k) be the Eulerian number, counting permutations of r symbols with k
descents, and Aᵣ(w)=ΣₖE(r,k)wᵏ. The standard generating identity

\[
 \sum_{\ell\ge1}\ell^r w^\ell=\frac{wA_r(w)}{(1-w)^{r+1}}
\]

and w=exp(2πiz) give the exact meromorphic identity

\[
 \boxed{-\frac{\sin^{2m}(\pi z)}{\pi^{2m}}
       \left(\pi\cot\pi z\right)^{(2m-1)}
   =\sum_{j=1-m}^{m-1}E(2m-1,j+m-1)e^{2\pi ijz}.} \tag{2}
\]

It is enough to derive it for Im z>0, where the series converges; continuation
then gives the identity everywhere. In detail πcot πz=−iπ(1+2Σₗ≥₁wˡ),
while sin^(2m)(πz)=(−1)^m(1−w)^(2m)/(2^(2m)w^m). Differentiating the
former 2m−1 times and multiplying gives precisely the minus sign in (2).

The cotangent summation and odd integrations by parts yield

\[
 \mathscr L_{n,m}=-\frac1{2\pi i(2m-1)!}
 \int_{n-i\infty}^{n+i\infty}(\pi\cot\pi z)^{(2m-1)}R_{n,m}(z)\,dz.
\]

Every apparent singularity at an integer ≤n is removable, including the
intermediate integration-by-parts boundary products: R has a zero of order
2m, while the total derivative orders in a boundary product sum to 2m−2,
leaving a zero of order at least one after the cotangent pole is counted.
At imaginary infinity the rational factor and cotangent derivatives give
convergence. Using (1)–(2), define

\[
 T_{n,m,j}=\frac{\kappa_{n,m}^2}{2\pi i}
   \int_{c-i\infty}^{c+i\infty}H_{n,m}(z)e^{2\pi ijz}\,dz,
       \quad0<c<n+1.
\]

The individual integrands are holomorphic in that strip. Their gamma factors
have exponential modulus e^(−2mπ|Im z|) times a polynomial, whereas the
Fourier factor grows at most e^(2π|j||Im z|). Since |j|≤m−1, a uniform
exponential decay remains. Thus line shifts and the decomposition are valid:

\[
 \boxed{\mathscr L_{n,m}=
 \frac1{(2m-1)!}\sum_{j=1-m}^{m-1}
 E(2m-1,j+m-1)T_{n,m,j}.} \tag{3}
\]

The weights are positive and sum to one. Conjugation gives Tₙ,ₘ,₋ⱼ=conj Tₙ,ₘ,ⱼ.
For m=1 there is one component. For m=2 the weights are (1,4,1)/6.
For m=3 they are (1,26,66,26,1)/120.

This is an exact combinatorial bridge: repeated differentiation of the
cotangent kernel produces the descent distribution of permutations. Increasing
Mellin-zero multiplicity does not just change a scalar estimate; it introduces
new Fourier components, hence new candidate geometric saddle cycles.

## 2. General action and a global vertical maximum theorem

In the strip 0<Re t<1 define

\[
 \psi_m(t)=(2m+2)t\log t+2m(1-t)\log(1-t)
                    -2(t+m)\log(t+m),\quad
 f_m(t)=\frac{(1-t)^m}{t^{m+1}(t+m)}.
\]

For each fixed a∈(0,1), Stirling is uniform on t=a+iy and gives

\[
 \kappa_{n,m}^2H_{n,m}(nt)=
 m(2\pi)^{m+1}n^{-m-1}f_m(t)
 e^{n(2m\log m+\psi_m(t))}(1+O(n^{-1})). \tag{4}
\]

The prefactor fₘ is O(|y|⁻²) at infinity. As before, every gamma argument
stays in the closed right half-plane with positive real part proportional to
n, which makes the Stirling remainder uniform.

For every m≥1 and 0<a<1, direct algebra gives

\[
 \boxed{\frac d{dy}\Im\psi_m'(a+iy)=
 \frac{2m(m+1)\{(3a+m-1)y^2+a(1-a)(a+m)\}}
 {(y^2+a^2)(y^2+(1-a)^2)(y^2+(a+m)^2)}>0.} \tag{5}
\]

Its imaginary derivative increases from −2mπ to 2mπ. Therefore for each
integer |j|≤m−1 the real action
Re(ψₘ(a+iy)+2πij(a+iy)) has exactly one global maximum on that vertical line,
where Im ψₘ′=−2πj, and has linear negative tails away from the maximum.
This is a global statement about a vertical line, not a claim that its maximum
is a complete complex saddle.

**Conditional component theorem.** If tⱼ satisfies 0<Re tⱼ<1 and
ψₘ′(tⱼ)+2πij=0 with the stated principal logs, then

\[
 \boxed{T_{n,m,j}=\mathcal C_m(t_j)n^{-m-1/2}\eta_j^n(1+O(n^{-1})),}
 \tag{6}
\]

where

\[
 \mathcal C_m(t)=m(2\pi)^{m+1/2}
     \frac{f_m(t)}{\sqrt{\psi_m''(t)}},\quad
 \eta_j=\left(\frac{m(1-t_j)}{t_j+m}\right)^{2m}. \tag{7}
\]

The square root has positive real part; (5) gives Re ψₘ″>0 there.
To prove (6), choose a=Re tⱼ, apply (4), and expand about the unique maximum
in (5). The tail gap and |fₘ|=O(y⁻²) justify localization. The resulting
complex Gaussian is √(2π/(nψₘ″)); multiplication by n/(2π) from dz/(2πi)
gives exactly (7). The exponential equality follows by collecting the terms
multiplying t in the action and using the saddle equation. No enumeration of
critical values alone is used to choose a saddle.

## 3. The Laurent critical equation reappears exactly

Exponentiating half the saddle equation gives

\[
 t_j^{m+1}=(-1)^j(1-t_j)^m(t_j+m).
\]

With x=−tⱼ, this is precisely

\[
 (m-x)(1+x)^m=\varepsilon x^{m+1},\qquad
 \varepsilon=(-1)^{j+m+1}.
\]

Furthermore (7) is exactly the Laurent critical value
[m(1+x)/(m−x)]^(2m). Thus each actual saddle in the contour strip maps to
one of the previously classified critical points. Conversely an algebraic
root has to satisfy the correct logarithm branch and lie in the permissible
strip before theorem (6) can be invoked. This qualification is essential.

For every m the j=0 component satisfies the theorem unconditionally. There is
one t₀∈(0,1) solving

\[
 t_0^{m+1}=(1-t_0)^m(t_0+m).
\]

The left-to-right ratio increases strictly from zero to infinity, since its
logarithmic derivative is (m+1)/t+m/(1−t)−1/(t+m)>0. The saddle equation
then holds with real logarithms. Its value is strictly between zero and one,
since 0<m(1−t₀)/(t₀+m)<1. Therefore the central component always decays
exponentially, with a positive real leading coefficient given by (7).
It does not follow that the entire positive-weight sum (3) is controlled by
this component: the other T components are complex and may decay much more
slowly.

## 4. Certified finite extension: pure ζ(7), ζ(9), ζ(11) forms decay and oscillate

For m=2,3,4,5, every nonpositive-real-part Laurent critical point other than
the single positive real point gives a saddle in 0<Re t<1. There is one for
each j=1−m,…,m−1. The largest component modulus is attained exactly at the
outer conjugate pair j=±(m−1), and that modulus is strictly less than one.

These finite assertions are proved by the exact-rational certificate
`../computations/certify_rodrigues_finite_saddles.py`, with output retained in
`../computations/rodrigues_finite_saddles_certificate.json`. Here is the proof
method, including why numerical root finding is not the certificate:

1. Approximate roots merely propose centers c with rational real and imaginary
   parts, rounded to 10^(−10). On the circle |z−c|=r=10^(−7), exact rational
   arithmetic verifies
   \[
    |P(c)|_1+\sum_{k=2}^{\deg P}|P^{(k)}(c)|_1r^k/k!
       <\max(|\Re P'(c)|,|\Im P'(c)|)r,
   \]
   where |u+iv|₁=|u|+|v|. Rouché's theorem therefore gives exactly one root
   in the disk. Disjoint disks equal in number to the polynomial degree
   account for every root. A real-centered disk has a real root because its
   unique root is fixed by conjugation.
2. For every upper-half saddle disk write t=a+ib. Exact interval bounds for
   atan(b/a), atan(b/(1−a)) and atan(b/(a+m)) enclose Im ψₘ′. The atan bounds
   use a 50-term alternating series after reciprocal/π/4 reductions; π is
   enclosed by Machin's identity 16atan(1/5)−4atan(1/239). The resulting
   interval isolates exactly one integral multiple of 2π. Since the polynomial
   equation already forces Re ψₘ′=0 and Im ψₘ′∈2πℤ, this certifies the
   complete saddle equation, including its logarithm branch.
3. The squared ratio in (7) gives rational bounds
   \[
    |\eta_j|=\left(\frac{m^2((1-a)^2+b^2)}{(a+m)^2+b^2}\right)^m.
   \]
   The lower bound for the outer pair exceeds every other component's upper
   bound, and its upper bound is below one. Separate atan bounds place arg η
   strictly between successive multiples of π, certifying that η is nonreal.

The resulting safe decimal enclosures, widened from the exact rational output,
are

| m | zeta value | dominant modulus enclosure |
|---:|---:|---:|
| 2 | ζ(5) | (0.0112132149, 0.0112132465) |
| 3 | ζ(7) | (0.0173831918, 0.0173832381) |
| 4 | ζ(9) | (0.0585022717, 0.0585024153) |
| 5 | ζ(11) | (0.3399158979, 0.3399166785) |

Let t★ be the upper saddle with j=1−m and η★ its value. Since the two outer
Eulerian coefficients equal one, (3) and (6) now prove

\[
 \boxed{\mathscr L_{n,m}=
 \frac{2n^{-m-1/2}}{(2m-1)!}\Re\{\mathcal C_m(t_\star)\eta_\star^n\}
 +O\bigl(n^{-m-3/2}|\eta_\star|^n\bigr),\quad m=2,3,4,5.} \tag{8}
\]

Other components have a strict exponential gap and are absorbed in this error.
Consequently each of these pure forms is exponentially small, has exact
limsup root rate |η★|, and takes both signs infinitely often. The same
rotation argument as in the ζ(5) note proves the sign assertion without
requiring an irrational-angle theorem. This strengthens the previously
proved analytic family: the conclusion now includes ζ(7), ζ(9) and ζ(11).

It remains far from an irrationality proof. The established sufficient
normalizer Dₘₙ^(2m+1) has exponential rate m(2m+1), while the analytic decay
rates are only about 4.05, 2.84 and 1.08 for these three values. These are
comparisons of available bounds, not impossibility results for actual reduced
denominators or alternative constructions.

## 5. A concrete boundary for this method

For m=6, numerical diagnostics put the analogous outer upper saddle at
1.034353+1.306981i, with associated modulus about 3.01799. This is outside
Re t<1, so the present vertical-line theorem cannot be applied to it. Moving
the line through z=n+1 crosses actual poles of the reflected components.
A different contour analysis with residues would be required. The numerical
candidate modulus exceeding one does not prove that the actual pure form
grows, and no such conclusion is made.

The all-m results are the exact Eulerian decomposition, global vertical
maximum theorem, conditional component asymptotic, and central-component
decay. The stronger full-form asymptotic is certified only for m=2,3,4,5.
The Eulerian connection exposes the extra cycles and identifies precisely
what is missing before one could claim a general odd-zeta recipe.
