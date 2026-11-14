# Independent audit of the higher-Rodrigues ζ(5) asymptotic

Reviewed: higher-rodrigues-laurent-geometry.md, principally Sections 4–6,
with its exact normalization traced back to
multiplicity-lifted-zeta-family.md. The earlier Laurent construction and
positive-period assertions were also checked.

**Verdict:** the stated \(m=2\) oscillatory asymptotic, sharp limsup, and
infinitely many terms of each sign follow from the displayed contour
argument. I found no unresolved mathematical gap. This is an audit of
these theorems, not a claim of an irrationality proof or historical novelty.
The all-\(m\) Fourier-contour extension in a separate later note is outside
this review.

One small missing explicit inequality was supplied during review:
\(1<r<3/2\) gives \(\lambda<10000<256^2\), proving algebraically that
\(|\eta|>1/256\). This justifies absorption of the real-saddle term without
relying on the displayed decimals.

## 1. Exact normalization and contour shift

For \(n\ge1\), the underlying rational function is

\[
F_n(z)=\binom{2n}{n}
 \frac{\prod_{j=1}^n(z-j)^2}{\prod_{k=0}^{2n}(z+k)}.
\]

It is exactly the gamma quotient in the draft. Its square \(R_n\) has
fourth-order zeros at \(1,\dots,n\), so the cotangent-derived contour
integrand has removable singularities there.

The original sum is \(-\tfrac16\sum_{k>n}R_n'''(k)\).
Closing an upward line to its right has clockwise orientation.
Three integrations by parts change sign three times, producing

\[
\mathscr L_n=-\frac1{12\pi i}\int
                  (\pi\cot\pi z)'''F_n(z)^2\,dz.
\]

Reflection gives

\[
F_n(z)^2=\binom{2n}{n}^2\frac{\sin^4\pi z}{\pi^4}
 \frac{\Gamma(z)^6\Gamma(n+1-z)^4}{\Gamma(z+2n+1)^2}.
\]

The sign from \(\sin\pi(z-n)\) disappears in its fourth power.
Also

\[
(\pi\cot\pi z)'''\frac{\sin^4\pi z}{\pi^4}
=-4-2\cos2\pi z.
\]

These identities give exactly the factor \(1/(6\pi i)\) in the reflected
integral. With the draft's \(T_{n,j}\) normalization,

\[
\mathscr L_n=\tfrac23T_{n,0}
                  +\tfrac16(T_{n,1}+T_{n,-1})
             =\tfrac23T_{n,0}+\tfrac13\Re T_{n,-1}.
\]

There are no poles of the reflected integrand in \(0<\Re z<n+1\).
For each fixed \(n\), uniformly on a closed horizontal interval inside
that strip, gamma asymptotics give

\[
H_n(x+iY)=O_n(|Y|^{-2}e^{-4\pi|Y|}).
\]

Multiplication by either Fourier factor or by the cosine loses at most
\(e^{2\pi|Y|}\). Thus horizontal sides tend to zero exponentially and
any two vertical lines inside the strip give the same integral. The
line may subsequently be chosen as \(c=an\); the contour equality holds
for each \(n\) before asymptotics are taken.

## 2. Stirling's formula is uniform on the whole chosen line

Write \(t=a+iy\), with fixed \(0<a<1\). Use
\(\Gamma(u+1)=u\Gamma(u)\) first. The three gamma arguments after this
rewrite are \(nt\), \(n(1-t)\), and \(n(t+2)\). Their real parts are
bounded below by a positive multiple of \(n\), and they remain in
\(|\arg u|\le\pi/2\). The logarithmic Stirling remainder is therefore
\(O(1/n)\) uniformly for all real \(y\), not just on compact intervals.
Exponentiating preserves a uniform relative \(1+O(1/n)\) error.

Independently counting the constants and powers gives

\[
H_n(nt)=(2\pi)^4n^{-2}
  \frac{(1-t)^2}{t^3(t+2)}e^{n\psi(t)}(1+O(n^{-1})).
\]

The central binomial square contributes \(16^n/(\pi n)\), so the draft's
factor \(16\pi^3n^{-3}\) is correct. All logarithms used here are analytic
principal logarithms: \(t\), \(1-t\), and \(t+2\) have positive real part.
The remaining amplitude \(f(t)\) is bounded on the line and is
\(O(|y|^{-2})\) at infinity; in particular its absolute value is integrable.

## 3. The full saddle, not just the maximum of the modulus

Direct differentiation gives

\[
\psi'(t)=6\log t-4\log(1-t)-2\log(t+2),\qquad
\psi''(t)=6/t+4/(1-t)-2/(t+2).
\]

I independently combined its real part over a common denominator:

\[
\Re\psi''(a+iy)=
\frac{12((3a+1)y^2+a(1-a)(a+2))}
 {(y^2+a^2)(y^2+(1-a)^2)(y^2+(a+2)^2)}>0.
\]

Thus \(\Im\psi'(a+iy)\) is strictly increasing in \(y\), from \(-4\pi\)
to \(4\pi\). For
\(g_j(y)=\Re(\psi(a+iy)+2\pi ij(a+iy))\),

\[
g_j'(y)=-\Im\psi'(a+iy)-2\pi j.
\]

For each \(j=0,\pm1\) this derivative has exactly one zero; its limiting
slopes at both tails have the required signs and nonzero magnitude.
This establishes a unique global maximum on each chosen vertical line.

The draft additionally chooses a line where the **complete complex
derivative** vanishes. For \(j=0\), \(t=2/3\) has \(\psi'(t)=0\).
For \(j=-1\), the upper root \(t_\star\) of \(2t^3-3t+2=0\) has
\(\Re t_\star=r/2\in(1/2,1)\). Since

\[
(1-t)^2(t+2)=t^3-3t+2,\qquad
\exp(\psi'(t)/2)=\frac{t^3}{(1-t)^2(t+2)},
\]

the algebraic equation gives the value \(-1\). Its modulus forces
\(\Re\psi'(t_\star)=0\). The strict range
\(0<\Im\psi'(t_\star)<4\pi\), furnished by the previous monotonicity,
then uniquely gives \(\psi'(t_\star)=2\pi i\).
There is no unresolved \(2\pi i\) ambiguity.

Finally \(\Re\psi''(t_\star)>0\), so the vertical quadratic expansion is
a decaying complex Gaussian. A curved steepest-descent path is not needed.

## 4. Uniform tails and the relative error

The global maximum argument suffices to make the saddle approximation
rigorous. Choose a fixed small neighborhood of its \(y\)-coordinate.
Strict monotonicity of \(g_j'\), continuity on compact sets, and its
nonzero limiting tail slopes imply a number \(\delta>0\) such that
\(g_j(y)\le g_j(y_\star)-\delta\) everywhere outside this neighborhood.
The integrable bound on \(|f(a+iy)|\), together with the uniform Stirling
error, bounds the full omitted integral by

\[
O\bigl(n^{-2}e^{n(g_j(y_\star)+\log16-\delta)}\bigr).
\]

This is exponentially smaller than the main Gaussian term.
Inside the neighborhood put \(y-y_\star=v/\sqrt n\).
The leading Gaussian is \(\exp(-\psi''(t_\star)v^2/2)\), with strictly
positive real quadratic coefficient. Its integral is
\(\sqrt{2\pi}/\sqrt{\psi''(t_\star)}\), using the root with positive
real part. The first Taylor corrections of order \(n^{-1/2}\) are odd
polynomials in \(v\) times an even Gaussian and integrate to zero on a
symmetric interval. The remaining Taylor terms and the uniform Stirling
error give relative \(O(1/n)\). Gaussian domination justifies these
integrations and extension of the local interval.

The orientation and measure are \(dz=in\,dy\), so the total leading
constant is

\[
8\pi^2\sqrt{2\pi}\,
\frac{f(t_\star)}{\sqrt{\psi''(t_\star)}}.
\]

For the real saddle, \(f(2/3)=9/64\) and
\(\psi''(2/3)=81/4\), giving exactly
\(\pi^{5/2}/(2\sqrt2)\).

## 5. Critical values, sign changes and sharp limsup

At a full saddle the exponent can be written without logarithmic
guesswork:

\[
\log16+\psi(t)+2\pi ijt
=\log16+4\log(1-t)-4\log(t+2)
  +t(\psi'(t)+2\pi ij).
\]

The last term vanishes. Thus the selected bases are \(1/256\) and
\(\eta=16((1-t_\star)/(t_\star+2))^4\), as stated.
The latter is the Laurent critical value for \(x=-t_\star\).
The exact cubic for the three plus-branch values has product one and
one positive real root \(\lambda\), giving \(|\eta|=\lambda^{-1/2}\).
The additional bound \(\lambda<10000<256^2\) proves that the real-saddle
term can be included in the smaller error.

The resulting assertion is therefore valid:

\[
\boxed{\mathscr L_n=
\frac{n^{-5/2}}3\Re(\mathcal C(t_\star)\eta^n)
 +O(n^{-7/2}|\eta|^n).}
\]

Here \(\mathcal C(t_\star)\ne0\): none of its amplitude factors vanishes,
and the second derivative has positive real part.

The sign argument does not need an unproved irrationality of
\(\arg\eta/\pi\). In the irrational-rotation case, density provides
subsequences on which the real part is bounded away from zero with either
sign. In the rational case the orbit has at least three points because
\(\eta\) is nonreal. Its sum is zero, and its real parts after multiplication
by a nonzero \(\mathcal C\) cannot all be zero: a nonzero real linear
functional vanishes on only two antipodal points of the unit circle.
Both signs therefore occur periodically. On either selected subsequence
the error is smaller than the main term by a factor tending to zero.

This proves infinitely many positive and infinitely many negative terms,
and also

\[
\limsup_n|\mathscr L_n|^{1/n}=|\eta|.
\]

It does not assert eventual nonvanishing at every index or a uniform
relative asymptotic after taking the real part near its zeros.

## 6. Checks of the preceding Laurent and positive-period statements

The constant-term extraction follows by matching the two \(z\)-exponents
and then independently matching each \(x_i\)-exponent. Differentiating
\(\log\Lambda_m\) gives \(z^2=C/P\) and
\(x_i=mz/(1+z)\), so the critical equations and values in the draft
follow, with all forbidden zero factors excluded.

On positive log coordinates the first softplus terms have positive
diagonal Hessian; the remaining softplus composition is convex.
The stated coercive bound follows from
\[
\max(w,0)+\max(S_--w,0)\ge(S_-+|w|)/2.
\]
This proves unique minimization and integrability of the inverse period.
I also differentiated the log-coordinate formula independently and
recovered the Hessian blocks and their determinant.

The support contains the zero vector and every positive coordinate unit
vector, so equality in the torus triangle inequality occurs only at
zero angle. This supplies the needed aperiodicity for the unique torus
Gaussian. For the inverse real period, the integrable \(n=0\) density
and a strict gap away from the unique minimum control the noncompact
tails. Thus the two positive asymptotics in Section 3 have no hidden
saddle multiplicity.

## 7. Independent computations

The new script check_higher_rodrigues_independent.py does not merely
compare the exact form with the claimed leading cosine. It:

1. Reconstructs the exact rational form from the double and simple
   partial-fraction coefficients of \(F_n^2\), checking cancellation of
   the total simple coefficient and the normalization at \(n=1\).
2. Checks the curvature numerator and algebraic branch equation exactly.
3. Integrates both Barnes components directly after removing their
   exponential scale, at \(n=1,3,10,30\).

The direct Barnes combination agrees with the independently reconstructed
rational form to relative discrepancies about \(1.5\cdot10^{-61}\)
at \(n=1\) and \(3.6\cdot10^{-50}\) at \(n=3\), at 65-digit working
precision. The quantities \(n(T_{n,j}/T_{n,j}^{\rm leading}-1)\) are
bounded in these tests; for \(j=-1\) they range from approximately
\(-0.633-0.037i\) at \(n=1\) to \(-1.058-0.138i\) at \(n=30\).
These quadratures are corroborating numerical evidence, not certified
enclosures. The all-\(n\) conclusions rest on Sections 1–5 above.

The available denominator bound still does not yield an irrationality
criterion. This sharper analytic theorem supplies a rigorously selected
oscillating cycle; it does not establish the actual growth of primitive
denominators.

## 8. Optional algebraic strengthening: each sign has density one half

The main proof above deliberately does not require irrationality of the
angle. Nevertheless that property can be proved, giving a stronger
corollary. The following argument was independently checked by the author
of the contour note after this reviewer proposed it.

Let \(K\) be the splitting field of \(f(x)=2x^3-3x-2\).
The polynomial

\[
f(x+1)=2x^3+6x^2+3x-3
\]

is Eisenstein at three, so \(f\) is irreducible. Its discriminant is
\(-216\), a nonsquare. Hence \(\operatorname{Gal}(K/\mathbb Q)=S_3\),
and its unique quadratic subfield is \(\mathbb Q(\sqrt{-6})\).

Every root of unity in \(K\) generates an abelian Galois extension of
\(\mathbb Q\). The corresponding quotient of \(S_3\) is abelian, so it
has order at most two. The quadratic cyclotomic fields are
\(\mathbb Q(i)\) and \(\mathbb Q(\sqrt{-3})\), neither equal to
\(\mathbb Q(\sqrt{-6})\). Thus the only roots of unity in \(K\) are
\(1,-1\).

The values \(\eta,\overline\eta\) are rational functions of roots of
\(f\), so their ratio belongs to \(K\). It is not one because the values
are nonreal. It is not minus one either: that would make the pair purely
imaginary, whose contribution to the second elementary symmetric sum of
\(\lambda,\eta,\overline\eta\) would be
\(\eta\overline\eta>0\). But their value polynomial gives that symmetric
sum as \(-17376/256<0\). Therefore

\[
\boxed{\eta/\overline\eta\text{ is not a root of unity},
\qquad \arg\eta/\pi\notin\mathbb Q.}
\]

Write \(\theta=\arg\eta\), \(\phi=\arg\mathcal C(t_\star)\).
The proved asymptotic becomes

\[
\frac{3n^{5/2}\mathscr L_n}
 {|\mathcal C(t_\star)|\,|\eta|^n}
=\cos(n\theta+\phi)+O(1/n).
\]

Irrational rotations are equidistributed. For every fixed
\(0<\epsilon<1\), all sufficiently large \(n\) with
\(\cos(n\theta+\phi)>\epsilon\) have positive \(\mathscr L_n\), and
positive \(\mathscr L_n\) requires
\(\cos(n\theta+\phi)>-\epsilon\). The corresponding limiting frequencies
are \(\arccos(\epsilon)/\pi\) and
\(\arccos(-\epsilon)/\pi\). Letting \(\epsilon\) decrease to zero gives
natural density one half for positive terms, and likewise for negative
terms. Zero terms are confined eventually to
\(|\cos(n\theta+\phi)|\le\epsilon\), whose density tends to zero with
\(\epsilon\).

Thus the positive and negative subsequences each have natural density
\(1/2\), and the zero subsequence has density zero. This uses a standard
equidistribution theorem and a small Galois argument; the earlier
infinitely-many-signs conclusion needs neither strengthening.
