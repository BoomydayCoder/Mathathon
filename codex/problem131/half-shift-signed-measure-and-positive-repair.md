# Half-shifted odd zeta: an exact positivity barrier and a rational repair

Status: theorem and proof written from the root agent's proposal, independently
checked by the integral agent. Exact certificates and numerical integral
checks are supplied. A further independent review is welcome. No historical
novelty or irrationality conclusion is asserted.

For odd s=2m+1>=3 define, initially at a>0,

\[
 G_s(a^2)=(s-1)a^{s-3}\zeta(s,a+1/2).
 \tag{1}
\]

The half shift removes the elementary endpoint corrections of the earlier
Hermite transform. It also changes positivity in a decisive way.

## 1. The pure transform is positive-Stieltjes only at s=3

For every odd s>=3 there is an exact signed-density representation

\[
 \boxed{G_s(z)=\int_0^\infty\frac{v_s(y)}{z+y^2}\,dy,\qquad
 v_s(y)=-\frac{y^{2m-1}}{(2m-1)!}
                  \frac{d^{2m}}{dy^{2m}}\tanh(\pi y).}
 \tag{2}
\]

It is valid throughout C\(-infinity,0], using the square root a with
positive real part in (1). Every absolute polynomial moment of v_s is finite.
There is no additional polynomial term or atom hidden in (2).

At s=3,

\[
 v_3(y)=2\pi^2y\tanh(\pi y)\operatorname{sech}(\pi y)^2>0.
\]

For s>=5, G_s(x)>0 for x>0 but G_s(x) tends to zero as x decreases to zero.
It therefore cannot be the Stieltjes transform of any nonzero positive
measure on [0,infinity): such a transform is decreasing and its limit at
zero is at least its strictly positive value at any fixed positive point.
The same argument excludes generalized positive Stieltjes representations
with nonnegative constant or atomic terms. This is a positivity obstruction
for this pure transform, not an obstruction to other zeta constructions.

**Stronger finite-correction obstruction.** For every odd s>=5, adding any
real rational function R(z) to G_s(z) still cannot produce a positive
Stieltjes transform. To prove this, first note that v_s must be negative
somewhere: otherwise (2) would be a nonzero positive Stieltjes representation,
contradicting the zero limit at the origin. By continuity it is negative on
an open y interval. Choose a point in that interval whose negative square
is not a pole of R; there are only finitely many excluded points. At that
point equation (3) gives

\[
 \operatorname{Im}G_{s,+}(-y^2)=-\frac{\pi v_s(y)}{2y}>0.
\]

A real rational function has real boundary value there, so the same positive
imaginary boundary value persists for G_s+R and, by continuity, at nearby
points in the upper half-plane. But every positive Stieltjes transform has
nonpositive imaginary part throughout that half-plane. This is impossible.
Any poles of R off the negative real axis would already violate analyticity
of a positive Stieltjes transform. Finite atoms or elementary rational
endpoint corrections therefore cannot remove the negative continuous part.
This does not exclude changing the transform or mixing another nonrational
Stieltjes function, as in the positive repair below.

## 2. Reflection, Cauchy reconstruction, and endpoint control

The cotangent partial-fraction identity, differentiated 2m times, gives

\[
 \zeta(2m+1,w)-\zeta(2m+1,1-w)
 =\frac{1}{(2m)!}\frac{d^{2m}}{dw^{2m}}[\pi\cot(\pi w)].
\]

At w=1/2+iy, cot(pi w)=-i tanh(pi y) and partial_w=-i partial_y. Since
(iy)^(2m-2)=(-1)^(m-1)y^(2m-2), the two boundary values of (1) satisfy

\[
 G_{s,+}(-y^2)-G_{s,-}(-y^2)
 =\frac{i\pi y^{2m-2}}{(2m-1)!}D_y^{2m}\tanh(\pi y)
 =-\frac{\pi i}{y}v_s(y).
 \tag{3}
\]

This is exactly the Stieltjes jump of the density v_s(y)dy: in the variable
x=y² its density is v_s(sqrt x)/(2sqrt x). The sign and the factor of two
in this change of variable are essential.

For completeness, the large-circle estimate can be obtained without an
unproved uniform Hurwitz asymptotic. The elementary Laplace representation is

\[
 \zeta(s,a+1/2)=\frac1{(s-1)!}\int_0^\infty
       \frac{t^{s-1}}{2\sinh(t/2)}e^{-at}\,dt.
 \tag{4}
\]

The function g_s(t)=t^(s-1)/(2sinh(t/2)) extends smoothly to zero and all
its derivatives are integrable at infinity. Repeated integration by parts
in (4) is therefore uniform for Re a>=0 as |a| tends to infinity. Its
first nonzero Taylor term is t^(s-2); the next term has degree s. This gives
G_s(a²)=a^(-2)+O(|a|^(-4)) uniformly in the closed right half-plane.
At a=0 the Hurwitz factor in (1) is analytic and finite, so G_s(z) is bounded
near z=0, and is O(|z|^(m-1)) for m>=2.

These estimates justify a keyhole Cauchy contour: the large circle vanishes,
the small circle at zero vanishes, and the cut jump is (3). The latter is
absolutely integrable. Indeed D_y^(2m)tanh(pi y)=O(y) at zero and decays
exponentially at infinity, so v_s(y)=O(y^(2m)) at zero and has exponential
decay times a polynomial at infinity. Cauchy reconstruction gives (2)
exactly. Boundedness at zero excludes an atom there; decay at infinity
excludes an additional polynomial or constant. No off-cut poles occur
because Re(a+1/2)>0.

## 3. All moments and explicit negative Gram certificates

The moment formula is

\[
 \boxed{M_{s,k}:=\int_0^\infty y^{2k}v_s(y)\,dy
 =\binom{s+2k-2}{2k}(-1)^k B_{2k}(1/2),\quad k\ge0.}
 \tag{5}
\]

In particular every moment is rational and M_(s,0)=1. To prove (5), expand

\[
 \frac{t}{2\sinh(t/2)}
 =\sum_{k\ge0}\frac{B_{2k}(1/2)}{(2k)!}t^{2k}
\]

at zero in (4), with any fixed finite number of terms. Repeated integration
by parts supplies the corresponding asymptotic with a controlled remainder.
The coefficient of z^(-k-1) in G_s is
binom(s+2k-2,2k)B_(2k)(1/2). Comparing with the Stieltjes expansion from
(2), whose signed density has all absolute moments finite, gives (5).

The first variance is

\[
 M_{s,2}-M_{s,1}^2
 =\frac{s(s-1)(-3s^2+31s+14)}{5760}<0\quad(s\ge11).
 \tag{6}
\]

Thus t-M_(s,1) already has negative squared norm for every odd s>=11.
The smaller cases also have exact polynomial certificates. For s=7 and 9,
respectively, the integer polynomials

\[
 1960t^2-13520t+16457,\qquad
 280t^2-4796t+11693
\]

have squared norms -1300901/4 and -19814465/4 under (5). For s=5, the
integer polynomial with coefficients in ascending order

\[
 (15737,-44451,34814,-10217,1242,-62,1)
\]

has squared norm

\[
 -683248827068443/1761607680.
\]

These are exact finite certificates that the rational moment functionals
are not positive, independent of the monotonicity argument in section 1.

## 4. A small positive repair for zeta(5)

Put u=pi y. Direct differentiation yields

\[
 v_5(y)=\frac43\pi^4y^3\tanh u\operatorname{sech}^2u
                           (1-3\operatorname{sech}^2u),
\]

and therefore

\[
 \frac{v_5(y)}{v_3(y)}=\frac23u^2(1-3\operatorname{sech}^2u).
\]

The density v5+v3/3 is strictly positive. Multiplication by the positive
quantity 3cosh²u/v3 reduces this claim to

\[
 (1+2u^2)\cosh^2u-6u^2>0.
\]

Since cosh²u>=1+u²+u⁴/3, writing v=u² gives the explicit certificate

\[
 (1+2u^2)\cosh^2u-6u^2
 \ge1-3v+\frac73v^2+\frac23v^3
 =\frac73(v-9/14)^2+\frac1{28}+\frac23v^3>0.
 \tag{7}
\]

Consequently

\[
 \boxed{4a^2\zeta(5,a+1/2)+\frac23\zeta(3,a+1/2)
  =\int_0^\infty\frac{v_5(y)+v_3(y)/3}{a^2+y^2}\,dy}
 \tag{8}
\]

is a positive Stieltjes transform with rational moments. Its total mass is
4/3; multiplication by 3/4 normalizes it to mass one if desired. At a=1/2
its target value is zeta(5)+(2/3)zeta(3), so positivity has been restored
by mixing in the lower odd zeta value. It has not isolated zeta(5).

## 5. Every fixed higher odd index admits a finite rational repair

There is a simple constructive generalization. Let a_(m,j) be the positive
integer coefficients already proved in the hyperbolic decomposition note.
The tanh derivative has the same coefficients with alternating signs:

\[
 -D_u^{2m}\tanh u
 =2\tanh u\sum_{j=1}^m(-1)^{j-1}a_{m,j}\operatorname{sech}^{2j}u.
\]

This follows from the directly checked identity

\[
 D_u^2[\tanh u\operatorname{sech}^{2j}u]
 =4j^2\tanh u\operatorname{sech}^{2j}u
 -2(j+1)(2j+1)\tanh u\operatorname{sech}^{2j+2}u.
\]

Set r=m-1>=1. Then

\[
 \frac{v_{2m+1}}{v_3}
 =\frac{u^{2r}}{(2m-1)!}
   \sum_{j=1}^m(-1)^{j-1}a_{m,j}\operatorname{sech}^{2j-2}u.
\]

The constant term is a_(m,1)=4^r>0. For every l>=1,

\[
 u^{2r}\operatorname{sech}^{2l}u
 \le4^l u^{2r}e^{-2lu}
 \le4^l\left(\frac{r}{le}\right)^{2r}
 <4^l\left(\frac{r}{2l}\right)^{2r}.
\]

Thus the explicit positive rational number

\[
 B_m=\frac1{(2m-1)!}
 \sum_{\substack{1\le l\le r\\l\text{ odd}}}
 a_{m,l+1}4^l\left(\frac{r}{2l}\right)^{2r}
 \tag{9}
\]

satisfies v_(2m+1)+B_m v3>0 for every y>0. All its moments are rational,
and its Stieltjes transform is G_(2m+1)+B_m G3. For example this coarse
construction gives B2=2 and B3=8; the special zeta(5) repair 1/3 is much
smaller. No optimality is claimed. This proves existence of a finite rational
positive repair for every fixed odd s>=5 using only the lower s=3 density.

The repair settles positivity, not arithmetic isolation. Rational moments
and positivity permit standard Gram and Padé constructions for these mixed
targets, but eliminating the lower-zeta term while retaining useful
positivity and denominator control is a separate unresolved problem.

## Reproduction and scope

[check_half_shift_positive_repair.py](../computations/check_half_shift_positive_repair.py)
checks the signed derivative recurrence, moment coefficients, all negative
Gram certificates, the exact positive square certificate, and the rational
general repair constants. Independent positive/signed integral quadrature
checks (2) and (8) at several parameters. Numerical checks corroborate the
analytic proof; no numerical positivity search is used in the proof.
