# Independent audit of the signed half-shift transform and positive repairs

**Verdict: accepted.** The signed resolvent, the pure-transform and finite
rational-correction obstructions, the moment certificates, and both
positive-repair constructions in
[the main note](half-shift-signed-measure-and-positive-repair.md)
are correct. No correction is required.

This audit also supplies a shorter independent derivation of the signed
measure and all its moments, using differentiation and integration by
parts instead of Cauchy reconstruction. It does not assert historical
novelty. The independent checker is
[check_half_shift_independent.py](../computations/check_half_shift_independent.py).

## 1. A direct alternative proof with no hidden atom or correction

The proved \(s=3\) resolvent from
[the elementary Padé proof](resolvent-legendre-pade-proof.md) is
\[
 R(a):=2\zeta(3,a+1/2)
   =-\int_0^\infty f''(y)h(a,y)\,dy,\qquad
 f(y)=\tanh(\pi y),\quad h(a,y)=\frac{y}{a^2+y^2}.
\]
Initially take \(a>0\). Put \(q=s-3=2m-2\), an even nonnegative
integer. Termwise Hurwitz differentiation gives
\[
 R^{(q)}(a)=(q+2)!\,\zeta(q+3,a+1/2),
 \qquad
 G_s(a^2)=\frac{a^q}{(q+1)!}R^{(q)}(a).
 \tag{A}
\]
The normalizing factorial is important: \(R\) already contains a
factor of two.

The kernel \(h\) is homogeneous of total degree \(-1\). If
\(\theta_y=y\partial_y\) and \(\theta_a=a\partial_a\), this gives
\((1+\theta_y)h=-\theta_a h\). The commuting Euler operators therefore
give the exact identity, valid for every integer \(q\ge0\),
\[
 \partial_y^q(y^q h)
 =\prod_{j=1}^q(\theta_y+j)h
 =(-1)^q\prod_{j=0}^{q-1}(\theta_a-j)h
 =(-a)^q\partial_a^q h.
 \tag{B}
\]
For our even \(q\), the sign is positive.

Differentiate the integral for \(R\), use (B), and integrate by parts
\(q\) times. The result is
\[
 G_s(a^2)
 =-\frac1{(q+1)!}\int_0^\infty
       \frac{y^{q+1}f^{(q+2)}(y)}{a^2+y^2}\,dy.
 \tag{C}
\]
This is exactly the stated density
\(-y^{2m-1}D_y^{2m}\tanh(\pi y)/(2m-1)!\).

All operations in this derivation are justified directly. On every
compact set of positive \(a\), derivatives of the rational kernel have
integrable majorants after multiplication by \(f''\).
At zero, \(y^q h=O(y^{q+1})\); every derivative of this factor appearing
in a boundary term still vanishes to order at least two.
At infinity, the derivatives of \(f\) appearing in boundary terms have
order at least two and decay exponentially, whereas the other factors
have at most polynomial growth. Thus every boundary term is zero.
For \(q=0\) the calculation is simply the original resolvent.

The density is absolutely integrable with every polynomial weight.
Both sides of (C), written in \(z=a^2\), are holomorphic on the slit
plane \(\mathbb C\setminus(-\infty,0]\). Equality on \(z>0\)
extends there by analytic continuation. This proves the complete
representation with no possible extra atom, polynomial, or rational
correction, independently of the contour argument.

The same integration by parts gives all moments:
\[
 M_{s,k}
 =\frac{(2k+s-2)!}{(s-2)!(2k+1)!}\,M_{3,k}.
 \tag{D}
\]
Here one transfers \(q=s-3\) derivatives from \(f^{(q+2)}\) onto
\(y^{2k+q+1}\). The boundary factors at zero again retain at least
two powers of \(y\), and all factors at infinity are exponentially
controlled. Substitution of
\(M_{3,k}=(-1)^k(2k+1)B_{2k}(1/2)\) gives precisely equation (5)
of the main note. In particular every signed measure has mass one.

## 2. Independent check of reflection and the uniform contour estimates

The original contour proof is also valid. In the odd Hurwitz reflection
formula there is a difference, not a sum. At \(w=1/2+iy\),
\(\cot(\pi w)=-i\tanh(\pi y)\) and
\(\partial_w=-i\partial_y\).
Combining this with
\((iy)^{2m-2}=(-1)^{m-1}y^{2m-2}\) gives
\[
 G_{s,+}(-y^2)-G_{s,-}(-y^2)
 =\frac{i\pi y^{2m-2}}{(2m-1)!}D_y^{2m}\tanh(\pi y)
 =-\frac{\pi i}{y}v_s(y).
 \tag{E}
\]
Thus the upper imaginary boundary value is
\(-\pi v_s(y)/(2y)\), as claimed. In the \(x=y^2\) variable the
density is \(v_s(\sqrt x)/(2\sqrt x)\); this accounts for the
otherwise easy-to-miss factor of two.

There is no unproved sector-uniform Hurwitz asymptotic. In the Laplace
integral let \(g_s(t)=t^{s-1}/(2\sinh(t/2))\).
It is smooth at zero, every derivative is integrable on the half-line,
and
\[
 g_s(t)=t^{s-2}-t^s/24+\cdots.
\]
Integrating its Laplace transform by parts \(s+2\) times gives a
remainder bounded by
\(\|g_s^{(s+2)}\|_1/|a|^{s+2}\) uniformly for
\(\operatorname{Re}a\ge0\). The endpoint terms of degrees below
\(s-2\) vanish; the degree \(s-1\) term is zero. After multiplication
by \(a^{s-3}/(s-2)!\), this proves
\[
 G_s(a^2)=a^{-2}+O(|a|^{-4})
\]
uniformly in the closed right half-plane, including its imaginary
boundary at large modulus.

Near \(a=0\) the Hurwitz factor is analytic and finite.
Consequently \(G_s(z)=O(z^{m-1})\) for \(m\ge2\), while \(G_3\)
is bounded. The density behaves as \(O(y^{2m})\) at zero and as a
polynomial times an exponential at infinity. These facts justify
the keyhole contour, endpoint limits and signed moment expansions
in the original argument.

## 3. The two positivity obstructions have the claimed scope

For \(s\ge5\), the transform is strictly positive on the positive axis
but tends to zero at the origin. Every nonzero positive Stieltjes
transform is decreasing on that axis, and its limit at zero is at least
any of its positive values. This immediately excludes a positive
representation for the pure \(G_s\).

The stronger rational-correction claim uses more than this real-axis
limit. Since the continuous signed density cannot be nonnegative
everywhere, it is strictly negative on an open interval.
At any point of that interval avoiding the finitely many poles of a
real rational \(R\), equation (E) gives a strictly positive upper
imaginary boundary value for \(G_s+R\). The rational function has
real boundary value there. Continuity from the upper half-plane
preserves this violation nearby, whereas a positive Stieltjes
transform has nonpositive imaginary part throughout that half-plane.
Off-cut poles of \(R\) would separately violate holomorphy.
Thus no finite real rational correction, including finite atoms,
can repair the pure transform.

This argument does not forbid the addition of another nonrational
transform. It therefore does not conflict with the positive repairs
in the same note. It also excludes only this specific normalized
half-shift transform, not all possible constructions for odd zeta
values.

## 4. Negative Gram witnesses and positive repairs

The exact moment formula gives
\[
 M_{s,2}-M_{s,1}^2
 =\frac{s(s-1)(-3s^2+31s+14)}{5760}.
\]
The quadratic factor is negative at \(s=11\) and decreasing for every
\(s\ge11\). Thus the linear witness works at every claimed odd index.
I independently evaluated the supplied quadratic witnesses for
\(s=7,9\) and degree-six witness for \(s=5\) using rational Bernoulli
numbers from their defining recurrence. All three displayed negative
norms agree exactly.

The ratio
\[
 \frac{v_5}{v_3}=\frac23u^2(1-3\operatorname{sech}^2u),
 \qquad u=\pi y,
\]
has the right sign and normalization. Multiplication of
\(v_5/v_3+1/3\) by \(3\cosh^2u\) gives
\((1+2u^2)\cosh^2u-6u^2\).
The lower Taylor bound for \(\cosh^2u\) is valid for every real \(u\)
because all omitted even coefficients are nonnegative.
The displayed square identity is exact and has positive constant
remainder \(1/28\). It proves strict positivity for all \(y>0\).
The repaired mass is \(4/3\), and at \(a=1/2\) its target is exactly
\(\zeta(5)+(2/3)\zeta(3)\).

For the general construction, differentiating
\(\tanh u\,\operatorname{sech}^{2j}u\) twice gives the stated
two-term recurrence. The positive coefficients obey
\[
 a_{m+1,j}=4j^2a_{m,j}+2j(2j-1)a_{m,j-1},
 \qquad a_{1,1}=1,
\]
so \(a_{m,1}=4^{m-1}>0\).
Only odd \(\ell=j-1\) give negative terms in the density ratio.
For each such term,
\[
 u^{2r}\operatorname{sech}^{2\ell}u
 \le 4^\ell u^{2r}e^{-2\ell u}
 \le4^\ell(r/(\ell e))^{2r}
 <4^\ell(r/(2\ell))^{2r}.
\]
The last strict inequality uses \(e>2\) and \(r\ge1\).
Summing these bounds gives the finite rational \(B_m\); all remaining
terms have positive sign. Thus the repair is strictly positive
throughout the stated range. The constants \(B_2=2,B_3=8\) check
exactly. This coarse general construction makes no optimality claim.

## 5. Computations and conclusion

The independent standard-library checker passes:

- 204 exact homogeneous derivative identities, including odd-order
  sign checks as well as the required even orders;
- 12 direct tanh derivative identities and 11 rational points with
  negative signed density;
- 280 integration-by-parts moment identities and the variance formula;
- all three explicit negative Gram certificates;
- the exact positive square identity and 11 rational repair constants;
- 65 exact positive Gram pivots for repaired measures, including ten
  consecutive pivots for the \(1/3\) repair.

The author's checker also passed its signed resolvent, reflection-jump
and repaired-resolvent quadratures. Finite checks support the proofs;
they are not extrapolated into all-degree positivity.

The established conclusion is a complete positivity barrier and
constructive positive repair for mixed targets. It gives no isolated
\(\zeta(5)\) irrationality proof. A later lower-bound construction
using an independent upper bound on \(\zeta(3)\) needs its own
inequality-direction and denominator audit.
