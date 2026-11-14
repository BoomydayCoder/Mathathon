# A genuine signed-measure model for ζ(7), and a quantitative obstruction

Status: independently constructed concrete model, elementary arithmetic
lemmas, and a proved obstruction to estimating its determinants by total
variation. Nonvanishing uses Proposition 5.4 of the primary paper below after
checking its hypotheses explicitly. No irrationality proof for ζ(7) is claimed.

Primary source read: Fischler–Lai–Sprang–Zudilin, *Many more odd zeta values are
irrational*, arXiv:2611.16156v1, 12 November 2026, §§2 and 5; local text
`codex/sources/fischler-lai-sprang-zudilin-2026.txt`. Their asymptotic theorem
assumes large s. We do **not** specialize that theorem to s=7; we apply only
their general signed-measure nonvanishing proposition, whose hypotheses are
verified below. The hand-picked denominator set {1,2,3} is not asserted to be
one of their totient-defined sets.

## 1. Elimination as a divided difference

Choose dilation weights

\[
 (w_1,w_2,w_3)=(45,-9,1).
\]

They satisfy

\[
 \sum_{b=1}^3w_b b^3=\sum_{b=1}^3w_b b^5=0,
 \qquad\sum_{b=1}^3w_b b^7=1080.
\]

This is Lagrange interpolation at the nodes 1,4,9:
w_b b³=1080/∏_{c≠b}(b²−c²). The divided difference annihilates constants
and linear functions of b² and takes the value one on (b²)². Thus the signs
and the coefficient 1080 are forced by a familiar interpolation construction.
Any nonzero weights annihilating ∑w_b b³ must change sign; positivity cannot
survive this particular method of eliminating the lower odd zeta values.

For even n≥2, set

\[
 R_n(t)=108^n(n!)^4
 \frac{(t+1/2)_n(t+1/3)_n(t+2/3)_n}{(t)_{n+1}^7}.
\]

For an admissible polynomial P, of degree at most 4n+5, define

\[
 L_n(P)=\sum_{b=1}^3w_b\sum_{j=1}^b\sum_{\nu\ge0}
 R_n(\nu+j/b)
 \left[P(\nu+j/b+n/2)+P(-\nu-j/b-n/2)\right]. \tag{1}
\]

The series converges absolutely because deg R_n=−4n−7. Reflection gives
R_n(−t−n)=−R_n(t). Partial fractions therefore remove even zeta values,
and the dilation weights remove ζ(3),ζ(5). Consequently

\[
 L_n(P)\in\mathbb Q+\mathbb Q\zeta(7).
\]

There is no differentiation of point masses in this construction.

## 2. The ordinary signed measure and its nonvanishing mechanism

On X=Z+{0,1/2,1/3,2/3}, define the periodic weights

\[
 \widehat w(x)=\begin{cases}
 37,&x\in\mathbb Z,\\
 -9,&x\in\mathbb Z+1/2,\\
 1,&x\in\mathbb Z+1/3\text{ or }\mathbb Z+2/3.
 \end{cases}
\]

Put R̂_n(x)=R_n(x−n/2), and give the atom at x the weight
sgn(x)ŵ(x)R̂_n(x) when |x|>n/2 and weight zero otherwise. This is an
ordinary, even signed measure μ_n, and (1) is exactly ∫P,dμ_n. It satisfies
the source's growth condition G_D with D=4n+6.

Let Θ={1/2,1/3,2/3}. Its weights, extended by zero off the support, satisfy
A_n(t)a(t+1)=B_n(t)a(t), where

\[
 A_n(z)=(z+1+n/2)^8(z-n/2)\prod_{\theta\in\Theta}(z+\theta-n/2),
\]
\[
 B_n(z)=(z-n/2)^8(z+1+n/2)\prod_{\theta\in\Theta}(z+\theta+n/2).
\]

Both have degree 12. Away from the finite gap this follows from the explicit
ratio R̂_n(z+1)/R̂_n(z); the extra factors vanish at every transition into
or out of the gap. The period-one lattice weights are essential here.

The associated Cauchy–Stieltjes equation has no rational solution. Indeed,
if r were one, the source's argument gives a period-one meromorphic function
g=(f−r)/R̂_n. At sufficiently large points T+1/3 and −T+1/3 its residues
are respectively +1 and −1, because r/R̂_n has only finitely many poles.
Periodicity requires the residues to agree, a contradiction. This verifies
the nonrationality hypothesis independently of any assumption about ζ(7).

**Concrete nonvanishing theorem.** For any even n and integer
1≤M≤2n−8, there is an m with

\[
 M\le m\le M+11,\qquad
 \det_{0\le i,j<m}L_n(t^{i+j})\ne0. \tag{2}
\]

This follows from primary Proposition 5.4 with d=12,D=4n+6:
2(M+11)−2<D is exactly the displayed restriction. Thus one can choose
m/n→c for any fixed 0<c<2, and also approach c=2. The mechanism is a bounded
gap between nonzero Hankel determinants, not positivity of every determinant.

The n=0 version is a useful warning. Its admissible moments are
μ₀=2160ζ(7), μ₂=μ₄=0, and odd moments zero. Hence its second and third Hankel
determinants vanish despite infinite support. The growth hypothesis supplies
only finitely many moments, and the bounded-gap theorem's size conditions
do not apply. Signed measures must not be treated as positive Gram measures.

## 3. An elementary denominator construction, without Bernoulli moments

The factor 108^n has a transparent role. The three degree-n polynomials

\[
 4^n\frac{(t+1/2)_n}{n!},\qquad
 3^{3n/2}\frac{(t+1/3)_n}{n!},\qquad
 3^{3n/2}\frac{(t+2/3)_n}{n!}
\]

are integer-valued at every integer. At primes not dividing the fractional
denominator, count roots in residue classes to supply v_p(n!); at the
exceptional prime the displayed power clears n+v_p(n!). Their product is
G_n(t)=108^n∏_{θ∈Θ}(t+θ)_n/(n!)³, so
R_n=(n!)⁷G_n/(t)_{n+1}⁷.

If a degree-at-most-n polynomial is integer-valued, then its Taylor coefficient
of order j at an integer is integral after multiplication by D_n^j, where
D_n=lcm(1,…,n). Expand it in the integer binomial basis and use Vandermonde's
identity; a degree-j coefficient of binom(u,k) has a denominator dividing
a product of j integers at most n. This property is preserved under products.

For P with this Taylor property, write
R_n(t)[P(t+n/2)+P(−t−n/2)]=∑_{ℓ=1}^7∑_{k=0}^n c_{ℓ,k}/(t+k)^ℓ.
The leading pole factor (n!)⁷/[k!(n−k)!]⁷ is integral. Expanding the remaining
unit factors proves

\[
 D_n^{7-\ell}c_{\ell,k}\in\mathbb Z.
\]

Summing over j/b gives b^ℓ(ζ(ℓ)−H_{bk}^{(ℓ)}). The even coefficients cancel
by reflection, the simple-pole coefficient sums to zero by decay, and ζ(3),ζ(5)
cancel by the chosen weights. The remaining rational tails have denominators
dividing D_{3n}⁷. Thus the affine polynomial representing L_n(P) satisfies

\[
 D_{3n}^7 L_{n,X}(P)\in\mathbb Z[X],\qquad
 L_{n,X}(P)|_{X=\zeta(7)}=L_n(P). \tag{3}
\]

This deliberately easy bound is not claimed optimal. In particular it does
not silently import the primary paper's stronger arithmetic lemma for a
different denominator set.

Use the integer-valued basis

\[
 U_i(t)=\frac{(t)_i}{b_i},\qquad
 b_i=(n!)^{\lfloor i/n\rfloor}(i-n\lfloor i/n\rfloor)!.
\]

It is a product of binomial polynomials of degrees at most n, so (3) applies
to U_iU_j. Therefore

\[
 \frac{D_{3n}^{7m}}{\prod_{i=0}^{m-1}b_i^2}
    \det_{i,j<m}L_{n,X}(t^{i+j})\in\mathbb Z[X]. \tag{4}
\]

The candidate has ordinary-measure nonvanishing and direct factorial/harmonic
arithmetic. Its difficulty moves to the real bound, as the next theorem shows.

## 4. A rigorous no-go for the total-variation estimate

Let H^{abs}_{n,m} be the Hankel determinant of the **positive** measure |μ_n|.
For admissible sizes m≤2n+3 with m/n→c∈(0,2], define the optimistic
absolute-value majorant

\[
 J^{abs}_{n,m}=\frac{D_n^{7m}}{\prod_{i=0}^{m-1}b_i^2}\,H^{abs}_{n,m}.
\]

Even this uses D_n rather than the larger D_{3n} actually proved above.
Nevertheless

\[
 \boxed{\liminf_{n\to\infty}\frac{\log J^{abs}_{n,m}}{n^2}
       \ge c\log54-c^2\log2>0.} \tag{5}
\]

**Proof.** Choose the m positive atoms x_j=n/2+j/2, 1≤j≤m. Their absolute
lattice weights are 9 or 37, hence at least one. Their Vandermonde square is

\[
 \prod_{i<j}(x_j-x_i)^2
   =2^{-m(m-1)}\prod_{i=1}^{m-1}(i!)^2.
\]

The Cauchy–Binet/Andreief expansion of a positive moment determinant contains
this configuration as a positive term. Also b_i≤i!, by the multinomial
coefficient i!/[(n!)^{⌊i/n⌋}(i mod n)!]∈Z. It follows that

\[
 J^{abs}_{n,m}\ge D_n^{7m}2^{-m(m-1)}
                          \prod_{j=1}^m R_n(j/2). \tag{6}
\]

The rational function R_n(t) is decreasing for t>0: the three numerator
reciprocal sums are each bounded by the corresponding denominator sum,
while that denominator occurs seven times. Stirling gives

\[
 \frac1n\log R_n(n)\longrightarrow\log(27/64).
\]

For c≤2 the same lower exponential rate applies to every j/2≤m/2;
if c=2 and m slightly exceeds 2n, replace n by m/2 in this estimate and use
continuity of the limiting rate. Finally
log D_n≥n log2−O(log n): the central binomial coefficient binom(2r,r)
divides D_{2r}, since each prime-power level contributes at most one to its
valuation, and binom(2r,r)≥4^r/(2r+1). Substituting these elementary bounds
in (6) proves (5), because 7log2+log(27/64)=log54.
For c≤2 the exponent is at least c log(54/4)>0. QED.

Thus the standard bound |det μ_n|≤det |μ_n| is quantitatively unusable for
this concrete single-ζ(7) family and normalization: its majorant itself grows
exponentially on the n² scale. This is stronger than merely obtaining an
upper estimate that fails to decay. It still says nothing about cancellation
in the actual signed determinant, or about substantially stronger arithmetic
normalizations. Those are precisely the remaining possible escapes.

## 5. Reproducible checks and next gate

`python codex/computations/check_signed_zeta7_model.py` verifies elimination,
reflection, odd-moment vanishing and the denominator bound for ten basis
entries at n=2. For example the first three even moments are

\[
\begin{aligned}
L_2(1)&=7309440\zeta(7)-107216919/16,\\
L_2(t^2)&=2332800\zeta(7)-20864043/16,\\
L_2(t^4)&=2332800\zeta(7)-12739959/16.
\end{aligned}
\]

The exact theorem (2) resolves nonvanishing for large n without invoking
positivity, while (5) identifies the next obstacle sharply. A useful further
attempt would need an analytic estimate preserving the periodic signed
cancellation—such as a determinant identity before taking absolute values—or
a different arithmetic construction. Merely switching from a positive measure
to this signed one and applying the paper's general nonvanishing theorem does
not complete a ζ(7) argument.
