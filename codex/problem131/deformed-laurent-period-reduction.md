# All deformed reciprocal periods live in four elementary polylogarithmic components

Status: independently reviewed; see
[deformed-laurent-period-independent-review.md](deformed-laurent-period-independent-review.md).
The elementary two-integration-by-parts proof for every q>0 is also
independently reviewed; see section 6 of that review. The Mellin proof is optional.
No new irrationality theorem or historical novelty is asserted. This extends
[laurent-deformation-arithmetic.md](laurent-deformation-arithmetic.md) beyond
its base period and complements the critical-point geometry.

For n≥0 and q>0 define

\[
 I_n(q)=\int_{(0,1)^3}
 \frac{[a(1-a)b(1-b)c(1-c)]^n}{(1-c+qabc)^{n+1}}\,da\,db\,dc
 =\int_{\mathbb R^3}\Lambda_q(e^u,e^v,e^w)^{-n-1}\,du\,dv\,dw.
\]

Both integrals are positive and finite. For 0<q<1 put L=log q and

\[
 F(q)=2\operatorname{Li}_3(q)-L\operatorname{Li}_2(q),\qquad
 G(q)=\operatorname{Li}_2(q)-L\operatorname{Li}_1(q),\quad
 \operatorname{Li}_1(q)=-\log(1-q).
\]

All logarithms in the formulas below are real, and all polylogarithms are
initially defined by their convergent power series.

## 1. An explicit all-n formula

Write Hᵣ=Σ_{h=1}ʳ1/h, H₀=0, and, for 0≤i≤n,

\[
 a_i=\binom ni^2\binom{n+i}{i}^2,\qquad
 b_i=2a_i(2H_i-H_{n+i}-H_{n-i}),\qquad j_i=n+1+i.
\]

Let P_{j,r}(q)=Σ_{h=1}^{j-1}qʰ/hʳ (the empty sum is zero). Then

\[
\boxed{\begin{aligned}
 I_n(q)=\sum_{i=0}^nq^{-j_i}\{&a_i[F(q)-2P_{j_i,3}(q)+L P_{j_i,2}(q)]\\
 &+b_i[G(q)-P_{j_i,2}(q)+L P_{j_i,1}(q)]\}.
\end{aligned}} \tag{1}
\]

In particular every n lies in the fixed space

\[
 I_n(q)\in\mathbb Q(q)F(q)+\mathbb Q(q)G(q)
                 +\mathbb Q(q)\log q+\mathbb Q(q). \tag{2}
\]

This is a containment, not an assertion that these four functions are
linearly independent over Q(q). The coefficient of F is especially concrete:

\[
 \alpha_n(q)=\sum_i a_iq^{-n-1-i}
     =q^{-2n-1}A_n(q),\quad
 A_n(q)=\operatorname{CT}\Lambda_q^n=\sum_i a_iq^{n-i}. \tag{3}
\]

Thus the same constant-term polynomial governs the leading polylogarithmic
component of the reciprocal period, with an explicit reversal and power of q.

## 2. Elementary proof: two Rodrigues integrations by parts

The main argument needs only real integrals, a geometric series, and partial
fractions. The Mellin calculation is retained as an optional alternative in
§6. The following parameter-explicit Beukers transformation was derived by
the root agent while seeking an elementary replacement of that calculation.

Set

\[
 L_n(x)=\frac1{n!}\frac{d^n}{dx^n}[x^n(1-x)^n],\qquad
 K(v)=\frac{-\log v}{1-v}
     =\int_0^\infty\frac{dt}{(1+t)(v+t)},\quad K(1)=1.
\]

For every q>0, the exact identity is

\[
 \boxed{q^n I_n(q)=\int_0^1\!\int_0^1
              L_n(x)L_n(y)K(qxy)\,dx\,dy.} \tag{4a}
\]

Here is the full change-of-variables proof. Insert the positive real integral
for K and integrate by parts n times in x. Rodrigues's factor and the nth
kernel derivative have canceling signs and factorials, yielding

\[
 q^n\int_{(0,1)^2}\!\int_0^\infty
 \frac{x^n(1-x)^n y^n L_n(y)}{(1+t)(qxy+t)^{n+1}}\,dt\,dx\,dy.
\]

Before integrating by parts in y, substitute t=yu. Its Jacobian cancels the
remaining y powers, giving

\[
 q^n\int_{(0,1)^2}\!\int_0^\infty
 \frac{x^n(1-x)^n L_n(y)}{(1+yu)(qx+u)^{n+1}}\,du\,dx\,dy.
\]

Now integrate by parts n times in y. This gives the positive integral

\[
 q^n\int_{(0,1)^2}\!\int_0^\infty
 \frac{[x(1-x)y(1-y)u]^n}
 {(1+yu)^{n+1}(qx+u)^{n+1}}\,du\,dx\,dy. \tag{4b}
\]

Finally u=(1−c)/(yc), with |du/dc|=1/(yc²). Since
1+yu=1/c and qx+u=(1−c+qxyc)/(yc), the u-dependent factors become
[c(1−c)]ⁿ/(1−c+qxyc)ⁿ⁺¹. Formula (4b) is exactly qⁿIₙ(q), proving (4a).

All boundary and convergence points can be checked before these operations.
The derivatives of xⁿ(1−x)ⁿ through order n−1 vanish at both endpoints,
and likewise in y. For fixed positive y,t the x integration by parts is
ordinary integration of smooth functions; after t=yu the y kernel is smooth
on its closed interval for fixed x,u>0. Initially absolute convergence is
bounded by Cₙ²J(q), where Cₙ=max|Lₙ| and J(q)=∫∫K(qxy)<∞. After the first
integration by parts use
qⁿxⁿyⁿ/(qxy+t)ⁿ≤1 and (1−x)ⁿ≤1 to dominate the absolute integrand by
Cₙ/[(1+t)(qxy+t)], whose integral is CₙJ(q). Substitution preserves this
absolute integrability. The second integration by parts is valid pointwise
in x,u; its resulting nonnegative y integral is bounded in absolute value
by the preceding absolute y integral. Integrating in x,u is therefore
legitimate and finite. No informal interchange of a conditionally convergent
integral is needed.

### Polynomial moments explain the squared rational function

Let

\[
 M_n(z)=\int_0^1x^zL_n(x)\,dx
 =(-1)^n\frac{\prod_{j=0}^{n-1}(z-j)}
                  {\prod_{j=1}^{n+1}(z+j)},\quad z>-1. \tag{4c}
\]

For z>n−1, n integrations by parts and an elementary beta integral prove
this formula. Expanding the polynomial Lₙ shows its left side is a rational
function of z, so the identity holds wherever both rational functions are
defined, in particular all real z>−1. Differentiation under the integral
also gives Mₙ′(z)=∫xᶻlog(x)Lₙ(x)dx.

For 0<q<1 expand
K(qxy)=Σ_{k≥0}qᵏxᵏyᵏ[−log q−log x−log y] in (4a). Its nonnegative
summands and the bound |Lₙ(x)Lₙ(y)|≤Cₙ² make termwise integration valid by
absolute domination. The resulting coefficient is
−log(q)Mₙ(k)²−2Mₙ(k)Mₙ′(k). Now

\[
 M_n(k)^2=R_n(k+1),\qquad
 R_n(t)=\frac{\prod_{j=1}^n(t-j)^2}{\prod_{j=0}^n(t+j)^2}.
\]

For k<n the square and its derivative vanish, because these are double
zeros. Shift k=n+ℓ, cancel the factor qⁿ in (4a), and define

\[
 r_n(z)=R_n(z+n+1)=
 \frac{\prod_{h=1}^n(z+h)^2}
      {\prod_{h=n+1}^{2n+1}(z+h)^2}.
\]

We have proved, without complex integration,

\[
 \boxed{I_n(q)=\sum_{k=0}^\infty
 q^k[-r_n(k)\log q-r_n'(k)].} \tag{5}
\]

The series converges absolutely for0<q<1. Its summands need not have a fixed
sign; positivity comes from (4b) or the original triple integral.

The partial fraction decomposition is

\[
 r_n(z)=\sum_{i=0}^n\left(\frac{a_i}{(z+j_i)^2}
                         +\frac{b_i}{z+j_i}\right). \tag{6}
\]

Taking the double-pole limits gives the displayed binomial formula for aᵢ.
The logarithmic derivative of the remaining factors gives
bᵢ/aᵢ=2(2Hᵢ−Hₙ₊ᵢ−Hₙ₋ᵢ). There is no polynomial part since rₙ=O(z⁻²).
Substitute (6) into (5), then use
Σ_{k≥0}qᵏ/(k+j)ʳ=q⁻ʲ[Liᵣ(q)−Pⱼ,ᵣ(q)], proving (1).

### The original Zudilin rational function appears without recurrence guessing

The squared polynomial moment and its shifted version are exactly

\[
 r_n(k)=R_n(k+n+1),\qquad
 R_n(t)=\frac{\prod_{j=1}^n(t-j)^2}{\prod_{j=0}^n(t+j)^2}.
\]

Since Rₙ and Rₙ′ both vanish at the integers1,…,n, equation (5) is

\[
 \boxed{I_n(q)=-q^{-n-1}\sum_{t=1}^{\infty}
      \frac{d}{dt}\big[q^tR_n(t)\big],\quad0<q<1.} \tag{5a}
\]

The derivative acts on the continuous variable t before evaluation at an
integer. At q=1 the identity becomes Iₙ(1)=−Σₜ₌₁∞Rₙ′(t).
Indeed Rₙ(t)=O(t⁻²), Rₙ′(t)=O(t⁻³); dominated convergence treats the latter
sum, and |log q|ΣqᵗRₙ(t)→0 treats the former. Thus this real-variable argument
identifies the positive Beukers period with the Zudilin rational-derivative
series directly, without first proving that both satisfy the same recurrence.
The square gives both the double-pole partial fractions and the double zeros
that make extension of the summation range harmless. The q deformation shows
why a logarithmic correction accompanies the derivative away from1.

This special-function bridge is classical in the Barnes-integral approach;
see [rodrigues-primary-source-comparison.md](rodrigues-primary-source-comparison.md)
for the precise Nesterenko/Zudilin source dictionary. The point here is the
explicit elementary design explanation and its all-n parameter bookkeeping.

## 3. The first deformation, with every extra term visible

For n=1 the pairs (aᵢ,bᵢ) are (1,−4),(4,4). Thus

\[
\boxed{I_1(q)=\frac{(q+4)F(q)+4(1-q)G(q)
                         +q(8\log q-12)}{q^3}.} \tag{7}
\]

A further direct integral check reduces n=1 to a single real integral:

\[
 I_1(q)=\int_0^1 u[(1+u)(-\log u)-2(1-u)]C_1(qu)\,du,
\]
\[
 C_1(v)=\int_0^1\frac{c(1-c)}{(1-c+vc)^2}\,dc
       =\frac{(1+v)\log v-2(v-1)}{(v-1)^3},\qquad C_1(1)=\frac16.
\]

The product-variable Jacobian gives the u factor: integrate a from u to1
with b=u/a. This one-dimensional integral also gives a stable independent
numerical check of (7).

## 4. Why q=1 removes the unwanted periods for every n

Since rₙ(z)=O(z⁻²), the coefficient of z⁻¹ in (6) vanishes:

\[
 \sum_{i=0}^n b_i=0. \tag{8}
\]

As q→1⁻, F→2ζ(3), G→ζ(2), and L→0; in particular
log(q)log(1−q)→0. Therefore (1) becomes

\[
\boxed{I_n(1)=2A_n(1)\zeta(3)
 -\sum_{i=0}^n\left[2a_iH_{n+i}^{(3)}+b_iH_{n+i}^{(2)}\right].} \tag{9}
\]

The ζ(2) term cancels by (8), while the explicit logarithmic terms vanish.
The remaining finite sum is rational. For n=0 this is 2ζ(3); for n=1 it is
10ζ(3)−12, exactly the Apéry/Beukers normalization.

This is stronger than merely observing that the base period has no log at1:
the same pole-residue identity cancels ζ(2) at every n. Away from1, the
simple-pole coefficients are weighted by different q⁻ʲ. In fact, for n≥1,
their weighted sum vanishes at **exactly one positive parameter**, namely1.
This strengthening was independently derived by the integral-review agent.

To prove it, set δᵢ=2Hᵢ−Hₙ₊ᵢ−Hₙ₋ᵢ, so bᵢ=2aᵢδᵢ. One has δ₀<0,
δₙ>0, and

\[
 \delta_{i+1}-\delta_i=\frac2{i+1}-\frac1{n+i+1}+\frac1{n-i}>0.
\]

Thus the nonzero coefficients bᵢ have exactly one sign change. Choose a
real c separating their negative and positive indices. For q>1 the weights
q⁻ⁱ decrease with i, so every nonzero summand satisfies
bᵢq⁻ⁱ<bᵢq⁻ᶜ. Summing and using Σbᵢ=0 gives Σbᵢq⁻ⁱ<0. For0<q<1 the
inequalities reverse. Consequently βₙ(q)=q⁻ⁿ⁻¹Σbᵢq⁻ⁱ has exactly one
positive zero, q=1, and has the asserted opposite signs on the two sides.
Its derivative there is −Σibᵢ=−Σ(i−c)bᵢ<0, so the root is simple.
If a coefficient bᵢ vanishes, choose c at that zero index; the strict
inequalities still hold for every nonzero summand. This elementary
one-crossing balance makes the effect of tilting the weights explicit.
For example β₁(q)=4(1−q)/q³. This is a uniqueness theorem for coefficient
cancellation in (1); it does not assert functional independence of its
polylogarithms or exclude special identities at other q.

The apparent singularities at q=1 in isolated polylogarithmic expressions
are removable for the integral. More strongly, Iₙ is holomorphic in Re q>0.
For q₀>0 and |q−q₀|<q₀, expand the denominator using
abc/(1−c+q₀abc)≤1/q₀. On each smaller disk the binomial series is dominated
by a fixed multiple of the integrable q₀ integrand. This proves local
holomorphy and all derivative interchanges directly. The formulas above use
only the real interval0<q<1 and its limit, so no cut convention for Li at
q>1 is needed.

## 5. A differential recipe, and the q→0 warning

With D=q d/dq,

\[
 DF=G,\qquad DG=-\frac{q\log q}{1-q},\qquad D\log q=1. \tag{10}
\]

Thus the four-component space in (2) is closed under D. The rational Mellin
quotient also gives an exact all-n differential ladder:

\[
 \boxed{(D+2n)^2(D+2n+1)^2 I_n=(D+n)^4 I_{n-1},\quad n\ge1.} \tag{11}
\]

For an elementary proof, any polynomial P satisfies
P(D){qᵏ[−r(k)log q−r′(k)]}
=qᵏ[−P(k)r(k)log q−(P r)′(k)]. Use this in the absolutely convergent
series (5) and the rational identity
rₙ(k)/rₙ₋₁(k)=(k+n)⁴/[(k+2n)²(k+2n+1)²] to prove (11).
Initial data
I₀=F/q and the positive integral select the intended solutions of this
fourth-order ladder. No claim is made that (11) is a minimal differential
or n-recurrence.

At q→0⁺, (5) gives the precise logarithmic divergence

\[
 I_n(q)=-r_n(0)\log q-r_n'(0)+O(q|\log q|),\qquad
 r_n(0)=\frac{(n!)^4}{((2n+1)!)^2}>0. \tag{12}
\]

Although (1) contains negative powers of q, their apparent stronger
singularities cancel. For n=1, I₁(q)=−log(q)/36−1/108+O(q|log q|).
The growing logarithm is consistent with positivity and shows why setting
q=0 inside the original reciprocal integral is not a finite deformation.

The conceptual recipe is now explicit: use a positive Laurent model to
control the small integral, expose its squared rational function through
Rodrigues polynomial moments, and read the arithmetic target from pole
multiplicities and residues.
At q=1 the constant-term polynomial, the residue cancellation, and the
special critical geometry align. Decay alone does not force this alignment.


## 6. Optional advanced proof: the Mellin transform

For real s with 0<s<n+1, apply Tonelli and the beta integral first to q,
then to a,b,c. The q integral contributes B(s,n+1−s), the a and b
integrals each contribute B(n+1−s,n+1), and the c integral contributes
B(n+1−s,s). Consequently

\[
 \mathcal M I_n(s)=\int_0^\infty q^{s-1}I_n(q)\,dq
 =\frac{\Gamma(s)^2\Gamma(n+1-s)^4}{\Gamma(2n+2-s)^2}
 =\frac{\pi^2}{\sin^2\pi s}\,r_n(-s). \tag{4}
\]

Absolute convergence proves the same identity for complex s in that strip.
On Re s=1/2, the right side is integrable in Im s: the rational factor is
O((1+|s|)⁻²) at infinity, while the sine factor decays exponentially.
Mellin inversion (equivalently Fourier inversion after q=eˣ) therefore gives
I_n(q) as the vertical contour integral of q⁻ˢ times (4).

For 0<q<1 close this contour to the left through Re s=−N−1/2.
All poles there are the double poles at s=−k, k≥0. On the new vertical
line the integral is O(q^{N+1/2}(N+1)⁻²), with a constant depending on n;
the horizontal segments vanish as their heights tend to infinity because
of the exponential sine bound. These bounds justify the residue sum.
At s=−k the residue is qᵏ[−rₙ(k)log q−rₙ′(k)], since the squared sine
has no simple-pole term. This recovers (5) independently.


The Mellin method is useful when a new model lacks the elementary Rodrigues
transformation. Here both routes lead to exactly the same rational function;
the real-variable route explains its square and initial double zeros more
directly.

The independent review also supplies a branch-free extension of (1) and (11)
to all q>0. Define

\[
 G(q)=\zeta(2)+\int_1^q\frac{\log u}{u-1}\,du,\qquad
 F(q)=2\zeta(3)+\int_1^q\frac{G(u)}u\,du.
\]

The integrands have removable singularities at1. These are real analytic on
the positive axis and agree with the polylogarithmic definitions on(0,1).
Together with analyticity of Iₙ, the real analytic identity theorem extends
the formulas. No choice of a complex logarithm on a polylogarithmic cut is
required.
