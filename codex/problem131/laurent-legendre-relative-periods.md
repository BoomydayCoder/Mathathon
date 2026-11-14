# The Laurent model gives Legendre P/Q, with a necessary boundary term

Status: exact one-variable reductions and divergence identity; exact algebraic
checks of the Apéry operator's symmetric-square and reciprocal structure.
The full two-variable exact-form certificate was a gap when this note was
written; it is now supplied by the [geometric divergence proof](geometric-divergence-proof.md),
with [completed independent review](geometric-divergence-independent-review.md).
The homogeneous versus relative-period distinction below remains essential.

## 1. The real generating function is not a homogeneous period solution

For

\[
 \Lambda(x,y,z)=\frac{(1+x)(1+y)(1+z)((1+x)(1+y)+xyz)}{xyz},
\]

write A_n=CT Λⁿ and
F_n=∫_{R³}Λ(e^u,e^v,e^w)^{−n−1}du,dv,dw. The previously established
transformations identify A_n with the Apéry integers and F_n with Beukers'
positive integral. Their initial values are

\[
 A_0=1,\quad A_1=5,\qquad F_0=2\zeta(3),\quad F_1=10\zeta(3)-12.
\]

Let θ=t d/dt and

\[
 \mathcal L=\theta^3-t(2\theta+1)(17\theta^2+17\theta+5)
                       +t^2(\theta+1)^3.
\]

Using the established Apéry recurrence for indices n≥1 gives

\[
 \mathcal L\sum_{n\ge0}A_nt^n=0,\qquad
 \boxed{\mathcal L\sum_{n\ge0}F_nt^n=-12t.} \tag{1}
\]

Even without invoking the recurrence at higher indices, the initial values
already prove that the second generating function cannot be annihilated by
the homogeneous operator: its coefficient of t after applying L is −12.
The real integral is over a noncompact chain. Any geometric certificate must
account for its boundary, or formulate the problem in relative cohomology.
The common homogeneous recurrence begins at n=1, not n=0. If a homogeneous
operator for both generating functions is desired, the fourth-order operator
(θ−1)∘L annihilates both; this does not remove the explanatory boundary term.

## 2. Integrating one variable exposes the two Legendre solutions

Put

\[
 p=(1+x)(1+y),\quad q=\frac{p^2}{xy},\quad\delta=q-p.
\]

Then the Laurent polynomial becomes

\[
 \boxed{\Lambda=pz+(p+q)+q/z=(1+z)(p+q/z).} \tag{2}
\]

The elementary constant-term identity for Legendre polynomials gives

\[
 \boxed{\operatorname{CT}_z\Lambda^n
  =\delta^n P_n\!\left(\frac{p+q}{\delta}\right).} \tag{3}
\]

Indeed expanding the last product in (2) gives
∑_{k=0}^n binom(n,k)²p^kq^{n−k}, the standard polynomial expression on the
right. Although its displayed form uses δ in the denominator, the expression
is a polynomial in p,q and extends to δ=0. Taking CT_x CT_y gives directly

\[
 A_n=\sum_{j=0}^n\binom nj^2\binom{n+j}{j}^2.
\]

For positive x,y we have q>p>0. Substitute z=√(q/p)e^v. Then
Λ=p+q+2√(pq)cosh v. The positive real z integral is the classical second
Legendre solution:

\[
 \boxed{\int_{\mathbb R}\Lambda(e^w)^{-n-1}\,dw
  =\frac{2}{\delta^{n+1}}
    Q_n\!\left(\frac{p+q}{\delta}\right).} \tag{4}
\]

Here Q₀(ξ)=½log((ξ+1)/(ξ−1)), ξ>1; no branch ambiguity occurs. The next
section gives an elementary recurrence and initial-value proof of (4), so
one need not assume an integral representation of Q_n.

Thus the compact constant-term construction and positive real integral select
**P_n and Q_n**, the familiar polynomial and logarithmic Legendre solutions.
This is an explicit structural explanation for having two different solutions
before the remaining x,y integrations. It does not yet prove that those final
integrations obey the third-order Apéry operator.

## 3. An exact one-variable divergence certificate, including its boundary

Let a=p+q, b=2√(pq), δ²=a²−b² and D(v)=a+b cosh v. Direct differentiation gives

\[
 b\frac{d}{dv}\left(\frac{\sinh v}{D(v)^{n+1}}\right)
 =-nD^{-n}+(2n+1)aD^{-n-1}-(n+1)\delta^2D^{-n-2}. \tag{5}
\]

For n≥1 the boundary term vanishes at both ends. If
J_n=∫_R D(v)^{−n−1}dv, equation (5) gives

\[
 (n+1)\delta^2J_{n+1}=(2n+1)aJ_n-nJ_{n-1},\qquad n\ge1.
\]

At n=0, however, sinh v/D(v) tends to ±1/b. Its boundary contribution is
2/b and the relation is instead

\[
 \delta^2J_1=aJ_0-2.
\]

An elementary substitution gives J₀=log(q/p)/δ=2Q₀(a/δ)/δ. The last relation
then gives Q₁(ξ)=ξQ₀(ξ)−1, and the recurrence proves (4) for every n. This
small calculation already displays the exact phenomenon responsible for the
inhomogeneous initial data: a boundary that vanishes at positive indices
survives at index zero.

The compact branch satisfies the corresponding polynomial recurrence

\[
 (n+1)C_{n+1}=(2n+1)aC_n-n\delta^2 C_{n-1},
 \quad C_n=\operatorname{CT}_z\Lambda^n.
\]

These identities contain variable coefficients a,δ depending on x,y. They
cannot simply be integrated and replaced by the scalar coefficients 34 and 5.
That remaining reduction requires a separate two-variable certificate.
The reviewed geometric divergence proof now supplies it, with a rational
three-variable lift and controlled real boundary terms. This does not justify
replacing the variable coefficients in these one-variable identities by scalars.

## 4. Compute the final −12 without evaluating ζ(3)

The established Beukers transformation writes

\[
 F_n=\int_0^1\!\int_0^1L_n(x)L_n(y)
          \frac{-\log(xy)}{1-xy}\,dx\,dy,
 \qquad L_1(x)=1-2x.
\]

See [`beukers-explained.md`](beukers-explained.md) for the complete derivation
and original references. Put K(u)=−log u/(1−u). Then

\[
 F_1-5F_0=\iint[-4-2x-2y+4xy]K(xy)\,dx\,dy.
\]

The term 4(xy−1)K(xy)=4log(xy) integrates to −8. For the other terms,
the substitution u=xy gives

\[
 \iint xK(xy)\,dx\,dy
 =\int_0^1(1-u)K(u)\,du=\int_0^1-\log u\,du=1.
\]

The y term is identical, yielding another −4. Thus F₁−5F₀=−12 without
knowing F₀ or evaluating any zeta value. This computes the inhomogeneous
constant directly from the double-integral geometry. It is not by itself a
full divergence certificate at general n; the separate reviewed geometric
proof now constructs that certificate and evaluates its surviving corner.

## 5. Two exact operator checks preceding the full certificate

First, L is a symmetric square. The second-order equation is

\[
 y''+\frac{2t^2-51t+1}{t(t^2-34t+1)}y'
       +\frac{t-10}{4t(t^2-34t+1)}y=0. \tag{6}
\]

If y₁,y₂ solve (6), then y₁²,y₁y₂,y₂² solve Lf=0. This follows by substituting
the standard symmetric-square equation
f'''+3af''+(a'+4b+2a²)f'+(2b'+4ab)f=0 and checking the coefficients.
The symmetric-square viewpoint already appears in Beukers' later exposition;
we independently verified the displayed algebra, not historical novelty.

Second, L has the reciprocal identity

\[
 \mathcal L_t\big[t^{-1}f(1/t)\big]
       =-t\,(\mathcal L_u f)(1/t). \tag{7}
\]

This is immediate from θ_t[t⁻¹f(1/t)]=t⁻¹(−θ_u−1)f(u). It explains why
resolvents 1/(1−tΛ) and 1/(Λ−t), related by parameter inversion and a scalar
factor, can share the same operator. It does not make the noncompact real
chain into a closed cycle or remove its boundary.

`python codex/computations/check_laurent_legendre_periods.py` checks the
constant-term formula, the divergence numerator, and the symmetric-square
identity exactly. The historical gap was an explicit x,y exact-form reduction
deriving L from the geometry, with controlled boundary terms for the real
chain. The reviewed geometric proof now closes it: unknown operator
coefficients are determined by a finite rational system, then a rational lift
and two edge fluxes give \(\mathcal LF=-12t\). Its sizeable certificates do
not yet constitute a short hand derivation. This note's one-variable reduction
remains useful for explaining why the real period has a boundary term; it
does not make that period a homogeneous solution.
