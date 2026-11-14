# Beukers's integral proof and its exact bridge to Zudilin

Status: complete derivation below, using the independently audited elementary lcm bound in the companion Zudilin exposition. Algebraic identities were independently derived in this session, and the Mellin bridge was independently found by the connections agent. This is an exposition/simplification of classical mathematics, not a claim of a new irrationality proof. Original source: F. Beukers, *A note on the irrationality of ζ(2) and ζ(3)*, Bull. London Math. Soc. 11 (1979), 268–272, DOI https://doi.org/10.1112/blms/11.3.268. Publisher metadata verified; the English publisher PDF was not retrieved, but a Russian translation linked from Zudilin’s bibliography was retrieved at https://www.math.ru.nl/~wzudilin/zw/Beukers79.pdf and saved as `codex/sources/beukers-1979-russian.pdf`. Pages 3–4 were visually inspected because the PDF text extraction has broken font encoding. The integral transformation agrees with the derivation below; the translation displays a stray (-1)^n at one intermediate line, which the direct calculation here avoids. Beukers's own later exposition was retrieved from https://webspace.science.uu.nl/~beuke106/caen.pdf and saved as `codex/sources/beukers-caen.pdf` (its title page dates the file May 9, 2003). Zudilin's original 2002 paper was read locally.

## The destination explains the construction

Let D_n=lcm(1,...,n), D_0=1. We seek a positive linear form I_n in 1 and ζ(3), whose coefficients become integers after multiplication by D_n³, and which decreases faster than D_n^(-3). Under ζ(3)=p/q, qD_n³I_n would then be an integer strictly between zero and one.

Define

    L_n(x) = (1/n!) (d/dx)^n [x^n(1-x)^n]
           = Σ_{k=0}^n (-1)^k binom(n,k) binom(n+k,k) x^k.

These are shifted Legendre polynomials, with L_n(0)=1. The name is optional: their two needed features follow directly from this formula. Their coefficients are integers, and integration by parts transfers their derivatives to a kernel without boundary terms. The endpoint zeros of x^n(1-x)^n have multiplicity n. Those zeros are precisely why n integrations by parts are legal without correction terms.

The integral we first choose is

    I_n = ∫_0^1∫_0^1 L_n(x)L_n(y) [-log(xy)]/(1-xy) dx dy.

Why this kernel? For n=0 its geometric-series expansion gives 2ζ(3). Multiplying by integer polynomials preserves controlled arithmetic. The Legendre choice will also permit transformation to an integral with a positive, exponentially small integrand.

## Arithmetic without mysterious cancellation

Put H_r^(s)=Σ_{j=1}^r 1/j^s, H_0^(s)=0, and let

    J_rs = ∫∫ x^r y^s [-log(xy)]/(1-xy) dx dy,   0≤r,s≤n.

Expanding (1-xy)^(-1) gives

    J_rs = Σ_{k=0}∞ [1/((k+r+1)²(k+s+1))
                         +1/((k+r+1)(k+s+1)²)].

All summands are nonnegative, so this interchange follows from monotone convergence. If r=s,

    J_rr = 2ζ(3)-2H_r^(3).

If r≠s, the summand equals

    [1/(k+r+1)² - 1/(k+s+1)²]/(s-r),

so the series telescopes in shifted indices:

    J_rs = [H_s^(2)-H_r^(2)]/(s-r).

Since |s-r|≤n, D_n³ J_rs is an integer for r≠s; and D_n³(J_rr-2ζ(3)) is an integer. Write c_k=(-1)^k binom(n,k)binom(n+k,k). Expanding L_n(x)L_n(y) now yields

    I_n = u_n ζ(3)-v_n,
    u_n=2Σ_{k=0}^n c_k² ∈ Z,    D_n³v_n∈Z.

The diagonal terms carry ζ(3); the off-diagonal terms are rational. This is the integral version of separating double poles from cross terms in a squared partial-fraction expansion.

## Turning an oscillating polynomial integral into a positive one

First use the elementary identity, for 0<a<1,

    ∫_0^1 dz/[1-(1-a)z] = -log(a)/(1-a).

Thus

    I_n = ∫_[0,1]^3 L_n(x)L_n(y)/Q dx dy dz,
    Q=1-z+xyz.

For fixed y,z in the open interval, integrate by parts n times in x. The derivative of Q^(-1) is

    ∂_x^n Q^(-1)=(-1)^n n! (yz)^n Q^(-n-1).

The integration-by-parts sign cancels (-1)^n, giving

    I_n=∫ x^n(1-x)^n L_n(y)(yz)^n/Q^(n+1) dx dy dz.       (A)

The remaining L_n(y) still changes sign. A direct second integration by parts would differentiate y^n/Q^(n+1), making an unpleasant sum. The following change of variable is designed to remove that y^n first:

    w=(1-z)/(1-z+xyz).

For fixed x,y>0 this reverses [0,1]. Set E=1-w+xyw. Then

    z=(1-w)/E,    Q=xy/E,    |dz/dw|=xy/E².

Consequently

    (yz)^n Q^(-n-1)|dz/dw| = (1-w)^n/(x^n E).

The x^n in (A) cancels, and we obtain

    I_n=∫ (1-x)^n(1-w)^n L_n(y)/E dx dy dw.

Now n integrations by parts in y have a simple kernel again. Since

    ∂_y^n E^(-1)=(-1)^n n!(xw)^n E^(-n-1),

we obtain the central identity

    I_n = ∫_[0,1]^3
          [x(1-x)y(1-y)w(1-w)]^n / [1-w+xyw]^(n+1)
          dx dy dw.                                      (B)

The interior integrand is positive, so I_n>0 for every n. The case n=0 follows directly; integrations by parts are then vacuous.

**Analytic justification.** The original kernel has finite integral 2ζ(3), and each L_n is bounded on [0,1]. Thus the signed original integral is absolutely integrable. The first integration-by-parts formula holds pointwise for y,z∈(0,1). Its positive x-kernel, integrated in x, equals ∫L_n(x)/Q dx and hence is bounded by ||L_n||∞∫dx/Q; multiplying by |L_n(y)| therefore supplies an integrable bound. Change of variables is legitimate on the open cube; its boundary has measure zero. The same argument after the transformation legitimizes the second integration by parts. Alternatively perform the pointwise integrations first and use these bounds to apply Fubini. This avoids silently assuming a singular corner contributes no boundary term.

## The sharp exponential bound using squares, not multivariable calculus

Write

    h(x,y,w)=x(1-x)y(1-y)w(1-w)/(1-w+xyw).

We prove its exact maximum is ρ=(√2-1)^4. Put t=√(xy), so 0≤t≤1. From x+y≥2√(xy),

    (1-x)(1-y)≤(1-t)².

The square [1-(1+t)w]²≥0 is exactly the inequality

    w(1-w)/(1-(1-t²)w)≤1/(1+t)².

Therefore

    h≤[t(1-t)/(1+t)]².

Put a=√2-1, which satisfies a²+2a=1. Then

    a²(1+t)-t(1-t)=(t-a)²≥0,

so h≤a⁴. Equality occurs at

    x=y=a,   w=1/(1+a)=1/√2.

This provides both the sharp bound and its motivation: first collapse x,y to their product, then optimize the fractional-linear w factor by completing a square, then complete another square in t.

From (B),

    0<I_n≤ρ^n ∫_[0,1]^3 dx dy dw/(1-w+xyw)
          =2ζ(3)ρ^n.

For n>0 the upper inequality is strict, since equality of h occurs only on a measure-zero set.

## The contradiction and elementary denominator growth

An adequate lcm estimate is D_n≤C·c^n for any fixed c<ρ^(-1/3)≈3.239 (with fixed C). One can use the elementary classical estimate D_n<3^n for n≥1, which Zudilin states. More usefully for a self-contained project, Section 10 of `zudilin-explained.md` proves log D_n≤c n+O(log²(n+1)), c=(6/5)[(log2)/2+(log3)/3+(log5)/5-(log30)/30]≈1.10555043. I independently audited its floor-weight identity, period-30 nonnegative table, factorial estimate, and geometric iteration. Since 3c+log ρ≈-0.208843<0, this slightly weaker bound suffices. Thus the combined notes give an elementary proof without an unproved lcm lemma.

Alternatively the prime number theorem gives log D_n=ψ(n)=n+o(n), which is more than sufficient but adds an advanced prerequisite. This alternative is a complete route using a standard stated theorem.

Assume ζ(3)=p/q. Then qD_n³I_n is a positive integer. With D_n<3^n,

    0<qD_n³I_n<2qζ(3)[27(√2-1)^4]^n →0,

since 27(√2-1)^4≈0.794805<1. Eventually this positive integer is below one, a contradiction. With PNT, the exponential rate is instead exp(3)ρ<1, so the same argument works.

## Exact identity with Zudilin: the Mellin transform reveals the same object

This section proves the equality directly, without invoking irrationality, comparing transcendental coefficients, or guessing a recurrence.

For t>0, polynomial integration gives a rational function

    M_n(t)=∫_0^1 x^(t-1)L_n(x)dx=Σ_k c_k/(t+k).

For t>n, n integrations by parts and the beta integral give

    M_n(t)=(-1)^n [(t-1)...(t-n)]/[t(t+1)...(t+n)].

Both sides are rational functions in t, so the identity, first proved on t>n, holds identically wherever defined, in particular at every positive t. Squaring gives exactly Zudilin's rational function:

    M_n(t)²=R_n(t).

Differentiating under the integral gives

    R_n'(t)=∫∫ (xy)^(t-1)log(xy)L_n(x)L_n(y) dx dy.

Summing over integer t≥1 is legitimate by domination with
||L_n||∞²[-log(xy)]/(1-xy), whose integral is finite. Hence

    -Σ_{t=1}∞R_n'(t)=I_n.

This also explains all the matching binomial coefficients: c_k are simultaneously polynomial coefficients and residues of the Mellin transform. Their squares are simultaneously diagonal terms of the double integral and double-pole coefficients of R_n.

## What this actually simplifies

Combining Zudilin's Lemma 1 with the Mellin bridge and (B) replaces his analytic half: Ball's positive series, its elaborate factorial bound, both recurrence certificates needed to identify it with F_n, and the initial-value comparison are unnecessary for the irrationality proof. The replacement uses two integrations by parts, one explicit fractional-linear substitution, and two completed squares.

This is a concrete reduction in calculation, but it recovers the classical Beukers approach. It is not a new proof of Apéry's theorem. It does not by itself explain the numerical coefficients 34,51,27,5 in the recurrence. Rather, it establishes that the recurrence is optional for irrationality, while the Mellin transform explains why the integral and rational-function approaches share it once it has been proved on either side.

## Advanced connection, with limitations

In his retrieved later exposition Beukers gives the generating-function differential equation for A(t)=Σ a_n t^n, a_n=Σ_k binom(n,k)²binom(n+k,k)², and its representation as the symmetric square of a second-order equation. The singularities include (1±√2)^4. The small singularity is exactly ρ above, and ρ^(-1)=17+12√2 is the large characteristic root of the limiting Apéry recurrence. This numerical coincidence has a precise verified algebraic component: (1+√2)^4=17+12√2 and its reciprocal is (√2-1)^4. A full geometric derivation of the recurrence from periods/modular parametrization remains outside this note; the source should be read critically, as its displayed sign for c_n-a_nζ(3) is inconsistent with n=0 and our positive integral convention.

## Reproducible checks and next steps

`codex/computations/check_beukers.py` verifies exact monomial arithmetic, the integer denominator assertions, the coefficient match with the rational-function formula, both square identities, and the Möbius change of variables. These finite checks supplement the proofs, not replace them.

Next: incorporate the companion’s verified elementary lcm estimate into a unified presentation; develop the recurrence by finite linear algebra from the Mellin rational function; explain the symmetric-square equation from a period model. A proposed generalization must preserve both the positive integral bound and a single desired zeta coefficient—adding dimensions alone does not establish irrationality of ζ(5).
