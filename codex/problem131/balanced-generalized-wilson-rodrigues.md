# A generalized Wilson Rodrigues formula for every balanced polynomial

Status: all-index construction and proof accepted in the
[independent root review](balanced-rodrigues-independent-review.md). The exact computations corroborate the formulas; the contour argument
below is the proof of orthogonality. No historical novelty is asserted.

This resolves the previously open explicit-construction and constant-term
questions in `hyperbolic-balanced-bernoulli-construction.md`. The relevant
connection is a Wilson divided difference with **2m+2 gamma parameters**, rather
than the classical four-parameter multiple-Wilson family. Increasing the number
of gamma factors makes one raising step increase the degree by m. Periodic
multipliers then turn the same adjoint identity into m simultaneous families
of orthogonality conditions.

## 1. Polynomial raising operators and the formula

Fix integers m,n>=1, and write N=2m+1 and D=mn+m+1. For an even polynomial f
and an integer s>=0 define

    A_s(x)=(x+s/2)(x+s/2+1)^N,
    B_s f(x)=[A_s(x)f(x+1/2)-A_s(-x)f(x-1/2)]/(2x).       (1)

The numerator is odd, so B_s f is an even polynomial. Define the composition
with the rightmost operator acting first:

    F_(m,n)=B_0 B_1 ... B_(n-1) 1.                         (2)

Then the monic balanced polynomial is

    p_(m,n)(-x^2)=(-1)^(mn) F_(m,n)(x)/(D)_n,              (3)

where (D)_n=D(D+1)...(D+n-1). Equations (1)-(3) are an explicit finite
Rodrigues-type construction, using rational arithmetic alone.

If f has degree 2r and leading coefficient c, expansion of its two shifts
shows that B_s f has degree 2(r+m) and leading coefficient

    c [r+(m+1)s+2m+1].                                   (4)

At stage s in (2), r=m(n-1-s), so the bracket is D+s. Consequently F has
degree 2mn and leading coefficient (D)_n, proving the monic normalization.
The first operator on an arbitrary even f simplifies to

    B_0 f(x)=[(x+1)^N f(x+1/2)+(1-x)^N f(x-1/2)]/2.       (5)

In particular the n=1 case is the previously known cotangent polynomial.

## 2. The meromorphic weights and their elementary Pearson identity

For s>=0 put

    w_s(x)=Gamma(s/2+x)Gamma(s/2-x)
           [Gamma(1+s/2+x)Gamma(1+s/2-x)]^N
           /[Gamma(2x)Gamma(-2x)].                        (6)

These are even meromorphic functions. On the imaginary axis x=iy, w_s is
positive for y>0. Gamma reflection and recurrence give

    w_0(iy)=4 pi^(2m+1) y^(2m+1)
            coth(pi y)csch(pi y)^(2m).                   (7)

For clarity, |Gamma(iy)|^2=pi/[y sinh(pi y)],
|Gamma(1+iy)|^2=pi y/sinh(pi y), and
1/|Gamma(2iy)|^2=2y sinh(2pi y)/pi; multiplying gives (7).

Gamma recurrence also gives the two exact ratios

    w_(s+1)(x+1/2)/w_s(x)=-(2x+1)A_s(x)/(2x),
    w_(s+1)(x-1/2)/w_s(x)=-(2x-1)A_s(-x)/(2x).            (8)

These are meromorphic identities, so removable special points cause no
restriction. If

    E v(x)=v(x+1/2)/(2x+1)-v(x-1/2)/(2x-1),

then (8) is precisely

    E(w_(s+1) f)=-w_s B_s f.                             (9)

## 3. A justified contour adjoint identity

Define the even-polynomial divided difference

    Dcal g(x)=[g(x+1/2)-g(x-1/2)]/(2x).                  (10)

Let g(x) be an even polynomial multiplied by an even trigonometric polynomial
of period one whose exponential growth is at most exp(2r pi |Im x|), where
0<=r<=m-1. Let f be any even polynomial. Then

    integral_(iR) g(x) w_s(x) B_s f(x) dx
      = integral_(iR) Dcal g(x) w_(s+1)(x) f(x) dx.        (11)

The orientation is upward on both vertical lines. Here and at subsequent
steps the trigonometric factor may be translated by a half-integer; it remains
even, period one, and has the same growth.

To justify (11), set U(x)=w_(s+1)(x)f(x)/(2x). It is holomorphic throughout
a neighborhood of the closed strip |Re x|<=1/2. At x=0 the reciprocal gamma
product in (6) has a zero of order two, while all numerator gamma factors
have positive arguments, so w_(s+1)=O(x^2) and U has a removable zero. For
s>=1 the numerator gamma factors have no poles in this strip. For s=0 the
only possible poles are the simple poles of Gamma(1/2-x) at x=1/2 and
Gamma(1/2+x) at x=-1/2; each is canceled by the simple zero of
1/[Gamma(2x)Gamma(-2x)]. The other numerator factors are regular there.
There are no other candidate singularities in the strip.

Uniform Stirling bounds in a fixed vertical strip give

    |w_s(c+iy)| <= C (1+|y|)^M exp(-2m pi |y|)

for |y| sufficiently large, uniformly in |c|<=1/2; the particular finite
power M may depend on m,s. Multiplication by g, its fixed shifts, f and 1/x
still leaves a polynomial times exp(-2(m-r) pi |y|). Since m-r>=1, all
vertical integrals converge absolutely and horizontal connector integrals
vanish as their heights tend to infinity.

Use (9), change variables by +/-1/2 in its two terms, and deform the resulting
vertical lines back to iR. No pole is crossed by the holomorphy just proved.
The left side of (11) becomes

    - integral_(iR) [g(x-1/2)-g(x+1/2)] U(x) dx,

which is the right side. This checks the sign as well as convergence. The
contour argument uses only gamma recurrence, gamma's pole locations, and its
standard fixed-strip decay; it invokes no general Wilson orthogonality theorem.

## 4. Periodic multipliers give every balanced condition

For j=1,...,m set r=m-j and

    h_j(x)=sin(pi x)^(2r).

This is even and has period one. For every function q for which the displayed
expressions are defined,

    Dcal(h_j q)(x)=h_j(x+1/2)Dcal q(x),

because h_j(x+1/2)=h_j(x-1/2). Iteration therefore gives

    Dcal^n[h_j(x)x^(2k)]
       =h_j(x+n/2)Dcal^n[x^(2k)]=0,  0<=k<n.            (12)

The polynomial divided difference lowers the degree in x^2 by one. The shifted
h_j factors are even for these half-integer shifts (they alternate between an
even power of sine and an even power of cosine), so every application of (11)
is justified by section 3.

Apply (11) n times to the composition (2). By (12),

    integral_(iR) h_j(x)x^(2k) w_0(x) F_(m,n)(x) dx=0
       for 1<=j<=m and 0<=k<n.                           (13)

On x=iy, h_j(iy)=(-1)^r sinh(pi y)^(2r), x^(2k)=(-1)^k y^(2k),
and (7) shows that (13) is a nonzero constant times

    integral_0^infinity y^(2k) p_(m,n)(y^2)
        y^(2m+1)coth(pi y)csch(pi y)^(2j) dy=0.           (14)

The integrand is even on the full imaginary axis, so restricting to y>0
introduces only a factor two. All endpoint integrals converge: at y=0 the
weight in (14) is O(y^(2m-2j)), with nonnegative exponent. This proves all mn
balanced conditions. The established balanced normality theorem identifies
(3) uniquely with the previously constructed p_(m,n). Its rationality is
also immediate directly from (1)-(3), and the previously proved simple
positive roots therefore apply.

## 5. The factorial constant and the first negative-square evaluation

Equation (5) gives B_0 f(0)=f(1/2). For s>=1,

    B_s f(s/2)=(s+1)^N f((s+1)/2),                       (15)

because A_s(-s/2)=0. Chaining (15) through (2) yields

    F_(m,n)(0)=(n!)^N.

Together with (3), this proves the formerly conjectured formula

    p_(m,n)(0)=(-1)^(mn)(n!)^(2m+1)/(mn+m+1)_n
             =(-1)^(mn)(n!)^(2m)/binom((m+1)n+m,n).      (16)

A second collapsing chain begins at x=1. For every s>=0,

    B_s f((s+2)/2)=(s+1)(s+2)^(N-1) f((s+3)/2),         (17)

since A_s(-(s+2)/2)=(-1)*0^N=0. Thus

    F_(m,n)(1)=n!((n+1)!)^(N-1),
    p_(m,n)(-1)/p_(m,n)(0)=(n+1)^(2m).                  (18)

Both identities hold for all m,n>=1. They arise from zeros of the raising
operator coefficient at successive half-integers, not from finite fitting.

## 6. What this establishes, and what it does not

The construction proves an explicit all-index Rodrigues composition, the
factorial constant term, and the negative-square ratio (18). It explains why
the ordinary n-fold difference trial missed the family: the correct adjoint
uses a quadratic-lattice divided difference and the two polynomial factors
in A_s, rather than an ordinary forward difference.

It does not by itself establish useful asymptotic denominator growth for an
odd-zeta irrationality argument, nor optimal approximation errors. The
composition is a concrete starting point for arithmetic analysis. Historical
identification of this generalized gamma Rodrigues family remains to be
checked against primary literature; no novelty claim is made here.

## Reproduction

`python codex/computations/check_balanced_generalized_wilson.py` checks 19
balanced systems against their original rational Bernoulli moments, 20
arbitrary-polynomial leading/endpoint cases, 24 gamma Pearson ratios, and
five nonzero contour-adjoint quadratures. All checks pass. The degree-15
case (m,n)=(3,5) and degree-12 case (2,6) are included.
