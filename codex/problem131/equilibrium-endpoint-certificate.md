# A rigorous endpoint certificate for the signed balayage system

Status: a complete existence-and-uniqueness proof for the stated endpoint equations inside a rational box, backed by integer-only interval quadrature. This validates the endpoint system supplied by the parent research investigation; it does not by itself prove that the associated measure is the equilibrium measure or establish its energy inequality. Those variational and energy claims require their separate arguments.

## 1. The exact system and certified box

Let c=3/40 and λ=37/40, and define, for u>=0,

    q(u)=-3 on [0,c],
          3 on (c,1],
          2 on (1,infinity).

For 0<a<b define

    I(a,b)=integral_0^infinity q(u)/sqrt((u²+a)(u²+b)) du,
    F(a,b)=λ+(1/2)integral_0^infinity
                  q(u)[sqrt((u²+a)/(u²+b))-1]du.

The following rational rectangle contains exactly one solution of I=F=0:

    0.00007119508 < a < 0.00007119517,
    0.75394957 < b < 0.75394967.                         (1)

All decimal endpoints here denote exact terminating rational numbers. The widths are9×10^-11 and10^-7, respectively.

The interval computation certifies the signs on the *entire* four faces, not merely at their corners:

| face | rigorous enclosure |
|---|---|
| a=a_min, all b in the box | I∈[-0.000001673334625,-0.000000415767693] |
| a=a_max, all b in the box | I∈[0.000000483172361,0.000001740739136] |
| b=b_min, all a in the box | F∈[0.000000032042900,0.000000032513736] |
| b=b_max, all a in the box | F∈[-0.000000034809551,-0.000000034338715] |

The displayed enclosures are outward-rounded from the exact dyadic intervals. Their separation from zero proves the signs.

## 2. One sign change gives monotonicity on the zero set

Write

    h_(a,b)(u)=1/sqrt((u²+a)(u²+b)),
    J_a=integral q(u)h_(a,b)(u)/(u²+a) du,
    J_b=integral q(u)h_(a,b)(u)/(u²+b) du.

At a point where I(a,b)=0, both J_a and J_b are strictly negative. For example, subtract the constant1/(c²+a) inside J_a. Its product with I is zero, so

    J_a=integral q(u)h_(a,b)(u)
                   [1/(u²+a)-1/(c²+a)]du<0.

Before c, q is negative and the bracket positive; after c, q is positive and the bracket negative. Both sets have positive measure. The same argument applies to J_b. The later change from3 to2 does not matter, because it preserves the positive sign.

Differentiation under the integral is justified uniformly on the compact box: a,b remain positive and the tails of the derivatives are integrable. It gives

    I_a=-J_a/2>0,       I_b=-J_b/2>0                    (on I=0),
    F_a=I/4,
    F_b=-(1/4)integral q(u)(u²+a)h_(a,b)(u)/(u²+b)du.

On the zero set, the last formula simplifies to

    F_b=(b-a)J_b/4<0.                                 (on I=0)    (2)

These derivative statements are asserted on I=0; they are not presumed to hold everywhere off that curve.

## 3. Existence and uniqueness without a two-dimensional black-box solver

For each fixed b in the certified interval, the signs on the two a-faces and the intermediate value theorem give at least one a with I(a,b)=0.

That a is unique. Indeed if a_2>a_1 and I(a_1,b)=0, then

    h_(a_2,b)(u)/h_(a_1,b)(u)
       =sqrt((u²+a_1)/(u²+a_2))

is strictly increasing in u. Subtract its value at u=c when integrating against q(u)h_(a_1,b)(u). On each side of c the product is positive, so I(a_2,b)>0. Hence two a-roots at the same b are impossible.

Because I_a>0 at each root, the implicit-function theorem supplies a continuously differentiable root curve a=a(b) throughout the interval. Its derivative is -I_b/I_a<0. Along this curve,

    d/db F(a(b),b)=F_a a'(b)+F_b=F_b<0,

using F_a=0 and (2). Thus F decreases strictly along the I-zero curve.

The entire-face bounds for F imply it is positive at b_min and negative at b_max on that curve. The intermediate value theorem gives exactly one joint root. It lies strictly inside the box, proving (1). This is uniqueness within the certified rectangle; no assertion of global uniqueness outside the rectangle is needed.

## 4. Make the infinite integrals finite and numerically stable

The script splits at u=c and u=1. For u>=1 set u=1/v. The two tail integrands become

    h_tail(v)=1/sqrt((1+av²)(1+bv²)),
    g_tail(v)=(a-b)/{
          sqrt(1+bv²)[sqrt(1+av²)+sqrt(1+bv²)]},

integrated over0<=v<=1 with weight2. The second formula follows by rationalizing the numerator before dividing by v². It is continuous at v=0; no cancellation of two nearly equal square roots is required.

On the finite u intervals, the second integrand is similarly evaluated as

    g(u)=(a-b)/{
          sqrt(u²+b)[sqrt(u²+a)+sqrt(u²+b)]}.

The first integrand is h(u)=1/[sqrt(u²+a)sqrt(u²+b)]. All radicands are strictly positive on the parameter box. The rationalized formulas are exact identities.

## 5. What makes the numerical certificate rigorous

`codex/computations/certify_equilibrium_endpoints.py` implements its own fixed-point interval arithmetic using integers at scale2^120. No floating-point arithmetic enters a sign assertion.

- Every rational input is enclosed using integer floor and ceiling division.
- Addition and subtraction are exact at the fixed scale.
- Products and quotients are rounded outwards using integer division; division rejects intervals containing zero.
- Square-root endpoints are obtained by integer `isqrt`, with the upper endpoint incremented unless the scaled radicand is a perfect square.
- The final decimal printing is also outward-rounded.

To certify integration error, the script propagates interval Taylor jets through degree4. The zeroth jet variable ranges over the whole subinterval, while its first coefficient is1. Formal multiplication, reciprocal and square-root coefficient recurrences then enclose f^(4)(u)/4! simultaneously for every u in that subinterval and every allowed parameter value on the selected box face.

For one subinterval of width h, Simpson's rule has error bounded by

    h^5 sup|f^(4)|/2880
       =h^5 sup|jet_4|/120.

The script encloses the Simpson evaluation itself and adds this rigorous symmetric error interval. It uses512 equal subintervals on[0,3/40],256 on[3/40,1], and128 on the transformed tail[0,1]. The finer mesh near zero resolves the small a-scale.

For each I-face, the accumulated Simpson error is less than5.857×10^-8. For each F-face it is less than1.160×10^-11. Parameter uncertainty and rounding widths are included separately in the final interval, not omitted from the error budget.

This calculation is reproducible with the Python standard library alone; it does not rely on SciPy stopping criteria or a floating-point root finder. Approximate roots were used only to choose a candidate rational box. The proof depends solely on the four validated signs and the analytic monotonicity argument.

## 6. What this certificate does and does not settle

The endpoint equations have a uniquely specified solution inside the small rectangle. This removes an existence gap that would remain if one only reported numerical roots.

To infer an equilibrium theorem, one must still verify nonnegativity and mass of the proposed measure, constancy of the appropriate potential on its support, and the off-support variational inequality. To infer a claimed energy or irrationality inequality, its functional must be evaluated with a separate rigorous error bound. Those are distinct mathematical steps, not consequences of a root finder or this endpoint certificate alone.

Run `python codex/computations/certify_equilibrium_endpoints.py` to reproduce all four whole-face sign assertions.
