# What the confluent Wilson transform gives toward a Rodrigues formula

Status: bounded constructive investigation following the independent audit.
The exact identities below are proved and checked. They do not yet prove the
all-balanced factorial constant term or an all-n Rodrigues formula. The
original monic polynomials and their conjectures are left unchanged.

## 1. The obvious integer shift along this path degenerates

Let w_a(y) be the unscaled gamma quotient with parameters(a,1−a,0,1).
Gamma recurrence gives, for y>0,

    w_(a+1)(y)/w_a(y)
      = |(a+iy)/(−a+iy)|² =1.

The same holds for the positive-parameter companion(a,1−a,1,1). Equivalently,
the denominator sinh²πy+sin²πa is periodic in a with period1. Thus the
particular paired shift a→a+1,1−a→−a gives no nonconstant Christoffel factor.
It cannot be inserted directly as the usual one-parameter polynomial raising
step. This is a statement about this specific shift, not an impossibility
claim for other Rodrigues constructions. A shift crossing a zero real gamma
parameter also does not license replacing the real-axis integral by an
unjustified formal Wilson contour integral.

## 2. An exact second-antiderivative hierarchy in the parameter

Define the polynomial numerator

    N_m[p](a)=sin²πa · T_m[p](sin²πa),

and E=t d/dt. Its Bernoulli expression proves

    N_m[t p]″(a)
      + N_m[(2E+2m+2)(2E+2m+3)p](a)
      = 2 µ3(t^m p).                                      (1)

The right side is a constant in a. To derive it, for one monomial t^e put
D=2m+2e+2. Then

    N_m[t^e]=(-1)^(m+e+1)[B_D(a)−B_D]/D,

and differentiation of the t^(e+1) expression twice gives
(-1)^(m+e+2)(D+1)B_D(a). The first term containing B_D(a) cancels against
the displayed Euler polynomial. The remaining Bernoulli constant is
(-1)^(m+e)(D+1)B_D=2µ3(t^(m+e)). Linearity proves(1).

For a balanced polynomial with n≥2, the right side vanishes: it is the j=1,
k=1 orthogonality condition. More generally, for0≤k≤n−2,

    N_m[t^(k+1)p]″
     +N_m[t^k(2E+2m+2k+2)(2E+2m+2k+3)p]=0.                (2)

This is a genuine differential constraint coupling the consecutive confluence
conditions. It keeps the Bernoulli endpoint correction rather than dropping
it. It does not close as a fixed scalar differential equation for N_m[p]
alone: the Euler polynomial acts on the unknown polynomial coefficients,
and is not a scalar on a general p. Deriving a closed raising operator from
(2) remains an additional task.

## 3. Both boundary cases of the constant-term conjecture are classical

For m=1, the ordinary polynomial is the monic Wilson polynomial with
parameters(0,1,1,1). In the standard finite hypergeometric formula, at t=0
only the zeroth term remains. Before monic normalization its value is(n!)³,
and its leading coefficient is(−1)^n(n+2)_n. Therefore

    p_(1,n)(0)=(-1)^n (n!)²/binom(2n+1,n).

At t=−1 the series stops after the first term, giving

    p_(1,n)(−1)/p_(1,n)(0)=1+n(n+2)=(n+1)².

For n=1 and arbitrary m, the separately proved cotangent polynomial gives

    p_(m,1)(0)=(-1)^m/(2m+1),
    p_(m,1)(−1)/p_(m,1)(0)=2^(2m).

These prove the two boundary cases of the conjectures. They do not extend
the factorial evaluation to m,n>1.

Nor is the general balanced polynomial just the scalar Wilson polynomial
of degree mn at the endpoint weight. A sharp small counterexample is
m=2,n=1:

    p(t)=t²−2t+1/5,
    µ3(t p)=0,  but  µ3(t² p)=2/15 ≠0.

Its first endpoint scalar weight is tµ3. An ordinary degree2 polynomial for
that weight would also need the second displayed integral to vanish. Thus a
successful confluent construction has to retain the genuinely multiple
conditions; simply specializing the scalar Wilson Rodrigues formula loses
them.

## 4. Primary-source scope and the next precise gap

The inspected primary source Beckermann–Coussement–Van Assche,
*Multiple Wilson and Jacobi–Piñeiro polynomials*, arXiv:math/0303093, §3.2,
Theorem3.3 varies one Wilson parameter with three fixed. The current path
varies the complementary pair(a,1−a), so that theorem is not a direct
substitution. Its §4.3.4 formal multiple continuous-Hahn limit likewise
varies one parameter in

    Γ(a+t)Γ(b_j−t)Γ(c+t)Γ(d−t).

The source explicitly distinguishes this formal family from a positive
multiple family with a varying conjugate pair. This confirms that its stated
Rodrigues/transform formulas do not automatically cover our paired path;
it is not evidence that no broader applicable theorem exists.

A possible continuation is a proved quadratic or Fourier–Jacobi transform
that converts the complementary-pair confluence into a known multi-parameter
Rodrigues system, with all parameter derivatives and endpoint factors retained.
Alternatively, one could derive a closed operator from(2), or evaluate the
confluent moment-determinant ratio for p(0). No such transform, closure, or
determinant evaluation was proved in this bounded step.

`check_confluent_wilson_raising.py` checks(1) on twelve arbitrary polynomial
cases, the exact scalar-weight counterexample, and the m=1 endpoint evaluations.
All checks pass. No additional finite fitting was used as a substitute for the
missing all-m,n derivation.

## 5. Targeted quadratic/continuous-Hahn follow-up

The usual gamma duplication formula makes the precise issue visible. For the
positive-parameter unscaled companion quotient w⁺ₐ,

    w⁺ₐ(y)=4π y²
       |Γ(a+iy)Γ(1−a+iy)Γ(1+iy)/Γ(1/2+iy)|².            (3)

Indeed |Γ(2iy)|²=|Γ(iy)Γ(1/2+iy)|²/(4π), and one Γ(1+iy)
cancels Γ(iy) up to iy. At the isolated point a=1/2 this reduces to
4πy²|Γ(1/2+iy)Γ(1+iy)|², an even continuous-Hahn weight times a polynomial.
Along the actual varying paired path, the residual gamma ratio in(3) remains.
This calculation does not supply a common continuous-Hahn intertwiner for
all parameters, and evaluating the isolated half point does not prove the
confluence at a=0. It is not a theorem forbidding a subtler quadratic or
Fourier–Jacobi transformation.

For a targeted later primary comparison I retrieved Jean Paul Nuwacu and
Walter Van Assche, *Multiple Askey–Wilson polynomials and related basic
hypergeometric multiple orthogonal polynomials*,
[arXiv:1904.01252](https://arxiv.org/abs/1904.01252). In §3.5, Theorem3.5,
the weights vary a_j while b,c,d are fixed. Its transform from multiple
little-q-Jacobi polynomials therefore does not directly provide the current
complementary-pair Rodrigues formula either. Local PDF/text files are
`codex/sources/nuwacu-vanassche-multiple-askey-wilson-190401252.*`.

This bounded source comparison and the exact paired-shift/duplication
calculations did not yield an all-m,n Rodrigues formula or a proof of the
factorial p(0) conjecture. They identify concrete hypotheses that a successful
transform must bridge. No historical novelty or general nonexistence conclusion
is drawn from the unsuccessful search.
