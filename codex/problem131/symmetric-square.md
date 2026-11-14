# Why a cubic recurrence hides a second-order differential equation

Status: all displayed identities below proved algebraically and independently checked by `codex/computations/check_symmetric_square.py`. The generating-function symmetric-square connection is established classical mathematics, explicitly discussed in F. Beukers's retrieved *Consequences of Apéry's work on ζ(3)*, pp. 2–3, https://webspace.science.uu.nl/~beuke106/caen.pdf. Our contribution here is a reproducible derivation, not a novelty claim. The modular and explicit elliptic parametrizations are now proved in the companion [apery-modular-elliptic.md](apery-modular-elliptic.md); they are not required for the differential-algebra argument here.

## 1. Translate the recurrence into an operator

Let A_0=1,A_1=5, and for n≥1 let

    (n+1)³A_{n+1}-(34n³+51n²+27n+5)A_n+n³A_{n-1}=0.

The standard Apéry formula is A_n=Σ_k binom(n,k)²binom(n+k,k)²; the recurrence can be proved by the telescoping certificate in the companion exposition. For the differential-equation derivation only the recurrence and initial data are needed.

Set A(z)=Σ_{n≥0}A_n z^n and θ=z d/dz. The advantage of θ is that θ(z^n)=nz^n, so polynomial coefficients in n turn directly into polynomial differential operators. With operator order as written (z multiplies after the θ operators act), the recurrence gives

    [θ³-z(2θ+1)(17θ²+17θ+5)+z²(θ+1)³]A=0.       (1)

Check coefficients explicitly: at z^m this is

    m³A_m-(2m-1)(17(m-1)²+17(m-1)+5)A_{m-1}
       +(m-1)³A_{m-2}=0,

which is the original relation with n=m-1. At m=1 it is A_1-5A_0=0. At m=0 it is 0=0.

Expanding θ²f=z²f''+zf' and θ³f=z³f'''+3z²f''+zf', and dividing (1) by z, gives

    z²(z²-34z+1)A'''
      +3z(2z²-51z+1)A''
      +(7z²-112z+1)A'
      +(z-5)A=0.                                 (2)

This is a third-order differential equation even though the discrete recurrence has order two. The cubic dependence on n becomes three derivatives; discrete order instead controls the degree in z. Confusing these two different notions of order obscures the connection.

## 2. A general symmetric-square identity

Suppose u solves

    u''+p(z)u'+q(z)u=0.                            (3)

Then Y=u² satisfies

    Y'''+3pY''+(p'+4q+2p²)Y'+(2q'+4pq)Y=0.        (4)

Proof: differentiate Y=u² three times, use u''=-pu'-qu and its derivative, then collect the coefficients of u²,uu',u'². More explicitly,

    Y'=2uu',
    Y''=2u'²-2puu'-2qu².

Differentiating this last expression and substituting (3) yields (4). Polarization shows that uv also solves (4) whenever u,v solve (3): apply the square identity to (u+v)²-u²-v².

If u,v are independent solutions on a simply connected region avoiding singularities, u²,uv,v² are independent. Indeed a linear relation would make the nonconstant function u/v satisfy a fixed nonzero quadratic polynomial where v≠0; a continuous nonconstant function on an open set cannot take values in that finite root set. Thus the third-order solution space is exactly their span. This construction is called the symmetric square because the degree-two symmetric products of a two-dimensional space form a three-dimensional space.

## 3. Discover the second-order equation from the third-order one

Normalize (2) so its leading coefficient is 1, and write the coefficients as a_2,a_1,a_0. Formula (4) forces

    p=a_2/3,
    q=(a_1-p'-2p²)/4.

There is no guessing in these choices. The decisive remaining condition is

    a_0=2q'+4pq.                                  (5)

For the Apéry equation, with d=z²-34z+1, direct rational simplification gives

    p=(2z²-51z+1)/(z d),
    q=(z-10)/(4z d).

The exact identities are

    3p=3(2z²-51z+1)/(z d),
    p'+4q+2p²=(7z²-112z+1)/(z²d),
    2q'+4pq=(z-5)/(z²d).

Each is checked symbolically with indeterminate z in the accompanying script. Multiplying (3) by zd produces the remarkably compact second-order equation

    (z³-34z²+z)u''+(2z²-51z+1)u'+(z-10)u/4=0.     (6)

The last identity (5) is the nonautomatic coincidence. Generic third-order equations do not satisfy it. It encodes structure beyond the mere fact that a sequence satisfies a recurrence.

## 4. The square root of the Apéry generating function

Seek a formal solution u(z)=Σ_{n≥0}v_n z^n with v_0=1. Equation (6) gives v_1=5/2 and

    (n+1)²v_{n+1}-(34n²+17n+5/2)v_n+(n-1/2)²v_{n-1}=0

for n≥1. This determines a unique formal series. Its first terms are

    u(z)=1+(5/2)z+(267/8)z²+(10225/16)z³+(1836275/128)z⁴+...

By the symmetric-square identity, u² solves (2); its constant term is 1. The coefficient relation (1) determines all later coefficients from that constant term, so

    A(z)=u(z)².

This is a formal-power-series identity; no assumption about convergence is needed for this step. Convergence near zero also follows from ordinary regular-singular differential-equation theory or direct bounds on the binomial sum. The script checks the coefficient identity independently through degree 30, supplementing the all-orders uniqueness proof.

## 5. Why the same quadratic irrational keeps appearing

The finite singularities visible in (6) are z=0 and the roots of d:

    ρ=17-12√2=(√2-1)^4,
    ρ^(-1)=17+12√2=(√2+1)^4.

The two numbers are also the roots of λ²-34λ+1=0, the characteristic equation obtained by retaining the leading-degree terms in Apéry's recurrence. Beukers's positive integral has maximum exactly ρ, proved by two square inequalities in the companion note.

These are three precisely checkable manifestations of the same algebraic numbers:

1. characteristic roots of the limiting recurrence;
2. nonzero finite singularities of the differential equation;
3. sharp exponential bound from the integral kernel.

The first two correspond directly under the generating-function translation; reciprocation exchanges the two roots because their product is 1. The integral calculation independently recovers the small root. The identity of the integral with Zudilin's linear form, via the Mellin transform, explains why these rates belong to the same approximation problem. Proving a complete asymptotic expansion, including polynomial factors and nonzero constants, requires additional analysis and is not claimed from singularity locations alone.

## 6. What is simplified, and what is not yet established

Equation (6) lowers the differential order and gives a systematic test for hidden structure: from any proposed third-order equation compute p and q, then test (5). It makes the two-dimensional origin of a three-dimensional solution space explicit, and suggests a search for periods of elliptic curves or modular parametrizations.

However, a symmetric-square factorization alone does not prove ζ(3) irrational. It supplies neither integral coefficients nor denominator bounds for the companion rational approximations, and it does not identify an approximation limit as ζ(3). Those arithmetic and identification steps remain essential. Nor does it automatically yield a ζ(5) construction. A serious advanced simplification should retain these requirements rather than replacing a complete elementary proof with a suggestive differential equation.

The companion [apery-modular-elliptic.md](apery-modular-elliptic.md) now proves the eta-product parametrization by two holomorphic weight-4 identities and a justified finite vanishing check. It also gives an explicit hypergeometric pullback and elliptic quartic period, and explains integrality and least-common-multiple denominators in the modular parameter. The present note remains the prerequisite differential reduction.
