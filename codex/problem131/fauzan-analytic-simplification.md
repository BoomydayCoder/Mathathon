# Replacing the sixteen-arcsine certificate by an analytic equilibrium construction

Status: source-driven construction developed in this session. The equilibrium algebra and global potential inequality below have proofs. A separate agent certified existence of the endpoint pair with exact integer-interval Simpson bounds; another certified the energy upper bound J<-1.38 with Arb. The integral agent independently checked the endpoint-interval algorithm, Carlson reduction and energy Jacobian. These are reproducible computational certificates, not a formal proof-assistant verification. No claim about the validity of the entire Fauzan preprint follows from this local analytic investigation. Historical novelty unresolved.

Primary source read: A. Fauzan, the supplied preprint saved as `codex/sources/fauzan-zeta5.pdf`, Section 6 and Appendix A. Its Section 6 reduces the real determinant bound to a logarithmic-energy comparison, and Appendix A certifies that comparison with sixteen nested arcsine measures. The parent agent proposed the signed-charge/balayage construction below; the integral agent independently derived its density, Stieltjes transforms, sign argument and endpoint monotonicity.

## 1. The exact optimization problem

Keep Fauzan's parameters

    α=3/40, λ=37/40,
    V(x)=2πsqrt(x)+∫_0^1 log(x+u²)du−6∫_0^α log(x+u²)du,
    C*=-2λ+12αλ(1−logα)+3λ²−2λ²log(2λ).

For a positive compactly supported measure ρ of mass λ, set

    Uρ(x)=∫log|x-t|dρ(t),
    I(ρ)=∫∫log|x-t|dρ(x)dρ(t),
    M=sup_(x≥0)[2Uρ(x)−V(x)].

The real determinant argument needs an upper bound

    J(ρ)=λM−I(ρ)+C*.

The preprint supplies J≤-1.3669955 through a rounded rational certificate. Its arithmetic bound for M_arith=200 is approximately A_200=1.34958769774, so any rigorous J<-A_200 is sufficient for that particular real/arithmetic comparison, assuming the separate arithmetic proof is correct.

Why this is the right functional: if σ also has mass λ, zero-mass logarithmic-energy negativity gives

    I(σ)-∫Vdσ≤∫(2Uρ−V)dσ−I(ρ)≤λM−I(ρ).

Thus an equilibrium measure, for which 2Uρ−V is constant on its support and no larger outside, gives the sharp continuum bound. No determinant computation is needed in this step.

## 2. Concrete simpler finite trials

Before using the analytic construction, I numerically optimized weights on small subsets of the source's nested intervals. For fixed endpoints the objective is convex in the weights on the fixed-mass affine space: its maximum-potential part is convex, and -I is convex on zero-mass perturbations. A grid-based SLSQP search gave the following indicative values:

| Number of components | Best tested objective J | Interpretation |
|---:|---:|---|
| 1 | +1.2312 | Far too weak for the tested interval |
| 2 | -0.9598 | Insufficient |
| 3 | -1.1526 | Insufficient |
| 4 | -1.2446 | Insufficient |
| 5 | -1.2885 | Insufficient |
| 6 | -1.3266 | Potentially sufficient with a larger arithmetic cutoff, not M_arith=200 |
| 7 | -1.3429 | Potentially sufficient with a larger arithmetic cutoff |
| 8 | -1.3451 | Still insufficient for M_arith=200 in this restricted search |
| 16 | -1.3718 | Slightly improves the source weights numerically |

These are exploratory grid maxima, not certified upper bounds or lower bounds on what a given number of components could achieve. Endpoints were selected from the source table rather than optimized continuously. They show that simply dropping most arcsine components loses useful margin. The exact equilibrium route is conceptually more promising than presenting a different unexplained table.

## 3. Express the field derivative as signed logarithmic charges

For x>0,

    V'(x)=∫_0^∞ q(u)/(x+u²)du,

where

    q(u)=-3 for 0<u<α,
          3 for α<u<1,
          2 for u>1.                                (Q)

Indeed 2∫_0^∞du/(x+u²)=π/sqrt(x), the derivative of 2πsqrt(x). Thus the field has one negative-to-positive sign change in its charge density, at u=α. This simple sign pattern will replace the numerous interval checks in Appendix A.

## 4. Sweep each charge onto one interval

Fix provisional endpoints 0<a<b. Let μ be the unit arcsine measure on [a,b], and define

    C_u=sqrt((u²+a)(u²+b)),
    dB_u(t)=C_u/(t+u²) dμ(t).

The arcsine Stieltjes transform gives ∫dμ(t)/(t+u²)=1/C_u, so B_u has mass one. It is the interval balayage of the point at -u², but the displayed formula is all that is needed.

Define a signed measure of mass λ by

    ρ=λμ+(1/2)∫_0^∞ q(u)(B_u−μ)du.                  (R)

The integral converges in total variation. Uniformly for t∈[a,b], C_u/(t+u²)-1=O(u^-2) as u→∞, and the integral near u=0 is bounded because a>0. Each B_u−μ has mass zero, so ρ has mass λ.

Its density is

    dρ(x)=F(x) dx/[πsqrt((x-a)(b-x))],
    F(x)=λ+(1/2)∫q(u)[C_u/(x+u²)−1]du.             (F)

For arbitrary endpoints this density might change sign. The endpoint equations will guarantee positivity.

## 5. Choose endpoints so the density vanishes at both edges

Set

    I0(a,b)=∫_0^∞ q(u)/C_u du,
    H(a,b)=λ+(1/2)∫_0^∞ q(u)[sqrt((u²+a)/(u²+b))−1]du.

Require

    I0(a,b)=0,   H(a,b)=0.                           (E)

Here H=F(b), while F(a)-F(b)=(b-a)I0/2. Hence (E) is exactly F(a)=F(b)=0.

For x>0 define

    J(x)=∫_0^∞ q(u)/[(x+u²)C_u]du.

Subtracting 1/(x+α²) times I0 from this integral gives

    J(x)=∫ q(u)/C_u [1/(x+u²)−1/(x+α²)]du<0.        (S)

Both the u<α part and the u>α part are strictly negative. The sign-change property (Q), rather than a numerical potential grid, proves this strict inequality everywhere.

Linear interpolation between F(a)=F(b)=0 gives the exact identity

    F(x)=(1/2)(x-a)(x-b)J(x).

Therefore on the interval,

    dρ(x)=-sqrt((x-a)(b-x))J(x)/(2π) dx>0.           (D)

This is a positive measure of the required mass. Its density is bounded and vanishes like a square root at both endpoints; in particular its logarithmic energy is finite.

## 6. Prove the global potential inequality

Let Gμ(z)=∫dμ(t)/(z-t)=1/sqrt((z-a)(z-b)), where the square root behaves like z at infinity. Partial fractions give

    G_(B_u)(z)=[1+C_u Gμ(z)]/(z+u²).

Thus outside [a,b],

    2Gρ(x)−V'(x)=2Gμ(x)F(x)
                 =sqrt((x-a)(x-b))J(x),             (P)

with the same branch convention. Inside the interval the real part of Gμ is zero, while the real part of G_(B_u) is 1/(x+u²); hence 2Uρ−V is constant there. This can equivalently be proved by the standard explicit arcsine/balayage potential formulas.

For 0<x<a, the square root branch in (P) is negative, while J<0, so (2Uρ−V)' is positive. For x>b it is negative. Consequently

    2Uρ(x)−V(x)≤M for all x≥0,

with equality exactly on [a,b]. Continuity includes x=0 and both endpoints. This proves the desired global inequality symbolically once (E) is solved.

The measure is therefore the true equilibrium measure: the zero-mass energy inequality from Section 1 gives maximality, and its strict version gives uniqueness. This conclusion is stronger than a trial measure certificate, but no explicit elliptic-function formula is required.

## 7. Endpoint existence and uniqueness reduce to one monotone scalar equation

At a zero of I0, differentiation under the integral gives

    ∂_a I0=-J(a)/2>0,
    ∂_b I0=-J(b)/2>0.

For fixed b beyond the minimum of V, I0(a,b)→-∞ as a↓0, whereas

    I0(b,b)=V'(b)>0.

In the actual certified region b∈[.7,.8], positivity of V′(b) needs no numerical minimization: its numerator π+arctan(1/sqrt b)-6arctan(α/sqrt b) is at least π-6α/sqrt b>0, using arctan t≤t. Thus there is a root a(b)∈(0,b). Every root has positive a-derivative, so there can be only one: two upward crossings would require another zero that is not an upward crossing. The implicit function theorem gives

    a'(b)=-J(b)/J(a)<0.

For H above,

    ∂_a H=I0/4,
    ∂_b H=-(1/4)∫q(u)(u²+a)/[(u²+b)C_u]du.

Along I0=0, these simplify to

    d/db H(a(b),b)=(b-a)J(b)/4<0.                    (M)

Therefore one pair of rigorously certified opposite signs in b proves existence and uniqueness of the endpoint pair. A two-dimensional root existence argument is unnecessary.

Nonrigorous numerical bracketing obtained independently in SciPy is:

| b | a(b) | H(a(b),b) |
|---:|---:|---:|
| .700 | .00007547460057 | +.03683203422 |
| .750 | .00007148890338 | +.002644375239 |
| .753 | .00007126549678 | +.000635071458 |
| .754 | .00007119139631 | -.000033680771 |
| .760 | .00007075061255 | -.004035602807 |
| .800 | .00006796984617 | -.03026560738 |

The numerical root is a≈.0000711951252378642, b≈.7539496182816022. These decimals are approximations, not endpoint definitions. Rigorous endpoints are the unique solution of (E). The exact interval script `codex/computations/certify_equilibrium_endpoints.py` certifies opposite signs on all four faces of the rational box a∈[.00007119508,.00007119517], b∈[.75394957,.75394967], proving existence there; the monotonicity argument gives uniqueness.

## 8. A simpler formula for the norm constant

The interval potential at the right edge satisfies

    Uμ(b)=log((b-a)/4),
    U_(B_u)(b)−Uμ(b)
       =2log[2sqrt(u²+b)/(sqrt(u²+a)+sqrt(u²+b))].

Hence

    Uρ(b)=λlog((b-a)/4)
       +∫q(u)log[2sqrt(u²+b)/(sqrt(u²+a)+sqrt(u²+b))]du,
    M=2Uρ(b)−V(b).                                  (N)

This avoids numerically integrating the logarithmic singularity log(b-x) against the density. Since 2Uρ−V=M on the support,

    2I(ρ)−∫Vdρ=λM,
    J(ρ)=(λM−∫Vdρ)/2+C*.                            (J)

The remaining moment ∫Vdρ is an ordinary integral with bounded density (D). Use x=(a+b)/2+(b-a)cosθ/2 to remove square-root endpoints; the measure factor becomes -(b-a)²sin²θ J(x)/(8π)dθ.

The parent agent's independent Gauss–Legendre calculations at 128,256,512 nodes gave

    mass≈.925000000077,
    M≈-6.660781123818,
    I(ρ)≈-2.122676765675,
    J(ρ)≈-1.385509783517.

This would improve the source's certified real exponent by about .0185 and would leave about .0359 margin against A_200. The high-precision approximations above are now supplemented by a validated enclosure: `codex/computations/certify_equilibrium_energy.py` uses the certified endpoint box and Arb arithmetic to prove J<-1.38. Its saved JSON reports J in [-1.3855±.0000964], safely below -A_200. Agreement of node counts was only exploratory; the interval enclosure supplies the certificate.

## 9. The scalar certificate and its scope

After the functional proof above, only scalar enclosures are required. They have now been supplied as follows:

1. `certify_equilibrium_endpoints.py` proves four entire-face signs with exact dyadic interval arithmetic, integer square roots and composite Simpson remainders bounded by interval Taylor jets. This locates the unique endpoint pair.
2. `certify_equilibrium_energy.py` evaluates the norm expression using Carlson RF/RJ formulas for J(x), analytic interval quadrature on six transformed angular segments, and a closed plateau formula. It proves the rational target J<-69/50=-1.38; -1.38+A_200<-0.0304. The exact rational arithmetic constant A_200 is supplied by the separate arithmetic argument.
3. Record a quantitative bounded-density constant if an explicit finite-K estimate is desired. For the limiting K² coefficient, a fixed finite bound suffices: circle regularization changes the potential by O(ε), because a bounded density assigns O(r) mass to an interval of radius r. Taking ε=K^-2 leaves only lower-order errors. The treatment of the unbounded integration domain is unchanged in principle, since V(x)=2πsqrt x+O(log x), while Uρ(x)=λlog x+O(1/x).

For endpoint integrals, the substitution v=1/u sends the tail u≥1 to [0,1]. For example,

    ∫_1^∞du/C_u=∫_0^1dv/sqrt((1+av²)(1+bv²)),

and the negative mass-integrand tail becomes

    (a-b)∫_0^1 dv/[sqrt(1+bv²)(sqrt(1+av²)+sqrt(1+bv²))].

These are smooth, bounded and monotone. The finite pieces of 1/C_u and 1-sqrt((u²+a)/(u²+b)) are also monotone. Thus even simple Darboux bounds can certify coarse endpoint signs; validated high-order quadrature can then bound the energy more efficiently.

The theoretical simplification is substantial: one sign change in q proves positivity and the global potential inequality, replacing a sixteen-component numerical shape and its many potential cells. It uses more potential-theoretic language, but its key identities reduce to arcsine integrals and partial fractions. The full claimed irrationality theorem still depends on all separate arithmetic and moment results, which are not certified by this note.


## 10. Two further reductions of the scalar calculations

The connections agent derived the following reductions, and the integral agent independently checked them.

**Closed plateau constant.** Put cap=(b-a)/4, V0=-2+12α(1-logα), and

    K(u)=log[4u²/(sqrt(u²+a)+sqrt(u²+b))²].

Then

    M=2λlog(cap)-V0+K(1)-6αK(α)-4λ.                 (N2)

Indeed M=2λlog(cap)-V0+∫qK. Since K'=2/u-2u/C_u, integration by parts across the jumps of q gives ∫qK=K(1)-6αK(α)-2∫q(1-u²/C_u). The endpoint equations make the last integral equal 2λ. This removes the quadrature in (N), and is the formula used in the completed Arb certificate.

**One-dimensional energy alternative.** The equilibrium force equation, multiplied by x and integrated symmetrically against ρ, gives the virial identity

    ∫xV'(x)dρ(x)=λ².

The identity follows by symmetrizing x/(x-y) with y/(y-x); their sum is one. The bounded square-root density justifies taking the symmetric principal-value limit. Direct integration of V gives

    V(x)-2xV'(x)=log(x+1)-6αlog(x+α²)-2+12α.

Consequently

    ∫Vdρ=2λ²+Uρ(-1)-6αUρ(-α²)+λ(-2+12α).          (VIR)

For s>0 define

    D_s=(sqrt(s+a)+sqrt(s+b))²/4,
    τ_s=(sqrt(s+b)-sqrt(s+a))/(sqrt(s+b)+sqrt(s+a)).

The arcsine potential and the Poisson-kernel form of B_u give

    Uρ(-s)=λlog D_s+∫q(u)log(1-τ_s τ_(u²))du.        (OFF)

For verification, write x=(a+b)/2+(b-a)cosθ/2. The measure B_u is a Poisson kernel with parameter -τ_(u²), while log(s+x)=log D_s+2 Re log(1+τ_s e^(iθ)). Averaging its absolutely convergent Fourier series gives U_(B_u)(-s)-Uμ(-s)=2log(1-τ_sτ_(u²)), yielding (OFF).

Thus (VIR) reduces the energy moment to two smooth one-dimensional u-integrals, avoiding Carlson functions if desired. A separate agent has now certified this second implementation as well: `certify_equilibrium_energy_elementary.py` and its JSON give ∫Vdρ∈[1.9159±.0000703] and J∈[-1.3855±.0000295], again proving J<-1.38 and overlapping the Carlson result. This independent implementation uses three smooth scalar integrals rather than elliptic functions.

## 11. Independent check of the Carlson formula used in the certificate

For x>0, T≤∞, put s=T/sqrt(a+T²), c²=1-s², with s=1,c²=0 at infinity. The substitution u=sqrt(a)tanθ gives

    ∫_0^T du/[(x+u²)sqrt((u²+a)(u²+b))]
      =1/(x sqrt b) [s RF(c²,1-(1-a/b)s²,1)
       -(a/x)s³ RJ(c²,1-(1-a/b)s²,1,c²+(a/x)s²)/3].

This follows by splitting cos²θ/[x cos²θ+a sin²θ] as x^-1[1-(a/x)sin²θ/(1-(1-a/x)sin²θ)]. All real arguments of the Carlson functions are nonnegative, with the required positive arguments strictly positive. Combining the values at T=∞,1,α with coefficients 2,1,-6 gives J(x). The certificate's angular substitution x=a+(b-a)sin²θ has density factor -(b-a)²sin²θ cos²θ J(x)/π, which matches (D) exactly.


## 12. General existence is now analytic

`general-equilibrium-family.md` proves existence and uniqueness of the endpoints for every r>3, 0<α<1 and λ>0. Thus the four-face numerical certificate is now only a location certificate for the unique endpoints, not a prerequisite for existence. The same note gives the explicit density bound max(r-3,3)(b-a)/(16a^(3/2)), sufficient for the circle regularization error.
