# A rigorous scalar certificate for the exact equilibrium candidate

Status: two rigorous Arb enclosures of the same energy objective, using independently certified endpoint bounds. The first reduces the density to Carlson elliptic functions; the second uses elementary potential and virial identities. Both certify an upper bound below -1.38. This verifies the numerical inequality required by the proposed analytic simplification; the equilibrium construction and endpoint existence are separate proved inputs in the companion files. No historical novelty claim is made, and this does not independently audit every other part of Fauzan's irrationality proof.

## 1. The object being certified

Use alpha=3/40, lambda=37/40, and

    V(x)=2πsqrt(x)+integral_0^1 log(x+u²)du
                         -6 integral_0^alpha log(x+u²)du.

The exact equilibrium candidate from `fauzan-analytic-simplification.md` has support [a,b], mass lambda, and density

    dρ(x)/dx= -sqrt((x-a)(b-x)) J(x)/(2π),
    J(x)=integral_0^infinity q(u)/[(x+u²)sqrt((u²+a)(u²+b))]du,

where q=-3 on [0,alpha], q=3 on [alpha,1], and q=2 on [1,infinity). Its endpoint equations are checked in `equilibrium-endpoint-certificate.md`. The certified rectangle is

    0.00007119508 <= a <= 0.00007119517,
    0.75394957 <= b <= 0.75394967.

The energy objective is

    U = lambda*M/2 - (integral V dρ)/2 + Cstar,
    Cstar=-2lambda+12alpha*lambda*(1-log alpha)
                    +3lambda²-2lambda²log(2lambda),

where 2U^ρ(x)-V(x)=M on the support. The relation follows from integrating the plateau identity: 2I(ρ)-integral V dρ=lambda M. Consequently lambda M-I(ρ)+Cstar equals the displayed objective. This is the same analytic quadratic coefficient used in the determinant estimate.

The exact endpoint proof is essential. Merely evaluating a density at approximate endpoints would not establish the required mass, positivity, or global potential inequality. The computations below enclose the functions over the entire certified rectangle, so they include the exact endpoints supplied by that proof.

## 2. The plateau constant needs no numerical quadrature

Define

    K(u)=log[4u²/(sqrt(u²+a)+sqrt(u²+b))²],
    V0=-2+12alpha*(1-log alpha).

The endpoint equations and integration by parts, as derived independently by the parent, integral, and connections agents, give

    M=2lambda log((b-a)/4)-V0+K(1)-6alpha K(alpha)-4lambda.

Direct interval evaluation on the certified rectangle gives

    M in [-6.660781 +/- 0.000000286].

All logarithms and square roots here have strictly positive real arguments.

## 3. First route: one scalar integral with a Carlson formula for the density

Put

    I_T(x)=integral_0^T du/[(x+u²)sqrt((u²+a)(u²+b))].

For finite T, set sigma=T/sqrt(a+T²), c²=1-sigma², and m=1-a/b. Substitution u=sqrt(a)tan(theta) gives

    I_T(x)=1/(x sqrt(b)) * {
        sigma RF(c²,1-m sigma²,1)
        -(a/x) sigma³ RJ(c²,1-m sigma²,1,c²+(a/x)sigma²)/3
    }.

For T=infinity use sigma=1,c²=0. RF and RJ are the standard Carlson symmetric elliptic integrals. To verify the formula, the substitution gives the incomplete first-kind elliptic integral minus (a/x) times the sine-squared third-kind contribution. Their Carlson forms are respectively sigma RF and sigma³ RJ/3. No numerical fitting is involved. The integral agent independently reconstructed this reduction.

Then

    J(x)=2 I_infinity(x)+I_1(x)-6 I_alpha(x).

Writing x=a+(b-a)sin²(theta) removes both square-root endpoint factors, so

    integral V dρ
      =-(b-a)²/π * integral_0^(π/2)
          sin²(theta)cos²(theta) V(x(theta))J(x(theta))dtheta.

The elementary closed form used for V is

    integral_0^T log(x+u²)du
      =T log(x+T²)-2T+2sqrt(x)atan(T/sqrt(x)).

`certify_equilibrium_energy.py` encloses this integral using six intervals with breakpoints 0,.01,.03,.1,.3,.7,π/2. Arb performs rigorous adaptive complex quadrature at160-bit precision. The callback explicitly rejects complex balls with Re x<=0; it passes the analyticity flag to square roots and logarithms, checks the arctangent argument remains in the right half-plane, and checks the variable fourth Carlson argument has positive real part. The other Carlson arguments are fixed nonnegative real balls, with only the integrable first argument equal to zero in the complete case.

The resulting enclosures are

    integral V dρ in [1.916 +/- 0.000305],
    U in [-1.3855 +/- 0.0000964] < -1.38.

The complex quadrature's imaginary enclosure contains zero, as required for the real integral. A requested tolerance is not itself treated as a guarantee: the code checks the actual returned finite ball and the final strict inequality.

## 4. Independent route: elementary potential and virial identities

There is a separate reduction avoiding elliptic functions. From the equilibrium derivative equation, symmetric integration gives the virial identity

    integral x V'(x)dρ(x)=lambda².

One can justify the principal-value symmetrization by excluding |x-y|<=epsilon, swapping x,y, and then letting epsilon decrease to zero. The symmetrized kernel is constant, and the nonatomic density makes the excluded diagonal mass vanish.

Direct differentiation of V gives

    V(x)-2xV'(x)=log(x+1)-6alpha log(x+alpha²)-2+12alpha.

Therefore

    integral V dρ=2lambda²+U^ρ(-1)-6alpha U^ρ(-alpha²)
                                      +lambda(-2+12alpha).

For s>0 define

    D_s=(sqrt(s+a)+sqrt(s+b))²/4,
    tau_s=(sqrt(s+b)-sqrt(s+a))/(sqrt(s+a)+sqrt(s+b)).

The balayage potential formula gives

    U^ρ(-s)=lambda log D_s
              +integral_0^infinity q(u)log(1-tau_s tau_(u²))du.

Combining the s=1 and s=alpha² terms leaves three smooth one-dimensional integrals: over [0,alpha], over [0,1], and over the tail [1,infinity). On the tail use v=1/u. The apparent v^-2 singularity is removable because tau_(u²)=O(u^-2). The code evaluates it stably through

    2F1(1,1;2;z)=-log(1-z)/z,

including its value1 at z=0. This use of a special-function evaluator merely handles an elementary removable logarithmic quotient. Its branch is checked by Re z<1; all original real arguments are in (0,1).

`certify_equilibrium_energy_elementary.py` gives

    integral V dρ in [1.9159 +/- 0.0000703],
    U in [-1.3855 +/- 0.0000295] < -1.38.

The second enclosure overlaps the first, as checked programmatically. This is an independent mathematical reduction of the integral, although both numerical certificates use the same Arb library. It is therefore not an independent implementation of interval arithmetic.

## 5. What this simplifies, and what it leaves intact

The exact equilibrium construction turns the global potential comparison into analytic identities and two endpoint equations. Once those equations are rigorously enclosed, the needed energy inequality follows from a handful of scalar integrals. The published sixteen-component arcsine approximation and its large finite partition need not be reproduced for this candidate analytic route.

The result is an upper bound for the analytic quadratic coefficient, not yet a replacement for the complete irrationality proof. Retaining the source's discrete-energy regularization, its lower-order error estimates, and its prime-adic arithmetic normalization requires checking their hypotheses for this new measure. The parent and integral agents address that transfer separately. In particular the source's special arcsine description cannot simply be deleted where it was used to establish logarithmic regularity.

This note certifies U<-1.38. It does not infer historical novelty, claim a formal Lean verification, or silently extend the result to other odd zeta values.

## Reproduction

Run

    PYTHONPATH=/workspace/scratch/odd-zeta-python python codex/computations/certify_equilibrium_energy.py
    PYTHONPATH=/workspace/scratch/odd-zeta-python python codex/computations/certify_equilibrium_energy_elementary.py

The scratch installation contains python-flint0.9.0, exposing FLINT/Arb interval arithmetic. Elsewhere install that package normally. Both scripts save JSON enclosures alongside themselves and run independently; the elementary script makes no Carlson-function calls. It checks overlap with the first route's saved outward-rounded real interval. The certified endpoint box is independently produced by `certify_equilibrium_endpoints.py`, which uses exact integer interval arithmetic rather than Arb.
