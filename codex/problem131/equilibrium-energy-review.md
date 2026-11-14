# Independent audit of the equilibrium energy certificate

Status: independently reviewed the Carlson reduction, interval parameter inclusion, analytic callback domains, plateau formula, objective assembly and exact comparison with the source arithmetic constant. Reran the rigorous energy certificate successfully. No substantive error was found in these reviewed components. This is an audit of the replacement analytic certificate, not acceptance of every arithmetic theorem in Fauzan's preprint or a standalone claim to have proved ζ(5) irrational.

Primary source consulted directly: `codex/sources/fauzan-zeta5.txt`, Sections6–7 and AppendixB.3. The proposed measure and its global potential inequality are developed in `fauzan-analytic-simplification.md` and independently reviewed by the measure agent. Endpoint existence is proved separately in `equilibrium-endpoint-certificate.md`.

## 1. Parameter inclusion is exact

The unique endpoint solution is certified to lie in

    a∈[0.00007119508,0.00007119517],
    b∈[0.75394957,0.75394967].

The energy code constructs Arb balls with centers and radii

    a=0.000071195125 ±0.000000000045,
    b=0.75394962 ±0.00000005.

Their mathematical intervals are exactly the same rectangle. Arb includes its own rounding errors outward. I checked programmatically that the two exact rational endpoints of each certified interval are contained in the corresponding Arb ball, and that b-a>0. All formulas are evaluated over these full balls.

Some formulas use the endpoint equations and need only be correct at their unique root. Evaluating those formulas over the full parameter rectangle still encloses their value at that root; it is not necessary to assert that the formulas represent an equilibrium measure at every off-root point in the rectangle. Ignoring correlations between a and b widens the enclosure and does not invalidate it.

## 2. The Carlson formula has the correct numerator and normalization

The code needs

    J_T(x)=integral_0^T du/[(u²+x)sqrt((u²+a)(u²+b))].

For real x>0, substitute u=sqrt(a)tan θ. Put

    s=sin θ_T=T/sqrt(a+T²),
    c²=1-s²,
    m=1-a/b,
    h=1-a/x.

The transformed measure is

    du/sqrt((u²+a)(u²+b))
      =dθ/[sqrt(b)sqrt(1-m sin²θ)],

while

    1/(u²+x)=(1/x)[1-(a/x)sin²θ/(1-h sin²θ)].

The standard defining integral forms of Carlson's symmetric functions give

    integral_0^θ dφ/sqrt(1-m sin²φ)
       =s R_F(c²,1-ms²,1),
    integral_0^θ sin²φ dφ/
        [(1-h sin²φ)sqrt(1-m sin²φ)]
       =(s³/3)R_J(c²,1-ms²,1,1-hs²).

Consequently

    J_T(x)={sR_F-(a/x)(s³/3)R_J}/[x sqrt(b)].         (1)

This is precisely `Jpart` in the script, including the coefficient a/x multiplying R_J. For T=infinity, s=1 and c²=0; the other arguments remain positive, and the defining integrals converge. The charge density is q=2+1_[0,1]-6·1_[0,α], so

    J(x)=2J_infinity(x)+J_1(x)-6J_α(x),

also matching the code.

## 3. The density integration factor is correct

At the endpoint solution, the measure density is

    dρ(x)= -sqrt((x-a)(b-x))J(x)/(2π) dx.

Under x=a+(b-a)sin²θ, 0<=θ<=π/2,

    sqrt((x-a)(b-x))=(b-a)sinθcosθ,
    dx=2(b-a)sinθcosθdθ.

Thus

    integral V dρ=-(b-a)²/π integral_0^(π/2)
                      sin²θcos²θ V(x)J(x)dθ.           (2)

This agrees with `density_integrand`. The factor is not twice or half the required value.

The code's closed form

    integral_0^u log(x+t²)dt
      =u log(x+u²)-2u+2sqrt(x)atan(u/sqrt(x))

has the correct derivative in u and vanishes at u=0. Therefore `V` represents exactly the source field

    V(x)=2πsqrt(x)+integral_0^1 log(x+u²)du
                          -6integral_0^α log(x+u²)du.

## 4. Analytic callback requirements were checked

Arb's validated integrator requires its callback to reject a complex integration ball on which the selected expression is not analytic. The script does so conservatively:

- The transformed x must have strictly positive real part. This excludes the branch cuts of sqrt(x) and log(x+u²), and avoids x=0.
- In `V`, both square root and logarithm receive the integrator's `analytic` flag. The argument u/sqrt(x) is also required to have strictly positive real part, which excludes the imaginary-axis cuts of atan. Thus calling atan without a separate analytic flag is valid on the explicitly checked domain.
- In `Jpart`, the first three Carlson arguments are fixed positive real constants, except that the first is0 for the convergent infinite-tail case. The varying fourth argument c²+(a/x)s² is required to have positive real part. Carlson R_J is analytic there, as follows directly from its defining integral with denominator t+fourth. R_F is constant in the integration variable.
- If these checks fail, the script returns a non-finite ball, forcing subdivision rather than silently accepting a branch-crossing disk.

The interval parameters a,b are real and positive. Their widths are propagated through all operations. Every family member is analytic on the accepted domain; interval overestimation of the whole family is conservative.

Finite tolerance settings are goals, not claims of unproved accuracy: the result is an Arb enclosure. The code checks that every segment result is finite, sums the balls, and tests the final interval inequality. Its complex imaginary enclosure contains0, as required for the real integral. Taking the real component retains a valid enclosure of its exact real value.

The six integration segments end at the enclosed value π/2. Arb supports ball-valued endpoints; this is also how its own documented examples integrate to π. Endpoint rounding is therefore part of the validated integration, not an ignored decimal approximation.

## 5. Independently derived plateau formula

The code uses a simplified M formula not obtained by sampling the maximum potential. Define

    κ(u)=log[4u²/(sqrt(u²+a)+sqrt(u²+b))²],
    V0=V(0)=-2+12α(1-logα).

From the balayage formula at the right endpoint,

    M=2λlog((b-a)/4)-V0+integral_0^infinity q(u)κ(u)du.

To verify this, subtract V(b)-V(0)=integral q log((b+u²)/u²)du from the right-edge potential formula. The two logarithms combine to κ exactly.

Now

    κ'(u)=2/u-2u/sqrt((u²+a)(u²+b)).

Integration by parts on the three constant-charge intervals gives

    integral qκ =κ(1)-6ακ(α)
             -2integral q(u)[1-u²/sqrt((u²+a)(u²+b))]du.

The boundary terms uκ(u) vanish at0 and infinity. At the endpoint solution, the mass equation and I0=0 make the last integral equal2λ. Therefore

    M=2λlog((b-a)/4)-V0+κ(1)-6ακ(α)-4λ.             (3)

This is exactly the code's expression. No sign or factor discrepancy was found.

On the support, 2Uρ-V=M. Integrating against a measure of massλ gives

    2I(ρ)-integral Vdρ=λM.

Hence the desired objective is

    λM-I(ρ)+C* = λM/2-(1/2)integral Vdρ+C*,         (4)

which verifies the final assembly.

## 6. Independent execution and certified result

I reran

    PYTHONPATH=/workspace/scratch/odd-zeta-python \
      python codex/computations/certify_equilibrium_energy.py

at160-bit precision. The script completed successfully and returned enclosing balls displayed as

    M=[-6.660781 +/-2.86e-7],
    integral Vdρ=[1.916 +/-3.05e-4]+[+/-1.63e-4]i,
    C*=[2.6530359903404886611292503361359534827071655799 +/-3.57e-47],
    objective=[-1.3855 +/-9.64e-5].

The actual Arb comparison `objective < -1.38` passed. These rounded displays are not substituted for the underlying ball comparison. The source JSON records the segment results, precision and parameters.

Thus, combining the endpoint theorem, the independently derived positive measure/global potential inequality, and this scalar enclosure, the comparison measure satisfies

    λM-I(ρ)+C* < -69/50 = -1.38.                    (5)

This is stronger than the preprint's rational certificate U=-2733991/2000000=-1.3669955.

## 7. Replacement statement for the real determinant estimate

Use the source definitions of F_K=S_KΔ_K and K∈40Z_>0, together with its moment representation (6.10), scaling estimate (6.13), and factorial normalization (6.15). The new comparison measure yields

    0<F_K(ζ(5)),
    log F_K(ζ(5)) <= -(69/50)K²+o(K²).              (6)

Here is why the different density does not create a regularization gap. It is bounded on its compact interval [a,b]⊂(0,1): J(x) is continuous there and the density vanishes at the endpoints. If L bounds the density, its mass in any interval of radius r is at most2Lr. The circle-average identity therefore bounds the potential regularization error by

    integral_0^ε [2Lr]/r dr=2Lε.

With ε=K^-2, the cross-energy error after multiplication by K² is O(1), and the circle self-energy contributes O(KlogK). All logarithmic integrability requirements in the source's zero-mass energy lemma remain satisfied.

The unbounded integration domain also remains controlled. For t>=2, the support and mass bound give the same coarse source estimate

    2Uρ(t)-V(t) <=(13/10)log t+6α³/t-2πsqrt(t).

For K>=2, adding sqrt(t)/K gives a decreasing upper bound on[2,infinity). Its value at2 is below-6.78 by elementary bounds log2<0.7, π>3, sqrt2>1.4. The certified M is greater than-6.661, so this tail bound lies strictly below M. On0<=t<=2, the perturbation is at mostsqrt2/K. Thus the source's retained integrable factor exp(-sqrt t) and its scaling argument continue to apply with only lower-order errors.

The K²logK terms cancel exactly as in the source, leaving the coefficient in(5). This proves(6) with the stated moment/scaling formulas as explicit inputs. An explicit all-K bound with the original numerical linear-error constants is not asserted; those constants depend on the comparison density. The asymptotic K² coefficient is fully adequate for the intended comparison.

## 8. Exact comparison with the source's A_200

AppendixB.3 gives the exact rational number

    A_200=127125602969131786927559/94195881588024216000000.

Direct rational subtraction gives

    69/50-A_200
       =2864713622341631152441/94195881588024216000000
       ≈0.030412302258295783 >3/100.                  (7)

Therefore, **if** the separate arithmetic theorem establishing integer polynomials Q_(K,200) and the normalization growth bound A_200 is accepted, the new analytic certificate gives

    limsup K^-2 log Q_(K,200)(ζ(5))
       <=A_200-69/50 < -3/100.

At K=40n the exact quadratic margin is approximately48.6596836n². It is more than enough to beat the factor b^(37n) under a rationality assumption ζ(5)=a/b.

This last implication is conditional on those arithmetic/integrality results and the source's degree/positivity statements. The present review verifies the analytic replacement and numerical margin; it does not independently certify all p-adic estimates in the preprint. Accordingly it is not presented as a completed independent proof of ζ(5) irrationality.

## Reproducibility

Run `check_equilibrium_energy_review.py` with python-flint available (the environment path above). It reruns the endpoint and energy certificates, verifies inclusion of the exact endpoint rationals, and checks the exact A_200 margin. The mathematical branch-domain and Carlson derivations are the proofs in this review, not conclusions inferred from matching numerical samples.
