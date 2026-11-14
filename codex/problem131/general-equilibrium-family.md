# An exact equilibrium family with one sign change in the external charges

Status: general theorem proposed by the parent agent and independently derived/checked by the integral agent. The proof below is analytic and establishes existence and uniqueness of the endpoint pair without numerical root certificates. The previously certified endpoint box is needed only to locate the particular Fauzan measure and enclose its energy. Historical novelty is unresolved.

## 1. The family and theorem

Fix real parameters

    r>3,  0<α<1,  λ>0,

and on x≥0 define

    V(x)=2πsqrt x+∫_0^1 log(x+u²)du−r∫_0^α log(x+u²)du.

At zero use the continuous value

    V(0)=-2+2rα(1-logα).

Let

    q(u)=3-r  for 0<u<α,
          3   for α<u<1,
          2   for u>1.

Equivalently q=2+1_(0,1)-r1_(0,α). Then V'(x)=∫_0^∞q(u)/(x+u²)du for x>0.

**Theorem.** There is a unique pair 0<a<b solving

    I0(a,b)=∫_0^∞ q(u)/sqrt((u²+a)(u²+b))du=0,          (E1)
    H(a,b)=λ+(1/2)∫_0^∞q(u)
                    [sqrt((u²+a)/(u²+b))-1]du=0.      (E2)

For these endpoints the density

    dρ(x)=-sqrt((x-a)(b-x))/(2π)
          ·∫_0^∞ q(u)/[(x+u²)sqrt((u²+a)(u²+b))]du dx,
          a<x<b,                                    (D)

defines a positive measure of mass λ. Its logarithmic potential satisfies

    2Uρ(x)-V(x)=M on [a,b],
    2Uρ(x)-V(x)<M on [0,∞)\[a,b].                   (P)

Consequently it is an exact equilibrium comparison measure for the logarithmic-energy problem, rather than a discretized trial measure.

The endpoint statement is proved in Sections 2–5. The measure and potential statement is proved in Section 6. All integrals in (E1)–(E2) converge absolutely for positive a,b; their integrands are O(u^-2) at infinity.

## 2. The field has exactly one minimum

Put Φ(y)=V(y²) for y>0. Direct differentiation gives

    Φ'(y)=2π+2arctan(1/y)-2r arctan(α/y),
    Φ''(y)=2[(rα-1)y²+rα-α²]/[(1+y²)(α²+y²)].       (F)

Also

    Φ'(0+)=(3-r)π<0,   Φ'(+∞)=2π>0.

If rα≥1, then Φ''>0 everywhere, so Φ' has exactly one zero. If rα<1, Φ'' is positive and then negative, with one sign change. On its decreasing tail Φ' is strictly above its positive limit 2π. Hence its unique zero lies on the increasing part. In both cases Φ has exactly one minimum, at y_0>0. Set t_0=y_0². Then

    V'(x)<0 for 0<x<t_0,   V'(x)>0 for x>t_0.

This step uses only elementary differentiation and monotonicity; no numerical location of t_0 is required.

## 3. For each b>t_0, solve the first endpoint equation uniquely

For 0<a≤b define

    C_u=sqrt((u²+a)(u²+b)),
    J(x;a,b)=∫_0^∞q(u)/[(x+u²)C_u]du.

At any zero of I0, for every x>0,

    J(x;a,b)=∫ q(u)/C_u
                    [1/(x+u²)-1/(x+α²)]du<0.        (S)

Indeed q is negative below α and positive above it, while the bracket has the opposite sign. Both regions have positive length. This is the basic sign mechanism behind the theorem.

At a zero of I0, differentiating its absolutely convergent integral gives

    ∂_a I0=-J(a;a,b)/2>0,
    ∂_b I0=-J(b;a,b)/2>0.                            (I)

Now fix b>t_0. As a↓0, the negative part near u=0 diverges logarithmically, while the positive part stays bounded, so I0(a,b)→-∞. At the other endpoint,

    I0(b,b)=∫q(u)/(u²+b)du=V'(b)>0.

The intermediate value theorem gives a root a(b)∈(0,b). Every root has strictly positive a-derivative by (I), so every zero is an upward crossing. There cannot be two such roots: returning from positive values after the first crossing to negative values before the second would force a zero that is not an upward crossing. Thus a(b) is unique.

For b≤t_0 there is no root with a<b. The same upward-crossing argument, together with I0(b,b)=V'(b)≤0, rules it out. When b=t_0 the boundary zero a=b has positive a-derivative, so I0 is negative just to its left; an earlier upward crossing would again require a forbidden downward crossing.

The implicit function theorem and (I) show that a(b) is continuously differentiable on b>t_0 and

    a'(b)=-J(b;a,b)/J(a;a,b)<0.                       (A)

## 4. The branch starts at the field minimum

As b↓t_0, monotonicity implies that a(b) has a limit a_*≤t_0. That limit is positive: for any fixed b_1>t_0, a(b)≥a(b_1)>0 whenever t_0<b<b_1. Continuity of I0 gives I0(a_*,t_0)=0.

Section 3 ruled out a<t_0 at b=t_0, so a_*=t_0. Therefore

    a(b)→t_0 as b↓t_0,
    0<a(b)<t_0 for b>t_0.                            (B)

When a and b coalesce at this positive value, the integrand in (E2) tends to zero and is dominated by an integrable function. Hence

    H(a(b),b)→λ>0 as b↓t_0.                          (H0)

## 5. The second equation crosses zero exactly once

Differentiate H under the integral:

    ∂_a H=I0/4,
    ∂_b H=-(1/4)∫q(u)(u²+a)/[(u²+b)C_u]du.

Along the branch I0=0, this gives

    d/db H(a(b),b)=(b-a)J(b;a,b)/4<0.                (H1)

It remains to show that H tends to -∞. Write

    J_T(a,b)=∫_0^T[sqrt((u²+a)/(u²+b))-1]du.

Using q=2+1_(0,1)-r1_(0,α),

    H=λ+J_∞+(1/2)J_1-(r/2)J_α.

The finite-interval terms are bounded because their integrand lies between -1 and 0. For the infinite term, substitute u=sqrt b·v:

    J_∞/sqrt b
       =∫_0^∞[sqrt((v²+a/b)/(v²+1))-1]dv.

By (B), a/b→0 as b→∞. For 0≤a/b≤1, the absolute value of the integrand is bounded by

    1-v/sqrt(1+v²),

which is integrable on [0,∞), with integral 1. Dominated convergence therefore gives

    J_∞/sqrt b→-1,
    H(a(b),b)→-∞.                                   (H∞)

Combining (H0), strict monotonicity (H1), and (H∞), there is exactly one b>t_0 with H=0. Its unique a(b) proves the endpoint existence and uniqueness theorem for every allowed r,α,λ. No interval-arithmetic topology argument is needed for this general existence result.

## 6. Positivity and the global potential inequality

Let μ be the unit arcsine measure on [a,b]. For u≥0 define the probability measure

    dB_u(x)=C_u/(x+u²)dμ(x).

Its total mass is one because the arcsine integral of (x+u²)^(-1) is C_u^(-1). The signed integral

    ρ=λμ+(1/2)∫q(u)(B_u-μ)du

converges in total variation and has mass λ. Its density relative to μ is

    F(x)=λ+(1/2)∫q(u)[C_u/(x+u²)-1]du.

The endpoint equations imply F(a)=F(b)=0. Subtracting the linear interpolation of those two zero values gives

    F(x)=(1/2)(x-a)(x-b)J(x;a,b).

On a<x<b, this is positive by (S), and it gives exactly density (D). The density is bounded and vanishes at the endpoints, so its logarithmic potential is continuous and its energy is finite. Quantitatively, with Q=max(r-3,3), x≥a gives |J(x)|≤Q∫_0^∞du/(u²+a)²=Qπ/(4a^(3/2)). Hence its density is at most Q(b-a)/(16a^(3/2)). A circle regularization of radius ε increases its logarithmic potential by at most twice this bound times ε, since ∫_0^ε log(ε/u)du=ε. This supplies an explicit lower-order regularization constant without another numerical integral.

Let Gμ(z)=1/sqrt((z-a)(z-b)), with the branch asymptotic to 1/z. Direct partial fractions give

    G_(B_u)(z)=[1+C_uGμ(z)]/(z+u²).

Thus outside the support,

    (2Uρ-V)'(x)=sqrt((x-a)(x-b))J(x;a,b),             (W)

where the square root is negative for x<a and positive for x>b. Since J<0, the potential difference increases to the support from the left and decreases on the right. Inside the support, the real part of the arcsine Stieltjes transform vanishes, giving (2Uρ-V)'=0. Continuity includes x=0 and the endpoints, proving (P).

For any positive compactly supported measure σ of mass λ and finite logarithmic energy, zero-mass energy negativity yields

    I(σ)-∫Vdσ≤λM-I(ρ).

Equality holds for σ=ρ. Thus the measure solves precisely the comparison problem needed in Fauzan's determinant argument. Standard strict negativity gives uniqueness of the equilibrium measure; the endpoint pair was already proved unique independently.

## 7. General closed formulas for the plateau and energy

The scalar reductions also extend unchanged. Put

    cap=(b-a)/4,
    V0=-2+2rα(1-logα),
    K(u)=log[4u²/(sqrt(u²+a)+sqrt(u²+b))²].

Integration by parts across the jumps r at α and -1 at 1 gives

    M=2λlog cap-V0+K(1)-rαK(α)-4λ.                  (M)

The endpoint equations were used to replace ∫q(1-u²/C_u) by 2λ. The general virial relation and field identity are

    ∫xV'(x)dρ(x)=λ²,
    V(x)-2xV'(x)=log(x+1)-rαlog(x+α²)-2+2rα.

Consequently

    ∫Vdρ=2λ²+Uρ(-1)-rαUρ(-α²)+λ(-2+2rα).           (V)

For s>0,

    Uρ(-s)=λlog D_s+∫q(u)log(1-τ_sτ_(u²))du,
    D_s=(sqrt(s+a)+sqrt(s+b))²/4,
    τ_s=(sqrt(s+b)-sqrt(s+a))/(sqrt(s+b)+sqrt(s+a)).

These formulas reduce the real-energy certificate to smooth one-dimensional scalar integrals. The original Fauzan choice r=6,α=3/40,λ=37/40 has independently certified energy bound below -1.38; the general theorem does not assert the same numerical value for other parameters.

## 8. Why the analytic mechanism transfers to every fixed odd zeta value

For odd s≥3, the moment-family notes in this repository use the positive density

    w_s(y)=2y^s/(s-1)! ·(d/dy)^(s-1)[1/(exp(2πy)-1)],  y>0.

Because s-1 is even, its exponential-series representation is positive. It tends to 1/π at zero and satisfies

    w_s(y)=C_s y^s exp(-2πy)[1+O_s(exp(-2πy))]

at infinity, with C_s=2(2π)^(s-1)/(s-1)!. In particular there is a fixed constant A_s such that

    w_s(y)≤A_s(1+y)^s exp(-2πy).

Thus the leading exponential field after y=Ksqrt t is always 2πsqrt t, independently of fixed s. If the same rational damping/root multiplicities and h/K→λ are used, their logarithmic contributions produce the same remaining field terms and the same equilibrium problem. Changing s alters polynomial weight factors, contributing only lower-order Klog K terms under the usual integrable-tail estimate.

This is a statement about the leading real determinant estimate for fixed s and a fixed auxiliary family. It does not automatically transfer the arithmetic coefficient bounds, prime normalizers, degree statements or Bernoulli denominators. Those must be rederived with s present. It also gives no uniform estimate as s grows with K; fixed s is essential to the lower-order assertion.

The conceptual separation is useful: the analytic equilibrium mechanism is universal across the odd-moment family, while the arithmetic determines whether that shared decay is strong enough to force irrationality.

## 9. Relationship to the numerical certificates

The general theorem removes endpoint existence as a numerical obligation. For a specific parameter triple, validated numerics now serve two narrower purposes:

1. locate the unique already-established endpoint pair in a rational box;
2. enclose the resulting scalar energy tightly enough to compare with the independently proved arithmetic exponent.

The existing exact dyadic endpoint certificate and the two independent Arb energy implementations remain useful and reproducible. Their scope is the particular parameter triple in the Fauzan construction, whereas Sections 1–6 above establish the measure's existence and structure for the whole stated family.
