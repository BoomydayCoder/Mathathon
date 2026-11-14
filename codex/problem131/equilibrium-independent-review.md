# Independent audit of the analytic equilibrium replacement

Status: signed-charge algebra, positivity, endpoint calculus, and the global potential inequality independently rederived. The endpoint certificate and both energy certificates were independently rerun and passed. This is a local analytic audit, not validation of every arithmetic or determinant claim in the Fauzan preprint.

## 1. Infinite charge: subtraction is essential

Fix alpha=3/40, lambda=37/40 and 0<a<b. Let q(u) be -3 on (0,alpha), +3 on (alpha,1), and +2 on (1,infinity). The measure q(u)du has infinite mass. Integrals of B_u or of the arcsine measure separately against q are therefore not finite measures; B_u-mu must remain grouped.

Let C_u=sqrt((u²+a)(u²+b)) and dB_u(x)=C_u/(x+u²)dmu(x), with mu unit arcsine on [a,b]. Each B_u has mass one. Uniformly in x on the interval,

    C_u/(x+u²)-1=((a+b)/2-x)/u²+O(u^-4).

Near u=0 the ratio is bounded because a>0. Thus the integral of q(u)(B_u-mu) converges in total variation and has mass zero. The resulting rho=lambda mu+(1/2)integral q(B_u-mu) has mass lambda before positivity is proved.

The same bound permits interchange with logarithmic potentials: the grouped density difference is bounded by an integrable function of u times the arcsine density. Arcsine measures have uniformly bounded absolute-logarithmic integrals on fixed compact sets of evaluation points. This also establishes continuity of the resulting potential.

The field has the convergent normalized representation

    V(x)=V(0)+integral q(u)log(1+x/u²)du,
    V(0)=-2+12alpha(1-log alpha).

The unnormalized infinite integral of q(u)log(x+u²) is not legitimate. Differentiation of the normalized expression is valid for x>0.

## 2. Positivity and the global potential inequality

Write

    I0=integral q/C_u,
    H=lambda+(1/2)integral q(sqrt((u²+a)/(u²+b))-1),
    F(x)=lambda+(1/2)integral q(C_u/(x+u²)-1).

Then H=F(b) and F(a)-F(b)=(b-a)I0/2. If I0=H=0, both endpoint values vanish. For x>0 define J(x)=integral q/[(x+u²)C_u]. Subtract 1/(x+alpha²) times I0. The integrand is strictly negative on both sides of alpha, because q changes sign once in the opposite direction to the decreasing kernel. Hence J(x)<0.

Rational interpolation between the endpoint values gives

    F(x)=(x-a)(x-b)J(x)/2,
    d rho(x)=-sqrt((x-a)(b-x))J(x)dx/(2pi)>0.

The density is bounded, extends by zero to a Hölder-1/2 function, and has finite logarithmic energy.

Important qualification: J(0) diverges negatively. Use the sign and derivative formulas only for x>0; the potential inequality at zero follows by continuity.

The arcsine Stieltjes transform and a partial-fraction calculation give

    G_mu(z)=1/sqrt((z-a)(z-b)),
    G_Bu(z)=(1+C_u G_mu(z))/(z+u²),

with the square-root branch asymptotic to z. Inside [a,b], the real boundary values give 2U_rho'-V'=0. Outside it,

    (2U_rho-V)'=sqrt((x-a)(x-b))J(x).

This derivative is positive on (0,a) and negative on (b,infinity). Thus the potential is constant on the support and no larger anywhere else on [0,infinity). Continuity covers 0,a,b.

This proves the needed comparison inequality. Uniqueness of the equilibrium measure additionally uses strict zero-mass energy negativity; it is not needed for the upper bound. Strictness follows, for example, from the Gaussian-energy representation and injectivity of Gaussian convolution.

## 3. Endpoint derivatives and executable certificate

At a zero of I0,

    I0_a=-J(a)/2>0, I0_b=-J(b)/2>0,
    a'(b)=-J(b)/J(a)<0.

Furthermore H_a=I0/4 and H_b=-I0/4+(b-a)J(b)/4. Along the zero branch, dH/db=(b-a)J(b)/4<0. Differentiation is justified uniformly on compact parameter boxes with a>0.

The statement I0(b,b)=V'(b)>0 needs a restriction on b: V' is negative near zero. The proposed box lies in its positive region. For a global statement, the one-sign-change kernel implies every zero of V' is an upward crossing; its small/large-x signs then give a unique minimum of V.

The independent rerun of certify_equilibrium_endpoints.py passed all entire-face signs:

    I at a_min, all b: [-0.000001673335,-0.000000415767],
    I at a_max, all b: [ 0.000000483172, 0.000001740740],
    H at b_min, all a: [ 0.000000032042, 0.000000032514],
    H at b_max, all a: [-0.000000034810,-0.000000034338],

on a in [0.00007119508,0.00007119517] and b in [0.75394957,0.75394967]. Opposite a-face signs give a root for every b; positive derivatives at all roots imply uniqueness and a continuous branch. Opposite H signs and strict branch monotonicity then give a unique pair in the box.

The custom interval arithmetic was inspected. Operations round outward; positive-denominator division checks all corner quotients; isqrt encloses square roots. Interval Taylor coefficient 4 bounds f''''/4! throughout each subinterval. The Simpson error h^5*abs(jet4)/120 equals h^5*sup|f''''|/2880. The reciprocal tail substitutions and signed weights are correct.

## 4. A quadrature-free plateau constant

Let cap=(b-a)/4 and K(u)=log[4u²/(sqrt(u²+a)+sqrt(u²+b))²]. Explicit balayage potentials give

    M=2lambda log(cap)-V(0)+integral qK.

Here K=O(u^-2) at infinity and K=2log u+O(1) near zero. Also K'=2/u-2u/C_u. Piecewise integration by parts gives

    integral qK=K(1)-6alpha K(alpha)-2integral q(1-u²/C_u).

The endpoint products uK vanish at zero and infinity. The endpoint equations give integral q(1-u²/C_u)=2lambda, with the convergent parentheses kept grouped. Thus

    M=2lambda log(cap)-V(0)+K(1)-6alpha K(alpha)-4lambda.   (M)

This formula was independently confirmed by the integral agent and is used in the validated energy scripts. It removes the plateau quadrature entirely.

## 5. One smooth charge integral for the energy

Symmetrizing the equilibrium Hilbert-transform equation gives the virial identity integral xV'(x)d rho=lambda². A rigorous derivation cuts away |x-y|<epsilon and uses x/(x-y)+y/(y-x)=1. The diagonal-strip mass vanishes in the limit. Hölder continuity of the density justifies convergence of the truncated Hilbert transforms.

Direct calculation gives

    V(x)-2xV'(x)=log(x+1)-6alpha log(x+alpha²)-2+12alpha,

and hence

    integral Vd rho=2lambda²+U_rho(-1)-6alpha U_rho(-alpha²)
                         +lambda(-2+12alpha).             (E)

For s>0 define

    D_s=(sqrt(s+a)+sqrt(s+b))²/4,
    tau_s=(sqrt(s+b)-sqrt(s+a))/(sqrt(s+b)+sqrt(s+a)).

The angular arcsine coordinate makes B_u a Poisson kernel of parameter -tau_(u²). Applying its Fourier moments to log(s+x) gives

    U_Bu(-s)-U_mu(-s)=2log(1-tau_s tau_(u²)),
    U_rho(-s)=lambda log D_s+integral q(u)log(1-tau_s tau_(u²))du.   (P)

Combining the two instances of (P) in (E) gives one smooth charge integral. Its tail is O(u^-2); its integrand is bounded near zero since tau_s,tau_0<1. Under u=1/v,

    tau_(1/v²)=(b-a)v²/(sqrt(1+av²)+sqrt(1+bv²))²,

so the removable endpoint is explicit. These identities were independently confirmed by the integral agent. A nonrigorous consistency calculation gave integral Vd rho approximately 1.91586900821353; the separate validated script now certifies an enclosure.

## 6. Carlson formula and analytic callbacks

In J_T(x)=integral_0^T du/[(u²+x)sqrt((u²+a)(u²+b))], substitute u=sqrt(a)tan(theta). With s=sin(theta_T), c²=1-s² and k²=1-a/b, the result is

    [s RF(c²,1-k²s²,1)
       -(a/x)s³ RJ(c²,1-k²s²,1,c²+(a/x)s²)/3]/(x sqrt(b)).

This exactly matches certify_equilibrium_energy.py. The algebraic decomposition is cos²theta/[x cos²theta+a sin²theta]=(1/x)[1-(a/x)sin²theta/(1-(1-a/x)sin²theta)]. Combining 2J_infinity+J_1-6J_alpha produces q correctly.

The change x=a+(b-a)sin²theta gives the measure factor -(b-a)²sin²theta cos²theta J(x)/pi, also matching the script. The elementary log/arctangent expression for V is correct.

Accepted analytic boxes satisfy Re(x)>0. This keeps the principal square root, log(x+u²), and u/sqrt(x) away from the arctangent cuts. RF's constant arguments are nonnegative real with the nonzero ones positive. RJ's fourth argument c²+(a/x)s² has positive real part. These checks suffice for analyticity on every box accepted by validated quadrature. NaN on a rejected box requests subdivision.

The endpoint rectangle need not consist entirely of equilibrium pairs: the formulas evaluated over the box enclose their values at the unique actual root, which is enough. The real projection of a complex enclosure contains the true real integral; the imaginary interval is checked to contain zero.

## 7. Independent executable reruns and result

Using PYTHONPATH=/workspace/scratch/odd-zeta-python, both energy scripts were independently rerun. JSON writing was suppressed for this read-only audit; every mathematical assertion remained active.

The Carlson certificate passed with

    M in [-6.660781 +/-2.86e-7],
    integral Vd rho in [1.916 +/-3.05e-4],
    objective in [-1.3855 +/-9.64e-5] < -1.38.

The separate elementary potential/virial certificate passed with

    integral Vd rho in [1.9159 +/-7.03e-5],
    objective in [-1.3855 +/-2.95e-5] < -1.38.

The objective is (lambda M-integral Vd rho)/2+C*, obtained from 2I(rho)-integral Vd rho=lambda M. The two analytic representations agree and both yield the required strict bound.

The latest elementary script was then rerun independently in its standalone form, with no import of the Carlson script and no elliptic-function calls; it passed with the same enclosure. Its charge weights are 3 times the integral on [0,1], minus 6 times the integral on [0,alpha], plus 2 times the reciprocal tail, exactly reproducing q. In that tail, log(1-z)/v² is evaluated as -tau_s*z0*2F1(1,1;2;z), with z=tau_s*z0*v². The identity 2F1(1,1;2;z)=-log(1-z)/z makes the value at v=0 removable and correct. The analytic callback requires Re(z)<1, excluding the hypergeometric branch cut [1,infinity); the square-root and logarithm calls use their analytic-domain checks. The virial constants and the two off-support potential coefficients match (E) and (P) above.

## 8. Scope

This construction gives a positive bounded-density measure of mass lambda, a global potential inequality, and a validated scalar objective below -1.38. Bounded density yields the circle-regularization estimate integral log(max(|x-t|,epsilon)/|x-t|)d rho(t)<=2||density||_infinity epsilon, so the source's asymptotic regularization argument applies with O(epsilon) error.

A completely explicit bound is available without further quadrature. Since |q(u)|<=3 and x>=a on the support,

    |J(x)|<=3 integral_0^infinity du/(a+u²)²=3pi/(4a^(3/2)),
    ||density||_infinity<=3(b-a)/(16a^(3/2)).

The endpoint enclosure makes this a finite certified constant. It is deliberately coarse, but sufficient for the asymptotic regularization error; optimizing it is unnecessary for the limiting exponent.

It does not independently establish the preprint's arithmetic denominator bounds, determinant integrality, moment identities, or claimed irrationality theorem. Those are separate proof obligations. This audit supports replacement of the local analytic sixteen-arcsine certificate only.
