# Independent review of the endpoint certificate

Reviewed artifacts: equilibrium-endpoint-certificate.md and certify_equilibrium_endpoints.py. This review concerns the two endpoint equations and their certified rectangle. It does not certify equilibrium-measure positivity, the variational inequality, or an energy bound.

## Result

No error found in the reviewed analytic proof or interval implementation. Independently reran the standard-library certificate: all four whole-face sign assertions passed. The rectangle contains exactly one root of the stated I=F=0 system by the analytic uniqueness argument, not merely a numerical root approximation.

The rerun returned these outward enclosures (12 decimal places):

    I, a=a_min: [-0.000001673335,-0.000000415767]
    I, a=a_max: [ 0.000000483172, 0.000001740740]
    F, b=b_min: [ 0.000000032042, 0.000000032514]
    F, b=b_max: [-0.000000034810,-0.000000034338]

Every interval excludes zero with the required sign. The parameters on the other coordinate are intervals spanning the complete face.

## 1. Integer interval operations

The represented interval is [lo/2^120,hi/2^120]. Each operation was checked against that representation:

- Rational construction uses floor for the lower endpoint and ceiling for the upper, including negative inputs.
- Addition and subtraction are exact integer endpoint operations.
- Multiplication takes the minimum and maximum of all four corner products, then divides by the scale with outward rounding. This handles mixed signs correctly.
- Division by a negative interval reduces to (-a)/(-b). For a positive denominator interval, all four endpoint quotients are considered and rounded outwards. Intervals containing zero are rejected.
- Square-root endpoints use floor/ceiling integer square roots of endpoint times scale. The perfect-square upper-endpoint test is correct.
- absmax returns an exact dyadic upper bound on absolute value. Decimal output is separately rounded outwards, so displayed signs are not created by rounding.

No binary floating-point constant enters the certificate arithmetic. Divisions such as 1/a and /120 dispatch to interval operations; exact decimal inputs are parsed through Fraction strings. The code has no overflow limit because Python integers are unbounded.

Additional independent exact-fraction checks of the interval operations and representative jets are in check_endpoint_interval_audit.py. These supplement the operation-by-operation proof review; random checks are not treated as a proof of enclosure.

## 2. Taylor jets bound derivatives on entire cells

The jet's zeroth coefficient is the interval covering the full cell [l,r], its first coefficient is exactly one, and its higher coefficients start at zero. For every fixed x0 in that cell, this encloses the Taylor coefficients of x0+h.

The multiplication recurrence is the usual convolution. Reciprocal coefficients satisfy

    c_0=1/a_0,
    c_k=-(sum_{j=1}^k a_j c_(k-j))/a_0.

Square-root coefficients satisfy

    c_0=sqrt(a_0),
    c_k=[a_k-sum_{j=1}^{k-1}c_jc_(k-j)]/(2c_0).

These are exact formal identities for normalized derivatives f^(k)/k!. Interval substitution preserves containment at every step. Radicands and reciprocal denominators remain strictly positive on the specified box; the code also checks the required signs before division.

The a,b parameter intervals have zero derivative coefficients, correctly holding parameters constant when differentiating in the integration variable. Repeated use of the same parameter interval may overestimate dependence, but cannot underestimate it. Thus jet_4 bounds f''''/24 simultaneously on the entire cell and for every parameter value on the face. It is not merely a midpoint derivative estimate.

## 3. Simpson error and parameter uncertainty

For a panel of total width h using its two endpoints and midpoint, Simpson error is at most h^5 sup|f''''|/2880. With normalized jet coefficient4, this is h^5 sup|jet_4|/120, exactly the factor used by the script.

The three Simpson function values are themselves intervals evaluated with full face-parameter uncertainty. The code then adds the symmetric dyadic interval [-err.hi,err.hi], using the upper endpoint of the error estimate. Cell and piece errors are therefore enclosed before signed weights are applied. Negative piece weight -3 reverses interval endpoints through ordinary interval multiplication; it does not invalidate the bound.

Each cell is an independent Simpson panel with its own midpoint, so the number of cells need not be even. The cuts at3/40 and1 are handled as separate smooth integrals with constant weights. Shared endpoints cause no analytic issue, since the jump of q at one point does not affect the integral.

## 4. Infinite-tail transformation

Under u=1/v, the I tail becomes

    1/sqrt[(1+av²)(1+bv²)].

For F, rationalizing the square-root difference cancels the v² Jacobian factor and gives

    (a-b)/{sqrt(1+bv²)[sqrt(1+av²)+sqrt(1+bv²)]}.

This is continuous at v=0. The code implements exactly these formulas and the weight2, with the outside factor1/2 and addition37/40 only for F. The finite-interval rationalization is also exact. Thus no unstable subtraction or singular endpoint is being silently approximated.

## 5. Existence and uniqueness audit

On I=0, subtracting the value at u=c=3/40 from 1/(u²+a) or1/(u²+b) gives an integrand strictly negative on both sides of the sole sign change of q. Therefore J_a,J_b<0 and I_a,I_b>0. Differentiation under the integral is valid because a,b stay bounded away from zero and the derivative tails are integrable.

For a_2>a_1, the ratio h_(a_2,b)/h_(a_1,b) is strictly increasing in u. Subtracting its value at c proves I(a_2,b)>0 whenever I(a_1,b)=0. Thus each horizontal slice has at most one I-root; the certified a-face signs give existence on every slice.

The implicit-function theorem then supplies the unique root curve a(b), with a'(b)=-I_b/I_a<0. The identities F_a=I/4 and F_b=(b-a)J_b/4 on I=0 show that F decreases strictly along that curve. The two certified b-face signs force exactly one F-zero. Since all face signs are strict, the joint root lies in the interior.

The claim is uniqueness inside the specified rectangle. The proof does not assert uniqueness in the entire positive quadrant, and no such claim is needed for this certificate.

## Scope remaining outside this review

Endpoint equations alone do not show the signed balayage candidate is a positive equilibrium measure. Nor do they certify an energy value. Those require independent measure, variational, and energy arguments. This endpoint certificate can be used as a rigorous input to those arguments, but cannot replace them.
