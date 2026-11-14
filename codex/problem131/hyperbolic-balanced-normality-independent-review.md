# Independent full audit of all-balanced hyperbolic normality

Status: the full proof in `hyperbolic-balanced-normality.md` passed this
independent audit. I found no remaining analytic, zero-count, or normality gap.
The conclusion is all balanced indices (n,…,n), not arbitrary unbalanced
indices, denominator estimates, approximation rates, or irrationality.
Historical priority is unresolved. The proof originated with the root agent
and was written in detail by the integral agent; this is a separate audit.

## 1. Global continuation and all endpoint constants

The real Stieltjes derivation is correct. With
ρ(v)=arcosh(1/√v), one has ρ′(v)=−1/[2v√(1−v)] and ∫₀¹ρ(v)dv=1.
Integration by parts against log(1+v/u) has zero endpoint terms: near0 it is
O(v|log v|), and near1 the density tends to0. Differentiating the resulting
integral gives

    I′(u)=−arsinh(1/√u)/(u√(1+u)),

exactly the derivative of arsinh(1/√u)². Both quantities tend to0 at infinity.
Thus

    F(z)=π^−2 ∫₀¹ ρ(v)/(z+v)dv

is a single-valued analytic function on C\[−1,0]. This integral, rather than a
global unchecked principal logarithm, specifies the continuation. Its expansion
F(z)=1/(π²z)+O(z⁻²) is uniform outside a large circle.

The explicit logarithm is legitimate locally by analytic continuation from
the positive half-line. At z=−v+i0, 0<v<1,

    log((1+√(1+z))/√z)=ρ(v)−iπ/2.

The lower limit is its conjugate. Therefore the upper/lower values are

    F±(−v)=G(v)−1/4 ∓ i√G(v),
    G(v)=ρ(v)²/π²=F(v/(1−v)).

All π factors and the Stieltjes sign are correct. In particular the upper
boundary has negative imaginary part, and the jump of F is −2i√G.

The endpoint bounds used in the contour proof are valid uniformly in a local
slit neighborhood. Near0, use the local branch of log z with bounded argument:
F(z)=O(1+|log|z||²). Near−1, the continued logarithm has finite value
±iπ/2 and a square-root correction, so F is bounded. Consequently the small
endpoint circles in R/W contribute O(ε|log ε|²⁽ⁿ⁻¹⁾) and O(ε). No pole or
point-mass endpoint contribution survives. The artificial negative-real-axis
branch of the principal expression cancels after squaring to the left of−1;
the Stieltjes representation already proves this globally.

On the actual cut, G=O(log²v) near0 and G=O(1−v) near1. The density in the
Cauchy formula is therefore O(1+|log v|²ⁿ⁻³) and O(√(1−v)), respectively.
These bounds are integrable, including after multiplication by any polynomial
in v. They justify endpoint limits and expansion at infinity.

## 2. Jump triangularity and a short Chebyshev description

The jump polynomial has an independent useful description. Since

    w−1/4+i√w=(√w+i/2)²,

one gets

    J_k(w)=(w+1/4)^(k−1)
           U_(k−1)((w−1/4)/(w+1/4)),

where U is the ordinary Chebyshev polynomial of the second kind. Equivalently,

    J_0=0, J_1=1,
    J_(k+1)=2(w−1/4)J_k−(w+1/4)²J_(k−1).

This proves degree k−1 and leading coefficient k directly by induction. It also
shows that J_k has k−1 distinct positive roots

    w_l=(1/4)cot²(lπ/(2k)),  l=1,…,k−1.

Those additional root facts are not needed for the AT proof, but explain the
otherwise fairly long binomial jump expansion. No new historical claim is
attached to this elementary Chebyshev identity.

The formal degree-drop argument is exact: if k is the highest nonzero P_k
with k≥1, then J(v,w) has leading coefficient kP_k(−v) in degree k−1.
Thus the coefficient tuple for the lower family cannot disappear.

The Möbius coefficient map is even stronger than merely invertible. For
polynomials of degree<m, put

    T_m P(u)=(1+u)^(m−1) P(−u/(1+u)).

It is an integer involution: T_m²=I. In the monomial basis its matrix is
triangular with diagonal (−1)^j, hence determinant (−1)^(m(m−1)/2).
Thus multiplying by (1+u)^(m−1) preserves the degree bound and a nonzero
coefficient tuple. The multiplier is positive on u>0. The variable change
v=u/(1+u) is an analytic diffeomorphism and preserves zero multiplicities.

## 3. Cauchy representation and the exact moment count

Assume n≥2 and the lower-n theorem. If R has N=mn positive zeros counted
with multiplicity, choose their monic product W, repetitions included.
Division introduces only removable singularities at the chosen positive
zeros; every other point off the cut was already analytic. Since F=O(1/z)
and deg P_k≤m−1,

    R/W=O(z^(m−1−mn))=O(z^−M−1),  M=m(n−1).

W(−v) is nonzero and has the single fixed sign (−1)^N throughout [0,1].
The audited contour and endpoint bounds give exactly

    R(z)/W(z)=(1/π)∫₀¹ √G(v) J(v,G(v))/[W(−v)(z+v)]dv.

The factor1/π follows from upper-minus-lower jump −2i√G J/W. There is no
missing factor2, π, or residual polynomial, and no other pole is crossed.
For |z|>2 the geometric expansion of1/(z+v) is uniform in v∈[0,1] and is
dominated by an integrable density. Comparing with O(z⁻ᴹ⁻¹) kills exactly
the coefficients z⁻¹,…,z⁻ᴹ. Hence there are exactly M required zero moments,
with test powers v⁰,…,vᴹ⁻¹. The proof neither loses one moment nor assumes an
extra asymptotic cancellation.

If all higher P_k vanish, the original polynomial has too few zeros. Otherwise
the jump/Möbius argument puts J(v,G(v)) in the lower rectangular family of
dimension M with a nonzero coefficient tuple. By induction it is nonzero
and has at most M−1 zeros, counted with multiplicity, hence at most M−1 sign
changes. Multiplying by the polynomial with simple roots at its sign changes
turns it into a nonzero function of one sign. Its degree is≤M−1, and the
weight √G/|W(−v)| is strictly positive. Its integral therefore cannot be zero,
contradicting precisely those M moments.

The same induction rules out an identically zero R with a nonzero coefficient
tuple: its zero jump would contradict independence of the lower family, unless
all higher P_k vanish, after which P₀ must also vanish. This can be stated
before selecting W if one prefers to avoid speaking of zero multiplicities
for an identically zero function. There is no circularity in the draft's
handling of this case.

This proves the full multiplicity-counting bound, not merely a count of
sign-changing or distinct zeros.

## 4. The normality and root lemma needs no complete-prefix assumption

The change u=csch²(πy), t=y²=F(u) is an analytic diffeomorphism on the positive
axis with nonzero derivative. Thus the exact orthogonality function space is

    span{t^k csch(π√t)^(2j−2): 1≤j≤m, 0≤k<n},

which becomes the proved rectangle span{u^j F(u)^k: j<m,k<n}. Its dimension
is N, and every nonzero member has at most N−1 zeros counted with multiplicity.

All products with a polynomial are integrable against the common positive
measure ω₁(y)dy. At0, the most singular factor comes from component j=m and
ω_m(y)=O(1); extra powers y²ᵏ only improve this. At infinity each component
has exponential decay times a fixed polynomial. No hidden Jacobian is needed
because the integrals remain written in y.

The escaping-extra-roots lemma in §5 is sound. To audit its limiting step:

1. Prescribe the r sign-change points of a polynomial v and N−1−r distinct
   extra roots tending to infinity. N−1 equations in N coefficients have a
   nonzero solution q_R. Linear independence means it is a nonzero function.
2. The full-space multiplicity bound forces all these N−1 roots to be simple
   and forbids any additional root. Thus, after choosing its overall sign,
   vq_R≥0 on every fixed compact positive interval for all sufficiently large R.
   Even-multiplicity roots of v cause no problem.
3. Normalize coefficient vectors, not function values, to Euclidean norm1.
   Compactness yields a convergent subsequence with norm-one limit. Independence
   makes the limiting function q nonzero, and local convergence gives vq≥0.
4. Since q is a fixed linear combination in the same space, its product with v
   is integrable and orthogonal **by linearity**. One never needs a global
   dominated-convergence assertion for the escaping-root sequence. The analytic
   nonzero product is positive on an interval, contradicting orthogonality.

This excludes an orthogonal nonzero polynomial of degree<N, so the N×N mixed
moment matrix is nonsingular. It also excludes fewer than N sign changes for
the resulting monic degree-N solution; hence all N roots are simple and lie
in (0,∞). No complete-prefix or arbitrary-unbalanced AT theorem was used.

For balanced indices every component has the same k range. The previously
proved triangular change to t^(m−j)μ_(2j+1) therefore preserves all conditions
and gives a nonsingular rational moment system. Its monic solution is rational.
The irrational scale factors in the original component weights do not enter
this rationality conclusion.

## 5. Independent checks and remaining scope

`check_hyperbolic_balanced_independent.py` uses the Chebyshev recurrence rather
than the main checker's binomial formula. It verifies the jump identities,
positive simple J roots through k=11, and the integer Möbius involution through
m=12. Independent60-digit quadrature verifies the Stieltjes representation
at three positive and one nonreal point. Both sides of the cut and the
analytic continuation to the left of−1 are checked separately. These checks
corroborate, rather than replace, the all-index proof audited above.

The all-balanced normality theorem is established by this proof. General
explicit coefficients, their denominator growth, suitable approximation
errors, and any irrationality consequence remain open. An unequal-degree AT
extension has not been inferred from the rectangular theorem.
