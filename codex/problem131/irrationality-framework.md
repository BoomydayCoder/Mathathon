# Irrationality recipes and unexpected connections

Status: the elementary lemmas, integral arguments, and Mellin/Legendre identities below are independently derived and proved here. Classical theorems are explicitly identified as inputs. This is not a novelty claim. Work is aligned with the user's revised research program; the exact official Mathathon statement is archived in [problem131-official.md](../sources/problem131-official.md).

## 1. The arithmetic target comes before the auxiliary function

**Integer squeeze criterion.** Suppose alpha is real and for arbitrarily large n we have rational A_n,B_n and a positive integer d_n with

    d_n A_n, d_n B_n in Z,
    0 < |d_n(A_n alpha-B_n)| -> 0.

Then alpha is irrational. Indeed alpha=p/q would make q d_n(A_n alpha-B_n) a nonzero integer of absolute value tending to zero.

The denominator A_n is allowed to be enormous. The relevant quantity is the *linear form* A_n alpha-B_n after clearing arithmetic denominators, not merely the rational-approximation error alpha-B_n/A_n. Numerical accuracy alone proves nothing about irrationality.

If d_n <= exp((delta+o(1))n) and |A_n alpha-B_n| <= exp((-rho+o(1))n), the useful inequality is rho>delta. Equality is inconclusive and requires polynomial factors or sharper estimates. Nonvanishing is an independent obligation.

For zeta(3), Zudilin provides d_n=D_n^3 with D_n=lcm(1,...,n). With D_n<3^n and exponential error (sqrt(2)-1)^(4n), the decisive product is

    27 (sqrt(2)-1)^4 = 27(17-12 sqrt(2)) < 1.

With a merely D_n^4 denominator bound the corresponding elementary estimate is 81(17-12 sqrt(2))>1, so this argument fails. Improving an arithmetic exponent by one can matter more than making the real-variable estimate look cleaner.

## 2. A complete elementary proof for exp(r), r nonzero rational

It suffices to consider r=a/b>0, a,b positive integers: rationality of exp(r) is equivalent to rationality of its reciprocal exp(-r).

Choose

    f_n(x)=x^n(a-bx)^n/n!,    0<=x<=r,
    F_n(x)=sum_{j=0}^{2n} (-1)^j f_n^(j)(x).

Three design decisions explain the construction.

1. Zeros of order n at both endpoints make all low-order boundary derivatives vanish.
2. Division by n! makes the interior exponentially small relative to a factorial, while the surviving derivatives remain integers.
3. Since d(exp x)/dx=exp x, a finite alternating derivative sum provides an exact antiderivative: (exp(x) F_n(x))'=exp(x) f_n(x).

**Endpoint integrality.** Expanding at zero gives

    f_n(x)=sum_{k=0}^n (-b)^k binom(n,k) a^(n-k) x^(n+k)/n!.

Thus each derivative at zero is zero or an integer multiple of (n+k)!/n!. At the other endpoint write x=r-y; then

    x^n(a-bx)^n=(a-by)^n y^n.

The same expansion shows all derivatives at r are integers too. Hence F_n(0),F_n(r) are integers. This reflection calculation is essential: rational endpoints alone do not imply integer derivatives.

**Positivity and smallness.** The integral

    I_n=integral_0^r exp(x) f_n(x) dx
       =exp(r) F_n(r)-F_n(0)

is strictly positive. Since x(a-bx)<=a^2/(4b),

    0<I_n<=r exp(r) [a^2/(4b)]^n/n! ->0.

If exp(r)=p/q, then q I_n is a positive integer tending to zero, contradiction. This proves exp(r) irrational for every nonzero rational r without any transcendence theorem.

**Why the standard proof for e does not automatically handle all r.** Multiplying the Taylor tail of exp(a/b) by n! b^n produces terms involving a^(n+j), and the simple estimate used for e need not tend to zero. Symmetric endpoint vanishing repairs the construction. This is a useful example of an auxiliary polynomial solving an arithmetic defect, not just accelerating convergence.

## 3. Sine equations, pi², and nonzero rational powers pi^r

The newer [direct pi² proof](pi-powers.md) uses a finite inverse on the fixed interval(0,1), with an independently checked integer boundary identity. The classical pi argument below motivates the same mechanism.

Suppose pi=a/b with positive integers a,b and use the same f_n on [0,pi]. Let

    G_n=f_n-f_n''+f_n^(4)-...+(-1)^n f_n^(2n).

Then G_n''+G_n=f_n and

    d/dx [G_n'(x) sin x-G_n(x) cos x]=f_n(x) sin x.

Consequently

    J_n=integral_0^pi f_n(x) sin x dx=G_n(pi)+G_n(0).

The endpoint arithmetic from section 2 makes J_n an integer. The integrand is positive in the interior, and

    0<J_n<=pi [a^2/(4b)]^n/n! ->0,

which is impossible. Thus pi is irrational.

This is the Niven-type mechanism. It has the same arithmetic polynomial as the exponential proof; only the differential equation and boundary identity change. In the exponential case the endpoint contributes the unknown exp(r); in the sine case the hypothesized rational endpoint is itself the unknown pi, while sin(pi)=0 and cos(pi)=-1 supply integer boundary data.

**Powers of pi require a separate argument.** Irrationality of pi does not imply irrationality of pi^2 (sqrt(2) is a counterexample to that inference for general numbers). A uniform proof for pi^r, every positive integer r, uses the classical transcendence of pi. If pi^r were algebraic, pi would be algebraic over the algebraic numbers because it solves X^r-pi^r=0, and transitivity of algebraic extensions would make pi algebraic over Q, contradiction. Therefore every positive integral power of pi is transcendental and in particular irrational. This is a consequence of Lindemann's theorem, not a direct elementary proof for each power. The same conclusion holds for every nonzero rational exponent r=p/q: algebraicity of pi^(p/q) would imply algebraicity of pi^p, and then of pi. The explicitly reconstructed integer-power irrationality proof in the companion note also gives irrationality for all such rational powers by raising to q and, if needed, taking a reciprocal. No claim for arbitrary real exponents is intended.

A still shorter advanced route to the section-2 result is Lindemann's theorem that exp(beta) is transcendental for nonzero algebraic beta. Its short application conceals substantial machinery. Neither proof length nor prerequisite burden should be silently traded for the other.

## 4. Unexpected link I: the exponential integral is a Padé approximation

Treat z as a variable and define

    P_n(z)=sum_{j=0}^n [(2n-j)!/(j!(n-j)!)] z^j,
    Q_n(z)=P_n(-z).

Repeated integration by parts yields the exact identity

    Q_n(z) exp(z)-P_n(z)
      =(-1)^n/n! integral_0^z exp(t) t^n(z-t)^n dt.        (1)

Proof: set f(t)=t^n(z-t)^n/n! and F=sum(-1)^j f^(j). The endpoint expansions give F(0)=(-1)^n P_n(z), F(z)=(-1)^n Q_n(z), and integrate (exp(t)F(t))'=exp(t)f(t).

On putting t=zu, the right side is

    (-1)^n z^(2n+1)/n! integral_0^1 exp(zu) u^n(1-u)^n du.

Thus Q_n exp(z)-P_n has a zero of order 2n+1 at zero, though P_n and Q_n each have degree n. Since Q_n(0)=(2n)!/n! is nonzero, P_n/Q_n is the diagonal [n/n] Padé approximant to exp(z).

The same identity simultaneously explains:

- approximation: 2n+1 matched Taylor coefficients;
- arithmetic: factorial ratios in P_n,Q_n;
- positivity/nonvanishing: the integral at positive real z;
- smallness: the factorial and maximum of u(1-u).

This is useful for discovery. Instead of guessing an integral, impose high-order contact Q exp-P=O(z^(2n+1)), solve a linear system, and search for a positive remainder representation. Conversely, endpoint vanishing constructs the Padé approximant without solving that system.

## 5. Unexpected link II: Zudilin's rational function is a Legendre Mellin transform

Let

    L_n(x)=1/n! (d/dx)^n [x^n(1-x)^n]
          =sum_{k=0}^n (-1)^k binom(n,k)binom(n+k,k) x^k.

This is a shifted Legendre polynomial with convention L_n(0)=1. For every polynomial q of degree<n,

    integral_0^1 L_n(x) q(x) dx=0.                       (2)

Proof: integrate by parts n times; endpoint derivatives of x^n(1-x)^n of orders below n vanish, and q^(n)=0. Thus orthogonality is a direct consequence of endpoint zeros.

Now take a Mellin transform. For t>n, repeated integration by parts and the beta integral give

    integral_0^1 x^(t-1)L_n(x) dx
      =(-1)^n [(t-1)...(t-n)]/[t(t+1)...(t+n)].          (3)

For completeness, the beta integral used here is elementary:

    integral_0^1 x^(t-1)(1-x)^n dx
      =n!/[t(t+1)...(t+n)],

obtained by n integrations by parts for t>0. Multiplying it by the falling factorial from differentiating x^(t-1) proves (3) initially for t>n. Expanding L_n makes the left side a rational function sum l_k/(t+k), so equality extends as a rational identity and in particular to all t>0.

**The discovery explanation.** Orthogonality (2) forces the Mellin transform to vanish at t=1,...,n. A polynomial of degree n has a Mellin transform with poles only at 0,-1,...,-n. These prescribed zeros and poles essentially force the rational product in (3). Zudilin's auxiliary function is therefore not an arbitrary product: it is the square of this Mellin transform.

Writing l_k=(-1)^k binom(n,k)binom(n+k,k), we obtain

    R_n(t)=[integral_0^1 x^(t-1)L_n(x)dx]^2,
    A_n=sum_{k=0}^n l_k^2
       =sum_{k=0}^n binom(n,k)^2 binom(n+k,k)^2.

A_n is the squared Euclidean norm of the coefficient vector of L_n. Equivalently it is the constant coefficient of L_n(z)L_n(z^-1). Do not confuse this with the function-space norm integral_0^1 L_n(x)^2 dx, which is a different quantity.

Differentiation in Mellin space inserts a logarithm. Consequently the normalization F_n=-sum_{t>=1}R_n'(t) from Zudilin satisfies the exact identity

    F_n=integral_0^1 integral_0^1
          L_n(x)L_n(y) [-log(xy)]/(1-xy) dx dy.          (4)

To justify the interchanges, L_n is bounded on [0,1], and the nonnegative majorant [-log(xy)]/(1-xy) is integrable. Expanding the geometric series evaluates its integral as 2 sum_{t>=1}t^-3=2 zeta(3). This also verifies F_0=2 zeta(3).

This supplies an exact bridge between the rational-function and integral proofs. The companion beukers-explained.md develops the subsequent integrations by parts and positive triple integral. Combining that analytic argument with Zudilin's partial-fraction arithmetic can eliminate the two recurrence certificates from a proof of irrationality. One must retain the full transformation and bound, rather than declare positivity of (4): L_n changes sign and positivity in (4) is not manifest.

## 5a. A complete companion recipe for zeta(2)

The same shifted Legendre polynomial gives a particularly transparent double-integral proof. Set

    K_n=integral_0^1 integral_0^1 L_n(x)(1-y)^n/(1-xy) dx dy.

Expanding the two integer-coefficient polynomials reduces the arithmetic to

    J_kl=integral integral x^k y^l/(1-xy) dx dy
        =sum_(t>=1) 1/[(t+k)(t+l)].

If k=l this is zeta(2)-H_k^(2). If k!=l it is (H_k-H_l)/(k-l). Since 0<=k,l<=n, D_n^2 clears the rational parts in both cases. Therefore

    K_n=a_n zeta(2)-b_n,   a_n in Z, D_n^2 b_n in Z.

Integrate by parts n times in x, using the Rodrigues definition of L_n. All endpoint terms vanish and

    (-1)^n K_n=integral integral
       [x(1-x)y(1-y)]^n/(1-xy)^(n+1) dx dy >0.           (6)

Indeed the nth x-derivative of 1/(1-xy) is n! y^n/(1-xy)^(n+1); the n! cancels the Rodrigues factor. Integration can first be done for 0<=y<1, and the resulting nonnegative integrand and explicit bound below justify passage to the full square.

To maximize h=x(1-x)y(1-y)/(1-xy), put s=sqrt(xy). Since x+y>=2s,

    h=xy(1-x-y+xy)/(1-xy)<=s^2(1-s)/(1+s).

Equality occurs when x=y=s. The derivative of the last expression is

    2s(1-s-s^2)/(1+s)^2,

so the maximum occurs at g=(sqrt(5)-1)/2. Since g^2=1-g and 1+g=1/g, the maximum is g^5=(5sqrt(5)-11)/2. Equation (6) gives

    0<|K_n|<=zeta(2) g^(5n).

The elementary lcm appendix in [zudilin-explained.md](zudilin-explained.md), section 10, proves without the prime number theorem that

    log D_n <= c n + O(log(n+1)^2),
    c=(14 log2+9 log3+5 log5)/25 ≈1.105550427521.

This closes the denominator dependency: the exponential rate is 2c+5 log g≈-0.19495827<0. An exact sign check avoids relying on decimals. Since 1/g=(1+sqrt(5))/2>8/5, it suffices to check

    2^28 3^18 5^10 < (8/5)^125,

which is equivalent to the true integer inequality 3^18 5^135 < 2^347. Thus D_n^2 |K_n| tends to zero and the integer squeeze proves zeta(2) irrational. The sharper classical D_n<3^n would also suffice, but is no longer needed as an unproved input.

The connection is precise: the golden ratio appears as a stationary point of a two-variable kernel, whereas sqrt(2) appears in the zeta(3) kernel. These constants encode the competition between endpoint zeros and the singular denominator. They are not decorative numerical coincidences.

Historical attribution: the construction is the classical Beukers mechanism; a primary citation is F. Beukers, A note on the irrationality of zeta(2) and zeta(3), Bull. London Math. Soc. 11 (1979), 268–272, as cited in Zudilin. This section is an independent derivation, not a quotation from an independently retrieved copy of Beukers's paper.

## 6. Unexpected link III: recurrences are discrete differential equations

For solutions A,B of a_n u_(n+1)=b_n u_n-c_n u_(n-1), the Casoratian

    W_n=A_n B_(n+1)-A_(n+1)B_n

obeys W_n=(c_n/a_n)W_(n-1). This is the discrete analogue of Abel's Wronskian identity. It turns a recurrence into an exact formula for successive approximation errors:

    B_(n+1)/A_(n+1)-B_n/A_n=W_n/[A_n A_(n+1)].

Summing these increments yields the limiting ratio and its error. For Apéry's normalization A_0=1,A_1=5,B_0=0,B_1=6,

    W_n=6/(n+1)^3.

This explains why one recurrence can simultaneously describe a large integer sequence and a tiny linear form. The arithmetic denominators and identification of the limit as zeta(3) still require independent arguments. The existing recurrence.md gives details.

The generating function adds another perspective. If A_n obeys Apéry's recurrence and A(z)=sum A_n z^n, with theta=z d/dz, coefficient matching gives

    [theta^3-z(2theta+1)(17theta^2+17theta+5)
             +z^2(theta+1)^3] A(z)=0.                  (5)

For n>=2, the z^n coefficient is exactly the recurrence at n-1; the n=1 coefficient says A_1=5A_0. This is a third-order differential equation despite the second-order recurrence. The shift order and derivative order count different things.

The coefficient of the highest derivative is proportional to

    1-34z+z^2,

whose roots are 17+-12 sqrt(2)=(sqrt(2)+-1)^4. These are the same numbers controlling exponential growth and decay in Apéry's recurrence. This coincidence is forced by the leading coefficients: insert u_n approximately lambda^n into the recurrence and obtain lambda^2-34lambda+1=0. A rigorous asymptotic theorem is needed to promote this heuristic to an asymptotic formula, but the algebraic relation is exact.

This offers a test for proposed simplifications: a representation with a different claimed exponential rate must reconcile that rate with the recurrence and its dominant singularities. Algebraic geometry or modular parametrizations can explain more structure, but no specific such identification is claimed or proved here.

## 7. Discrete integration by parts and why creative telescoping helps

The central certificate in Zudilin has the form

    linear combination of R_(n-1),R_n,R_(n+1)
      =S_n(t+1)-S_n(t).

Summing over t leaves boundary terms, just as integrating a derivative does. This is the discrete counterpart of the repeated integration by parts that made the exponential and Legendre constructions work. Endpoint zeros are doing the same job in both settings: they kill boundary obstructions.

A rational summand has rational ratios under n-shifts and t-shifts. After dividing a proposed certificate by the summand, the problem becomes a rational identity. Clearing denominators converts it to polynomial coefficient equations. Thus an identity that looks mysteriously tailored can often be *discovered* by finite-dimensional linear algebra and then proved by a short exact certificate. It is important to show both the search ansatz and the final verification, and to check boundaries separately.

Primary methodological reference retrieved: Petkovšek, Wilf, Zeilberger, A=B (1996), Chapter 6, https://www.math.rutgers.edu/~zeilberg/AeqB.pdf . Local PDF and text are in codex/sources/. The book's introductory worked discussion explicitly distinguishes finding a summand recurrence, summing its telescoping side, and using initial values to identify the sum.

## 8. A practical research recipe

Start with a constant's series or integral representation. Then:

1. Specify the exact rational linear form and the denominator budget needed for the integer squeeze.
2. Force endpoint vanishing, moment cancellation, or prescribed poles to meet that budget.
3. Seek a second representation where sign and exponential decay are visible.
4. Connect the representations by an explicit transform, integration by parts, or a common recurrence with initial values.
5. Compare decay with denominator growth before polishing the proof.
6. When generalizing, identify each new unwanted constant and prove its cancellation; do not assume symmetry removes it.

For higher odd zeta values, the primary obstacle is often that a construction yields a small linear form in several odd zeta values. Such a form can prove at least one is irrational but does not isolate zeta(5). A successful isolated-value proof needs either exact cancellation of every other coefficient or another rigorous elimination mechanism that preserves arithmetic and decay. Smallness of a many-constant linear form is not by itself progress to a specific value's irrationality.

## References and verification boundaries

- W. Zudilin, An elementary proof of Apéry's theorem, arXiv:math/0202159 (2002), primary paper retrieved and read through the local extraction. Its normalization is F_n=2 A_n zeta(3)-2 B_n relative to the common A_0=1 normalization.
- M. Petkovšek, H. S. Wilf, D. Zeilberger, A=B (1996), primary authors' PDF retrieved. See above.
- I. Niven, A simple proof that pi is irrational, Bull. Amer. Math. Soc. 53 (1947), 509, DOI 10.1090/S0002-9904-1947-08821-2: bibliographic reference, original not successfully retrieved in this session (AMS returned 403; attempted Project Euclid PDF URL returned HTML). The proof in section 3 is fully written and checked independently, not quoted from an unseen PDF.
- Lindemann's transcendence theorem is invoked as an established theorem, not proved or independently audited here. No claim is made that the elementary integral arguments prove transcendence.
- Primary-source provenance of the exact hybrid Mellin presentation has not been surveyed. Its ingredients are classical; independent derivation is not historical novelty.

## 9. What changes between irrationality and transcendence, and why odd zeta is harder

The companion [pi-powers.md](pi-powers.md) now gives explicit classical proofs for every positive power of pi and for transcendence of exp(r), r nonzero rational; both were independently audited in pi-powers-review.md. These expand the underlying machinery rather than invoke transcendence as a black box.

The change in the auxiliary polynomial is revealing. For irrationality of exp(r), two endpoints and positivity suffice. A hypothetical polynomial relation among powers of exp(r) introduces finitely many endpoints jr. Give each nonzero endpoint a zero of order p, give zero a zero of order p-1, and divide by (p-1)!. All but one boundary contribution disappear modulo p; the remaining one is nonzero. Thus an elementary congruence replaces positivity when the signs of several coefficients are uncontrolled.

The exponential multiplication law is essential:

    exp(alpha)exp(beta)=exp(alpha+beta).

It converts products of hypothetical algebraic-conjugate relations into another finite sum of the same kind. This makes polynomial symmetry compatible with the analytic family being approximated.

Zeta values do not have this multiplication law. Instead, for integers a,b>1, absolute convergence and splitting a double sum into n>m, n<m, and n=m give

    zeta(a)zeta(b)=zeta(a,b)+zeta(b,a)+zeta(a+b),
    zeta(a,b)=sum_(n>m>=1) n^(-a)m^(-b).

This is a rigorous and unexpected link to multiple zeta values: multiplication introduces a second summation depth. For example,

    zeta(2)zeta(3)=zeta(2,3)+zeta(3,2)+zeta(5).

The identity is not an irrationality proof for zeta(5), because it introduces other constants whose required rational relations are not supplied. It explains why blindly copying the finite exponential symmetrization does not automatically isolate an odd zeta value. This is a limitation of that proposed transfer, not a theorem that no useful transfer can exist.

A concrete research direction is to search for finite families of iterated-integral or multiple-zeta relations whose unwanted coefficients cancel while preserving the denominator/decay inequality. Every cancellation needs an exact identity. A relation that merely replaces zeta(5) by new multiple zeta values has not crossed the arithmetic gate.
