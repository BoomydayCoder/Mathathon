# An elementary involution proof of the Apéry character rule

Status: the character rule is proved by a rational differential-equation symmetry, Lucas congruences and elementary constant-term identities. Elliptic curves and CM are an optional interpretation, not dependencies of the main proof. Exact checks pass in `../computations/check_apery_cm_character.py`; independent review of the new involution argument is requested. No historical novelty claim: the result belongs naturally to established work on Apéry congruences.

## 1. The precise theorem

For a prime p≥5, put

Aₙ=Σₖ binom(n,k)²binom(n+k,k)²,

Pₚ(t)=Σₙ₌₀ᵖ⁻¹Aₙtⁿ ∈ Fₚ[t],   Δ(t)=1−34t+t².

The independent elementary proof in [apery-digit-square-review.md](apery-digit-square-review.md) establishes

Pₚ(t)=Δ(t)^ε Qₚ(t)²,   ε∈{0,1},   Qₚ(0)=1,

where Qₚ is squarefree and coprime to Δ. The missing classification is

**ε=(1−(−6/p))/2.**

Thus Pₚ is a square when −6 is a quadratic residue modulo p, and Δ times a square when −6 is a nonresidue. The proof below explains exactly why the integer −6 appears. Equivalently, Pₚ is a square for p≡1,5,7,11 (mod 24), and Δ times a square for p≡13,17,19,23 (mod 24).

## 2. An auxiliary integer sequence and an exact bridge

Let

bₙ=Σₖ binom(n,k)²binom(2k,k),   B(x)=Σₙ≥₀bₙxⁿ.

Its first coefficients are 1,3,15,93,639. The recurrence is

(n+1)²bₙ₊₁−(10n²+10n+3)bₙ+9n²bₙ₋₁=0.           (1)

Here is an explicit finite certificate, so this identification does not rely on matching initial coefficients. With f(n,k)=binom(n,k)²binom(2k,k), set

G(n,k)=f(n+1,k) k³(3k−4n−4)/(n+1)²,

using the usual zero convention for binomial coefficients outside their range. Then

(n+1)²f(n+1,k)−(10n²+10n+3)f(n,k)+9n²f(n−1,k)
 =G(n,k+1)−G(n,k).

The rational identity is checked symbolically in the script; its endpoints also hold directly. Summing k=0,…,n+1 gives (1), since G(n,0)=G(n,n+2)=0. This proof uses n≥1, with b₀=1,b₁=3 checked separately.

Consequently B is the normalized analytic solution of

x(x−1)(9x−1)B''+(27x²−20x+1)B'+(9x−3)B=0.

The independently proved quadratic pullback in [apery-modular-elliptic.md](apery-modular-elliptic.md) therefore gives the all-orders integer-series identity

A(T(x))=(1−x)B(x)²,   T(x)=x(1−9x)/(1−x).       (2)

The sequence bₙ also has the Lucas property. Write n=pN+r,k=pK+s with 0≤r,s<p. Lucas's theorem applies directly to binom(n,k). If 2s≥p then binom(2k,k) vanishes modulo p; otherwise its first-digit factor is binom(2s,s) and its remaining factor is binom(2K,K). Summing gives bₚN₊ᵣ≡b_Nb_r mod p.

Let Bₚ(x)=Σₙ₌₀ᵖ⁻¹bₙxⁿ. Applying Lucas to both sides of (2), including the already proved Lucas property of Aₙ, gives

**(1−x)^{p−1}Pₚ(T(x))=Bₚ(x)² in Fₚ(x).**       (3)

Indeed B(x)/B(xᵖ)=Bₚ(x), A(t)/A(tᵖ)=Pₚ(t), and T(x)ᵖ=T(xᵖ). This is the exact bridge from the Apéry digit polynomial to an elliptic Hasse polynomial.

## 3. Two elementary special values of Bₚ

Define the Laurent polynomial

Λ(X,Y)=(1+X+Y)(1+X⁻¹+Y⁻¹).

The multinomial theorem gives

CT(Λⁿ)=Σ_{a+b+c=n}binom(n;a,b,c)²=bₙ.

Therefore

Bₚ(x)≡CT[(1−xΛ)^{p−1}] mod p,                 (4)

because (−1)ⁿbinom(p−1,n)≡1. Eliminating the variable Y produces the quartic

Φₓ(X)=[X−x(X²+3X+1)]²−4x²X(X+1)².

More precisely, write 1−xΛ=a−bY−c/Y with

a=1−x(3+X+X⁻¹), b=x(1+X⁻¹), c=x(1+X).

For m=(p−1)/2, the identity

binom(m,j)(−4)ʲ≡binom(2j,j) mod p

shows CT_Y(a−bY−c/Y)^{p−1}≡(a²−4bc)^m. Since Φₓ=X²(a²−4bc), taking CT_X in (4) gives

**Bₚ(x)=[X^{p−1}]Φₓ(X)^{(p−1)/2}.**             (5)

At x=1 the quartic degenerates to Φ₁(X)=(X²−1)². Formula (5), valid as a polynomial identity even at this degenerate value, therefore gives

**Bₚ(1)=1.**

Indeed the relevant coefficient is (−1)^m binom(p−1,m)≡1. At the other end of the polynomial, Lucas's theorem and the same central-binomial identity give

bₚ₋₁≡Σₖ₌₀ᵐbinom(2k,k)
     ≡Σₖ₌₀ᵐbinom(m,k)(−4)ᵏ=(−3)^m=(−3/p).

Thus Bₚ has degree exactly p−1, with leading coefficient (−3/p), and has value 1 at x=1. Neither calculation uses an elliptic-curve theorem.

## 4. A rational involution forces a finite-field functional equation

Define

M(x)=(x−1/9)/(x−1).

It satisfies M(M(x))=x and T(M(x))=T(x): it exchanges the two sheets of the quadratic map T. The second-order equation for B has the exact symmetry

**if y solves (1)'s differential equation, then y(M(x))/(1−x) does too.**

This can be verified by substituting into the displayed differential equation. For reproducibility, write its normalized coefficients as

P(x)=(27x²−20x+1)/[x(x−1)(9x−1)],

Q(x)=(9x−3)/[x(x−1)(9x−1)].

For h=1/(1−x), put P*=P(M)M'−M''/M' and Q*=Q(M)(M')². Direct rational identities are

P*−2h'/h=P,

Q*−P*h'/h+2(h'/h)²−h''/h=Q.

The exact script checks these identities. Their coefficients only require denominators invertible for p≥5, so the symmetry remains valid in characteristic p.

The polynomial Bₚ solves this differential equation over Fₚ. To see this without truncation-boundary assumptions, use B(x)=Bₚ(x)B(xᵖ), whose second factor has zero derivative modulo p; applying the differential operator to this identity yields B(xᵖ) times the operator applied to Bₚ. The factor B(xᵖ) is invertible.

Now set

R(x)=(1−x)^{p−1}Bₚ(M(x)).

It is a polynomial of degree at most p−1. It is also a solution of the differential equation: it equals (1−x)ᵖ times the transformed solution Bₚ(M)/(1−x), and the extra factor has zero derivative.

Any polynomial solution of degree at most p−1 is determined by its constant coefficient, because the coefficient recurrence has leading factors (n+1)²≠0 for n=0,…,p−2. Consequently R=cBₚ for some scalar c. Evaluate at x=1. Only the top-degree term of Bₚ survives on the left, giving

R(1)=bₚ₋₁(8/9)^{p−1}=bₚ₋₁.

Since Bₚ(1)=1, we obtain the exact finite-field identity

**(1−x)^{p−1}Bₚ(M(x))=(−3/p)Bₚ(x).**            (6)

This is the decisive extra symmetry. It converts the empirical sign into a quadratic character by an all-orders argument.

## 5. Read the character at the two fixed points

The fixed points x₀ of M are the two roots of Dₓ=9x²−18x+1. They are distinct for p≥5 and are ordinary points of the differential equation: none is 0,1 or1/9. At either point,

(1−x₀)²=8/9,   M'(x₀)=−1,

so (1−x₀)^{p−1}=(8/9)^{(p−1)/2}=(2/p).

There are exactly two possibilities.

* If Bₚ(x₀)≠0, evaluate (6) at x₀. It gives (2/p)=(−3/p), hence (−6/p)=1.
* If Bₚ(x₀)=0, this zero is simple. Indeed a root of multiplicity m<p at an ordinary point of a second-order equation has lowest coefficient proportional to m(m−1), forcing m=1. Differentiate (6) and evaluate at x₀. The term differentiating the prefactor vanishes because Bₚ(x₀)=0, so −(2/p)Bₚ'(x₀)=(−3/p)Bₚ'(x₀). Therefore (−6/p)=−1.

The alternatives are exhaustive, proving

**Bₚ(x₀)=0 if and only if (−6/p)=−1.**

The bridge (3) gives Pₚ(T(x₀))=0 if and only if Bₚ(x₀)=0, since x₀≠1. Also

Δ(T(x))=Dₓ(x)²/(1−x)²,

so the critical values T(x₀) are exactly the roots of Δ. The established square-factor theorem now yields

**Pₚ=Δ^{(1−(−6/p))/2}Qₚ².**

This completes the main proof. It uses finite sums, Lucas's theorem, rational differentiation, uniqueness of polynomial solutions below degree p, and quadratic characters. It does not use complex multiplication, Deuring's theorem, modular-form reduction or numerical prime-pattern extrapolation.

## 6. Optional geometric explanation: Hasse invariants and CM −24

The proof above is shorter than the geometric route that led to it. The geometry still explains why the same involution and discriminant arise, so we preserve the independently derived interpretation here. None of the following elliptic facts is required for Sections 1–5.

### The Hasse interpretation

Consider the curve Cₓ: V²=Φₓ(X). Its quartic discriminant is

256x⁶(x−1)³(9x−1).

When x is distinct from 0,1,1/9 and the characteristic is not 2 or 3, its smooth projective completion has genus one and has the rational point (X,V)=(0,x); hence it is an elliptic curve.

For completeness, (5) has the usual Hasse interpretation without invoking a mysterious coefficient rule. The holomorphic differential is ω=dX/V. Since V^{p−1}=Φₓ^m, write ω=Φₓ^m dX/Vᵖ. The Cartier operator keeps exactly the powers X^{jp−1}. The degree of Φₓ^m is 2p−2, so only X^{p−1} contributes. Thus

Cartier(ω)=Bₚ(x)^{1/p}ω

over an algebraic closure. An elliptic curve is supersingular exactly when Cartier vanishes on its one-dimensional space of holomorphic differentials. Equivalently, it is supersingular exactly when its geometric p-torsion has only the identity point. These standard equivalent criteria are the only Hasse-theoretic input used here. We have therefore proved

**Bₚ(x)=0 ⇔ Cₓ is supersingular.**               (H)

### Identify the exceptional curves with CM discriminant −24

The critical points of T are the roots of

Dₓ=9x²−18x+1.

They map to the two roots of Δ. The modular coordinate already proved in the companion note is

x(τ)=η(τ)⁴η(6τ)⁸/[η(2τ)⁸η(3τ)⁴]
     =q−4q²+10q³+⋯.

At τ₀=i/√6, its value is the smaller root x₋=1−2√2/3. This follows from Fricke invariance, the absence of interior ramification in the x-coordinate, and the exact positive product bound in the companion note.

We must identify the elliptic curve, not merely guess that the critical point is a CM point. Its explicit j-invariant is

j(Cₓ)=(9x⁴−12x³+30x²−12x+1)³/[x⁶(x−1)³(9x−1)].  (J)

The binary-quartic invariant formula proves (J); the script checks its numerator, denominator and discriminant algebraically. Moreover,

**j(C_{x(τ)})=j(6τ).**                          (JM)

Here is a finite proof of this modular-function identity. Both sides are meromorphic functions on X₀(6). In the cusp order c=1,2,3,6 their pole orders are at most 1,2,3,6. For the left side this follows from the coordinate values x=1/9,∞,1,0 and (J). For j(6τ) it follows from the standard modular transformation and cusp widths, or from its expression E₄(6τ)³/η(6τ)²⁴. Neither side has interior poles.

At infinity the exact expansion of the left side, computed from eta products, is

q⁻⁶+744+196884q⁶+O(q⁷),

which equals j(6τ) through q⁶. Their difference thus has a zero of order at least 7 at infinity, and possible poles of total order at most 1+2+3=6 elsewhere. A nonzero meromorphic function on a compact curve has equal total numbers of zeros and poles. The difference must vanish. This justifies the finite certificate; it is not extrapolation from a long series match.

At τ₀, (JM) identifies C_{x₋} with a complex elliptic curve whose period ratio is 6τ₀=i√6. Its lattice Z+Z·i√6 is preserved by multiplication by √−6. The maximal order in Q(√−6) is exactly Z[√−6], of discriminant −24, so this elliptic curve has CM by that maximal order. The other exceptional curve C_{x₊} is its algebraic conjugate over Q, and therefore has the same CM endomorphism field. Equivalently, reducing (J) modulo Dₓ gives

j(Cₓ)=4978368−2560896x,

so the two j-invariants satisfy

J²−4834944J+14670139392=0.

This last polynomial is a consistency check, not a replacement for the explicit modular identification (JM).

Both exceptional curves have good reduction at every prime above p≥5. The critical polynomial has distinct roots, and neither can reduce to 0,1 or1/9: substituting these values into Dₓ gives 1,−8 and−8/9. Their quartic discriminants therefore remain nonzero.

### A direct proof of the corresponding CM reduction criterion

The usual general name is Deuring's CM reduction criterion. For our situation we can prove the needed implication using elementary properties of elliptic endomorphisms rather than importing its full quaternion-algebra formulation.

The standard specialization fact required is this: after enlarging a number field so that the endomorphisms are defined, an elliptic curve with good reduction specializes its characteristic-zero endomorphisms to endomorphisms of the reduction, preserving composition, degree and the identity morphism. For α in a CM order, its degree is its field norm and α+ᾱ=[Tr α]. This is the familiar good-reduction endomorphism lemma; it is stated in the readable reference below, Section 21.5. We also use finiteness of the ideal class group of an imaginary quadratic field.

Let K=Q(√−6), O_K=Z[√−6], and let p≥5.

**If p is inert, the reduction cannot be ordinary.** On an ordinary elliptic curve over an algebraic closure of characteristic p, the geometric p-torsion is the cyclic group Z/pZ. The reduced endomorphism √−6 would act on this group by some scalar a∈Fₚ. Its equation α²=[−6] would imply a²=−6 in Fₚ, impossible when (−6/p)=−1. Thus the reduction is supersingular.

**If p splits, the reduction cannot be supersingular.** Write (p)=𝔭𝔭̄. By finiteness of the ideal class group there is an h≥1 with 𝔭ʰ principal, say (β)=𝔭ʰ. Then Norm β=pʰ, and p does not divide Tr β: β vanishes modulo 𝔭 but β̄ does not, so their sum is nonzero modulo 𝔭. In a supersingular reduction there is no nontrivial geometric p-torsion. Every endomorphism of degree a power of p is therefore purely inseparable: its separable kernel would otherwise contain a nontrivial point of p-power order. Both reduced β and β̄ have degree pʰ, so both have zero tangent map. Their sum [Tr β] would have zero tangent map as well, implying p divides Tr β. This is a contradiction.

Hence the good reduction of either exceptional CM curve is supersingular exactly when (−6/p)=−1. This direct proof makes explicit which part of the classical CM criterion is needed; it does not assume the digit-polynomial character rule.

## 7. Source and verification boundaries

Primary modular source: F. Beukers, *Irrationality proofs using modular forms*, Astérisque 147–148 (1987), 271–283; downloaded and cited in [apery-modular-elliptic.md](apery-modular-elliptic.md). The actual period, quartic, Hasse coefficient identity and finite modular-function certificate above are derived explicitly here.

Historical source for CM reduction: M. Deuring, *Die Typen der Multiplikatorenringe elliptischer Funktionenkörper*, Abhandlungen aus dem Mathematischen Seminar der Universität Hamburg 14 (1941), 197–272, DOI [10.1007/BF02940746](https://doi.org/10.1007/BF02940746). The publisher's nominal PDF endpoint returned an HTML landing page, saved as `../sources/deuring-1941.html`; the full original paper was **not** read. Its bibliographic metadata was verified, but no theorem-number quotation from that unread paper is asserted.

Readable standard-theorem reference actually retrieved: A. V. Sutherland, MIT 18.783, Fall 2023, [Lecture 21, Section 21.5](https://math.mit.edu/classes/18.783/2023/LectureNotes21.pdf), for specialization of endomorphisms, and [Lecture 13](https://math.mit.edu/classes/18.783/2023/LectureNotes13.pdf), for ordinary/supersingular elliptic curves. Local files are `../sources/sutherland-cm-lecture21.pdf` and `../sources/sutherland-ordinary-supersingular-lecture13.pdf`, with extracted text. These are expository references, not misidentified as the primary Deuring paper. The optional geometric section supplies the particular reduction argument instead of relying on a source statement not retrieved.

The exact script proves its polynomial/rational identities and checks the Hasse identity at small primes and the character pattern for all primes 5≤p<150. The main theorem follows from the elementary involution proof, independently of the standard elliptic facts used only in the optional interpretation. The finite prime tests do not replace that proof. Its final independent review remains a separate checkpoint; historical novelty is unresolved.
