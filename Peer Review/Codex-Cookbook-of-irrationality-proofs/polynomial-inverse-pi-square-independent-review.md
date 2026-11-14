# Independent audit: finite polynomial inverses and the direct π² proof

Status: **accepted** for the claims reviewed below. This review covers
[polynomial-inverses-and-boundary-pairings.md](polynomial-inverses-and-boundary-pairings.md)
and Section 2 of [pi-powers.md](pi-powers.md). It does not repeat the separate review of the higher-power transcendence argument. No historical novelty is asserted.

## Finite inverse and degree-raising lemma

Differentiation is nilpotent on each bounded-degree polynomial space, so both displayed inverse sums terminate exactly. The signs in (1+D)^(-1) and (a+bD²)^(-1) are correct. For the latter, a!=0 is essential and is stated.

The degree-raising argument is valid for a **linear** operator B, as intended in the note. Every nonzero polynomial of degree below n has its degree raised by one, so B is injective there and its image contains no nonzero constant. Dimension then gives the claimed direct sum and the unique remainder A. The annihilating functional E is uniquely determined by its value on the constant complement and zero on the image.

With the explicitly supplied k satisfying Bk=1, the function g=r+Ak satisfies Bg=P_n. Orthogonality removes r from <P_n,g>. Coercivity makes <Bg,g> strictly positive because g is nonzero. This both proves A!=0 and determines the sign of <P_n,k>. Cauchy–Schwarz gives the exact stated upper bound; symmetry of B is not additionally required for this argument beyond whatever is used to establish its actual quadratic form.

The Apéry specialization uses n>=1, so its coercivity constant n(n+1) is positive. The coefficient functional E_n annihilates Bx^j because
(j+1)^2 C_(n,j+1)=[n(n+1)-j(j+1)]C_(nj).
Thus its normalization is exactly the same sum of squares as the kernel's arithmetic coefficient. The note correctly keeps kernel existence, boundary regularity, arithmetic integrality, and growth as separate obligations.

For maximum literal clarity, the general lemma should say “linear operator B,” and its Apéry specialization should state “n>=1.” These are the intended hypotheses already used by the formulas, not additional mathematical gaps.

## Fixed-interval π² argument

Assume pi²=a/b with positive integers a,b. Expanding
f_n=x^n(1-x)^n/n! shows that every nonzero endpoint derivative is an integer:
its order is n+k and its value at zero is
(-1)^k binom(n,k)(n+k)!/n!. The reflection x->1-x gives the same property at one. Therefore every term in the displayed G_n has integral endpoint values.

The finite sum cancels exactly to
bG_n''+aG_n=a^(n+1)f_n.
Dividing by b and using pi²=a/b gives
G_n''+pi²G_n=pi²a^n f_n, with the correct exponent of a.

The derivative of G_n' sin(pi x)-pi G_n cos(pi x) is
(G_n''+pi²G_n)sin(pi x). Its upper-minus-lower boundary value is
pi[G_n(1)+G_n(0)], because cos(pi)=-1 and cos(0)=1. Consequently
G_n(0)+G_n(1)=pi a^n integral f_n sin(pi x), with positive sign and no missing factor of pi or b.

The right side is strictly positive on every n>=1 and at most
pi(a/4)^n/n!, which tends to zero with a fixed under the rationality assumption. The left side is an integer. This is the required contradiction, proving irrationality of pi² directly. The method does not require Euler's zeta(2) evaluation, Parseval, or an already known irrationality theorem for pi.

The associated exact checker verifies the finite inverse and integer endpoint polynomials for free parameters a,b at n=1,...,12. The all-index proof is the finite cancellation and endpoint coefficient argument above.
