# Independent audit: the Apéry digit polynomial is a square or discriminant times a square

Status: complete proof from the already established Apéry recurrence and the independently proved Lucas/reflection identities in apery-digit-arithmetic.md. The characteristic-zero differential identity is rederived below; no reliance on modular-form pattern fitting. Historical novelty unresolved. The Legendre-symbol rule mentioned at the end is only an experimental proposal here.

## Theorem

Let p>=5 be prime, A_n=sum_k binom(n,k)^2 binom(n+k,k)^2,

P_p(t)=sum_{r=0}^{p-1} A_r t^r in F_p[t],
Delta(t)=1-34t+t^2.

Then exactly one of the following holds:

P_p=Q_p^2, or P_p=Delta Q_p^2,

with Q_p in F_p[t] and Q_p(0)=1. All roots of Q_p are simple, and gcd(Q_p,Delta)=1. The degrees of Q_p are respectively (p-1)/2 and (p-3)/2. Thus each ordinary root of P_p has multiplicity exactly two; its only possible simple roots are the two roots of Delta, which either both occur or both do not occur.

## 1. Derive the exact first integral before reducing modulo p

The Apéry recurrence gives, with theta=t d/dt,

L=theta^3-t(34theta^3+51theta^2+27theta+5)+t^2(theta+1)^3,
L[A]=0 for A(t)=sum A_n t^n in Q[[t]].

The recurrence/binomial-sum identification is proved in the project's telescoping exposition; the present audit starts from that established input. Define

a=t^2 Delta,
b=2t(1-51t+2t^2)=a',
c=t(t-10),
I[y]=a(2yy''-y'^2)+b yy'+c y^2.

Differentiate, using b=a'. The y'^2 terms cancel and yield the exact polynomial differential identity

I'[y]=y[2a y'''+3a' y''+(a''+2c)y'+c'y]=2y L[y]/t.

For verification, L[y]/t has coefficients

t^2 Delta y'''+3t(1-51t+2t^2)y''+(1-112t+7t^2)y'+(t-5)y,

and a''+2c=2-224t+14t^2, c'=2t-10. Thus I'[A]=0 in Q[[t]], and I[A](0)=0 since a,b,c vanish at t=0. Consequently I[A]=0 exactly. All coefficients of A and I are integers, so the identity reduces modulo every prime.

**Important audit point:** proving only I'[A]=0 after reduction would be insufficient, because series in t^p have derivative zero in characteristic p. The characteristic-zero identity must be established first, as above.

## 2. Transfer the identity to the digit polynomial

Lucas factorization gives A(t)=P_p(t)G(t) with G(t)=A(t^p) in F_p[[t]]. Derivatives of G vanish in characteristic p, so A'=P_p'G and A''=P_p''G. Every term in I is homogeneous quadratic in y and its derivatives. Therefore

0=I[A]=G^2 I[P_p].

The series G has constant term 1, so is invertible; I[P_p]=0 as a polynomial.

Also A_{p-1}=A_0=1 modulo p and A_{p-1-r}=A_r modulo p by the reflection theorem proved independently in apery-digit-arithmetic.md. Hence P_p has degree exactly p-1, constant term 1, and is reciprocal:

t^{p-1}P_p(1/t)=P_p(t).

## 3. Analyze multiplicities over an algebraic closure

Let rho be a root of P_p of multiplicity m. Then rho is nonzero, since P_p(0)=1, and 1<=m<=p-1. Write x=t-rho and P_p=h x^m+O(x^{m+1}), h!=0.

If Delta(rho)!=0, the coefficient of x^{2m-2} in I[P_p] is

a(rho) h^2 [2m(m-1)-m^2]
= a(rho) h^2 m(m-2).

The terms b PP' and c P^2 have strictly larger order. Since a(rho)!=0 and 1<=m<p, the vanishing forces m=2 in F_p, hence m=2 as an integer in this range.

If Delta(rho)=0, its root is simple: the discriminant of Delta is 1152=2^7*3^2, nonzero modulo p>=5. Now a has a simple zero at rho and b(rho)=a'(rho)=rho^2 Delta'(rho)!=0. The coefficient of x^{2m-1} is

h^2 a'(rho)[m(m-2)+m]
= h^2 a'(rho) m(m-1).

The c P^2 term has larger order. Thus m=1, again using 1<=m<p.

This treatment also covers m=1 correctly: the derivatives' leading expressions may have zero coefficients, but the displayed order coefficient remains valid and is decisive.

## 4. Pair the exceptional roots and descend the square

The roots alpha,beta of Delta are distinct, nonzero, and alpha beta=1. Reciprocity of P_p makes their multiplicities equal, so either neither is a root or both are simple roots. Let epsilon be 0 or 1 accordingly. The polynomial R=P_p/Delta^epsilon belongs to F_p[t], because divisibility by the monic polynomial Delta over the algebraic closure is equivalent to zero remainder over F_p. Every root of R has even multiplicity exactly two.

Factor R over F_p into monic irreducibles. Finite fields are perfect, so irreducibles are separable; their exponents are exactly two. Hence R=s H^2 for H in F_p[t] and s in F_p^*. Since R(0)=1, H(0)!=0 and s H(0)^2=1. Defining Q=H/H(0) gives R=Q^2 and Q(0)=1. The root multiplicities show Q is squarefree and coprime to Delta. Degree counting gives the asserted degrees. QED.

## 5. Consequences and limits

Combining P_p=Delta^epsilon Q_p^2 with the digit-series equation from apery-digit-arithmetic.md gives

A(t)^{p-1}=Delta(t)^{-epsilon} Q_p(t)^{-2} in F_p[[t]].

Thus the same Delta that controls the characteristic-zero symmetric-square equation also controls the parity of multiplicities of the finite digit polynomial. This is a precise arithmetic/differential bridge. One may motivate it by the elliptic branch points in the broader project, but no elliptic parametrization is required for this proof.

Parent-agent computations suggested epsilon=(1-(-6|p))/2 for primes checked, where (-6|p) is a Legendre symbol. This audit does **not prove** that rule. The proof above establishes only epsilon in {0,1}, the square/Delta-square shape, and exact root multiplicities. A Legendre-symbol classification needs an additional argument.

No historical novelty claim: this theorem may be known in the literature on Apéry congruences, Picard–Fuchs equations, or Hasse invariants. That literature has not yet been systematically checked here.
