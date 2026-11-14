# Powers of pi: elementary integral proofs versus explicit transcendence machinery

Status: independently derived arguments, with all logical steps displayed; sections 3 and 5 independently audited with no gap found (see [pi-powers-review.md](pi-powers-review.md)). This is mathematical review, not formal verification. The proof for arbitrary positive integral powers is a classical Hermite–Lindemann-type construction, not a historical novelty claim. It is more substantial algebraically than the short Niven integral proof.

## 1. Audit of the comparison proofs

The proof for exp(a/b), a,b positive integers, in irrationality-framework.md is correct with

    f_n(x)=x^n(a-bx)^n/n!.

No additional power of b is needed: at x=a/b-y the numerator becomes (a-by)^n y^n. Thus derivatives at both endpoints are integers. With F_n=sum(-1)^j f_n^(j), the derivative of exp(x)F_n is exp(x)f_n. Both the sign and factorial smallness are correct. Negative nonzero rational exponents follow by taking reciprocals.

The exact Padé normalization is

    Q_n(z)exp(z)-P_n(z)
       =(-1)^n/n! integral_0^z exp(t)t^n(z-t)^n dt,
    P_n(z)=sum_{j=0}^n (2n-j)!/[j!(n-j)!] z^j,
    Q_n(z)=P_n(-z).

In particular the leading remainder coefficient is (-1)^n n!/(2n+1)!. For n=1, (2-z)exp(z)-(2+z) begins with -z^3/6, which confirms the easily missed minus sign. Exact symbolic tests for n=0,...,8 agree with the general endpoint proof.

The zeta(2) monomial integrals are

    J_kk=zeta(2)-H_k^(2),
    J_kl=(H_k-H_l)/(k-l) for k!=l.

The numerator and denominator in the off-diagonal formula have the same sign. The n-fold integration-by-parts sign is (-1)^n. The maximum of x(1-x)y(1-y)/(1-xy) is g^5, g=(sqrt(5)-1)/2, and the shared proven lcm rate c=(14log2+9log3+5log5)/25 is sufficient because 2c+5log g<0. The framework now cites the self-contained lcm appendix instead of relying on the unproved sharper bound D_n<3^n.

## 2. A short direct proof for pi squared

The Niven method also proves pi squared irrational by ordinary calculus alone.
Here is the fixed-interval version, which avoids both Parseval and zeta(2).
The original1947 one-page paper proves pi irrational using the same finite
inverse of D²+1; its scan is now archived at [niven-pi-1947.pdf](../sources/niven-pi-1947.pdf),
with URL/hash provenance. This adaptation and method are classical, not a
historical novelty claim.

Suppose pi²=a/b with positive integers a,b. For n>=1 define

    f_n(x)=x^n(1-x)^n/n!,
    G_n(x)=sum_(j=0)^n (-1)^j a^(n-j)b^j f_n^(2j)(x).

Every derivative of f_n at0 is an integer: its powers start at degree n,
and each nonzero endpoint derivative has a factorial j!/n!. Symmetry
f_n(1-x)=f_n(x) gives the same at1. Thus G_n(0),G_n(1) are integers.
The finite derivative sum telescopes:

    b G_n''+a G_n=a^(n+1) f_n,
    G_n''+pi² G_n=pi² a^n f_n.

Multiply by sin(pi x) and integrate the elementary Green identity:

    derivative[G_n' sin(pi x)-pi G_n cos(pi x)]
       =(G_n''+pi² G_n)sin(pi x).

At0 and1, sine vanishes and cosine is respectively1 and-1. Consequently

    G_n(0)+G_n(1)=pi a^n integral_0^1 f_n(x)sin(pi x)dx.

The right side is strictly positive and at most

    pi (a/4)^n/n! ->0,

because0<x(1-x)<=1/4. This is a positive integer tending to zero,
a contradiction. Therefore pi² is irrational.

**Why the construction works.** On polynomials D² is nilpotent. The finite
sum defining G_n is the exact inverse of bD²+a, after multiplication by the
integer a^(n+1). Its boundary pairing with the solution sin(pi x) converts
that algebraic identity into a small positive integer. The auxiliary
polynomial's endpoint zeros and factorial scaling are chosen to make every
boundary derivative integral. This is the same design issue seen in the
exp(r) proof, expressed through a second-order differential equation.

An alternative connection remains useful: the completed Beukers-type proof
for zeta(2), combined with Euler's pi²=6zeta(2), gives another proof. Parseval
for x on[-pi,pi] supplies Euler's identity from sine coefficients
2(-1)^(n+1)/n. That alternative uses a standard analysis theorem which the
direct argument above does not need.

Neither argument by itself proves pi^4 irrational: irrational numbers may
have rational squares. The uniform nonzero-rational-power statement below
uses the stronger classical construction and the stated reduction.

## 3. Every positive integral power: an explicit proof

**Theorem.** For every positive integer r, pi^r is irrational.

The following proof avoids citing the transcendence of pi as a black box. It uses polynomial symmetry, elementary algebraic-integer facts, and integration of entire functions along finitely many line segments. These are undergraduate algebra/analysis tools, but the construction is the core of a classical transcendence proof, not a lighter Niven-style shortcut.

**Scope for rational powers.** The integer-power theorem immediately implies
that pi^r is irrational for every nonzero rational r. Write r=p/q with
q>0 and p!=0. A rational pi^(p/q) would make pi^p rational after taking
its qth power; for p<0, take its reciprocal to contradict the positive
integer-power theorem. This does not assert the result for every real r:
for example r=log(2)/log(pi) gives pi^r=2. Thus the intended uniform
power statement has an explicit nonzero-rational range.

### Step A: manufacture an exponential relation stable under algebraic symmetry

Suppose pi^r=a/b with positive integers a,b. Let alpha_1,...,alpha_d be all roots, with d=2r, of

    h(X)=X^(2r)-(-1)^r(a/b)^2.

One root is i pi. Therefore

    product_{j=1}^d (1+exp(alpha_j))=0,

since one factor is 1+exp(i pi)=0. Expand the product over all subsets S of {1,...,d}. Write beta_S=sum_{j in S} alpha_j. Separate the zero subset sums, whose number is a positive integer k (the empty subset already contributes one). Enumerate the nonzero subset sums with multiplicity as beta_1,...,beta_m. Then

    k+sum_{j=1}^m exp(beta_j)=0.                         (A)

The polynomial product_S (X-beta_S) has rational coefficients. Indeed each coefficient is a symmetric polynomial with integer coefficients in alpha_1,...,alpha_d; the elementary symmetric polynomials in the alpha_j are rational, being coefficients of h. Removing the factor X^k leaves a monic rational polynomial with roots beta_1,...,beta_m and nonzero constant coefficient.

Choose a positive integer c clearing all its coefficients, and put

    g(X)=c product_{j=1}^m (X-beta_j) in Z[X].

Its leading coefficient is c and g(0)!=0. Repeated beta_j cause no problem.

**Motivation.** Multiplying over conjugate roots creates a relation whose algebraic arguments can be summed symmetrically. This converts algebraic integers arising later into ordinary integers. A single relation 1+exp(i pi)=0 would not by itself supply that arithmetic conversion.

### Step B: use endpoint vanishing to make an integer nonzero modulo a prime

Choose any sufficiently large prime p not dividing k c g(0), and define

    f_p(X)=c^(mp) X^(p-1) g(X)^p/(p-1)!,
    F_p(X)=sum_{ell=0}^{(m+1)p-1} f_p^(ell)(X).

At zero, derivatives below p-1 vanish. The derivative of order p-1 is

    c^(mp)g(0)^p.

Every higher derivative at zero is an integer divisible by p, since ell!/(p-1)! is divisible by p for ell>=p. Consequently

    F_p(0) in Z,
    F_p(0) ≡ c^(mp)g(0)^p mod p.                       (B)

At each beta_j, the polynomial g(X)^p has a zero of order at least p, so derivatives of order less than p vanish. We claim every remaining f_p^(ell)(beta_j) is p times an algebraic integer.

Here is the denominator check. The number c beta_j is an algebraic integer because it satisfies the monic integer polynomial

    c^(m-1)g(Y/c)=Y^m+g_(m-1)Y^(m-1)+c g_(m-2)Y^(m-2)+...+c^(m-1)g_0.

Let H(X)=X^(p-1)g(X)^p, an integer polynomial of degree D=(m+1)p-1. For ell>=p,

    f_p^(ell)(beta_j)/p
      =sum_{v>=ell} H_v binom(v,ell)
          [ell!/(p(p-1)!)] c^(mp) beta_j^(v-ell).

The bracketed factor is an integer. Moreover v-ell<=D-p=mp-1. Thus

    c^(mp) beta_j^(v-ell)
       =c^(mp-v+ell)(c beta_j)^(v-ell)

is an algebraic integer. The claim follows.

The sum sum_j F_p(beta_j)/p is rational: it is a rational-coefficient polynomial evaluated and summed over all roots of the rational polynomial g, hence is rational by Newton's identities (or symmetry). It is also an algebraic integer, being a sum of algebraic integers. A rational algebraic integer is an integer. Therefore

    sum_j F_p(beta_j) belongs to p Z.

Combining this with (B),

    N_p=kF_p(0)+sum_j F_p(beta_j) is an integer,
    N_p ≡ k c^(mp)g(0)^p mod p !=0.                    (C)

This proves nonvanishing before any analytic estimate.

**Motivation.** The zero at 0 has multiplicity p-1, while all other endpoint zeros have multiplicity p. After factorial division, precisely one endpoint derivative survives modulo p. This is the arithmetic analogue of arranging a positive integrand: it prevents the small quantity from being zero.

### Step C: show the same integer tends to zero

Since F_p-F_p'=f_p,

    (exp(-z)F_p(z))'=-exp(-z)f_p(z).

Integrate along the straight segment from 0 to beta_j and multiply by exp(beta_j):

    F_p(beta_j)-exp(beta_j)F_p(0)
      =-exp(beta_j) integral_0^(beta_j) exp(-z)f_p(z) dz.

Sum over j and use (A). This gives

    N_p=-sum_j exp(beta_j) integral_0^(beta_j) exp(-z)f_p(z) dz.  (D)

All segments and all beta_j are fixed as p varies. Choose M>=1 bounding |z| and |beta_j| on their union, and G>=1 bounding |g(z)| there. Then

    |f_p(z)| <= [c^m M G]^p/(p-1)!.

The right side of (D) is bounded in absolute value by

    m M exp(2M) [c^m M G]^p/(p-1)! ->0.

There are arbitrarily large primes avoiding the finitely many divisors of k c g(0), so (C) and (D) give nonzero integers tending to zero. Contradiction. Hence pi^r is irrational. QED.

### What this proof does and does not simplify

No Lindemann theorem was invoked. Nevertheless, the same argument works if i pi is any algebraic number satisfying exp(i pi)=-1, by taking its conjugates instead of the special 2r-root polynomial. It therefore essentially proves transcendence of pi. The conceptual cost has been exposed, not removed.

For a single low power, the zeta(2) integral proof plus Euler's identity is much shorter. For all positive powers at once, either invoke transcendence transparently or teach this classical Hermite–Lindemann construction. It would be misleading to label the latter a new elementary shortcut.

## 4. Reference status

- Ivan Niven, “The Transcendence of π”, American Mathematical Monthly 46 (1939), pp. 469–471. Publisher DOI 10.1080/00029890.1939.11998903; JSTOR DOI 10.2307/2302515. These metadata were retrieved from Crossref in this session. The publisher PDF returned HTTP 403 and the attempted JSTOR PDF returned HTML. Subsequently, Niven's primary text was retrieved and read as Appendix A (PDF pages 19–22) of Jonathan M. Borwein and Scott T. Chapman, *I Prefer Pi: A Brief History and Anthology of Articles in the American Mathematical Monthly*, https://www.carmamaths.org/resources/jon/31415.pdf . The authors state that the appendix reproduces Niven's proof exactly. The local source is codex/sources/pi-monthly-anthology.pdf, with text extraction alongside it. The independently reconstructed proof above matches Niven's subset-sum symmetrization, exceptional p-1 endpoint multiplicity, and prime-divisibility argument. Our prefactor c^(mp) is one extra harmless factor c compared with his c^(mp-1); our algebraic-integer explanation expands his symmetric-polynomial denominator argument.
- Ivan Niven, “A simple proof that pi is irrational”, Bull. Amer. Math. Soc. 53 (1947), 509: see retrieval limitations in irrationality-framework.md.
- The elementary lcm estimate used for zeta(2) is fully proved in zudilin-explained.md section 10; its applicability here was separately checked by the exact inequality 3^18 5^135<2^347.

The primary-source comparison is now complete at the structural level. A further pedagogical simplification removes algebraic integers entirely: if g(X)/c is monic rational, its elementary symmetric root functions have denominators dividing c. Every symmetric integer polynomial of total degree <=s can be expressed, by the fundamental theorem of symmetric polynomials, as an integer polynomial in those elementary symmetric functions, each monomial using at most s factors. Multiplication by c^s therefore clears all denominators. Apply this with s=mp-1 to the sum over beta_j of each derivative polynomial after dividing by p. This is Niven's direct route and uses only polynomial symmetry. The algebraic-integer route above is often conceptually cleaner for readers who already know the definition.

## 5. A further connection: the same endpoint multiplicities prove exp(r) transcendental

The machinery above can also strengthen the elementary exponential argument. This is a classical Hermite-type proof, included to explain the conceptual bridge between irrationality and transcendence rather than to claim a new result.

Fix nonzero rational r and suppose exp(r) is algebraic. It is nonzero, so it satisfies an integer polynomial relation with nonzero constant term:

    a_0+sum_{j=1}^s a_j exp(jr)=0,   a_0!=0.

Terms with zero coefficients may be discarded. Put beta_j=jr for the remaining positive indices, and choose g=c product_j(X-beta_j) in Z[X] with c>0. All beta_j are nonzero and g(0)!=0.

Use the same f_p and F_p from section 3 with m equal to the number of retained beta_j, for primes p not dividing a_0 c g(0). The endpoint arithmetic gives F_p(0)≡c^(mp)g(0)^p mod p, and each F_p(beta_j) is p times an algebraic integer. Since beta_j is rational, F_p(beta_j) is rational too, hence belongs to p Z individually. Therefore

    N_p=a_0F_p(0)+sum_j a_j F_p(beta_j)

is a nonzero integer modulo p. Summing the same integral identity with weights a_j and using the exponential relation shows |N_p|<=C_1 C_2^p/(p-1)! for fixed constants. This tends to zero and contradicts integrality. Hence exp(r) is transcendental.

This proof needs no polynomial conjugates or complex segments: the beta_j are rational and integration is along real segments (possibly negative). It reveals precisely what is added to the short irrationality proof: instead of two endpoints with equal zero orders, use finitely many endpoints and one exceptional order p-1, so a modular argument prevents cancellation among several terms. Positivity is replaced by nonvanishing modulo a prime.
