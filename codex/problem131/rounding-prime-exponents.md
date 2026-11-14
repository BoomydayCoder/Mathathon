# Rational superfactorial weights do not require a power of the determinant

Status: independently checked by the simplify agent; see
[the review](prime-exponent-rounding-independent-review.md). This is elementary
normalization bookkeeping, not a historical-novelty claim. It removes an
unnecessary fixed power from the PNT-free construction without weakening its
certified upper bound.

## 1. The shortest explanation: take the primitive part

**Lemma.** Let F be a nonzero polynomial over Q, q a positive integer, and c a
positive rational number. Suppose c F^q belongs to Z[X]. Write F=gP, where g>0
is rational and P is a primitive integer polynomial (the gcd of its coefficients
is one). Then c g^q is a positive integer. Consequently, for every real theta,

    |P(theta)| <= |c F(theta)^q|^(1/q).                 (1)

Moreover, P and F have the same degree and the same sign at every real point.

**Proof.** Gauss's lemma says that P^q is primitive. A rational multiple d of a
primitive integer polynomial has integral coefficients only if d is an integer:
an integer linear combination of its coefficients equals one, so the same
combination after multiplication equals d. Apply this with d=c g^q. It is
positive, hence at least one. Now

    |c F(theta)^q|^(1/q)=(c g^q)^(1/q)|P(theta)|

proves (1). No estimate on primes, factoring, or asymptotics is involved. QED.

This is a useful design principle. If clearing rational weights is easiest
after raising an auxiliary polynomial to a fixed power, do so in the arithmetic
certificate. The actual small integer polynomial may still be its primitive
part, with its original degree. The power need not survive into the proof.

## 2. An explicit primewise version: round DOWN

Let F be a nonzero rational polynomial and let v_p^G(F) be the minimum p-adic
valuation of its coefficients. Suppose finitely supported rational numbers b_p
satisfy

    v_p^G(F)+b_p >= 0                                (2)

for every prime p. Define the positive rational number

    m = product_p p^floor(b_p).

Because v_p^G(F) is an integer, (2) implies

    v_p^G(F) >= ceil(-b_p) = -floor(b_p).

Therefore mF belongs to Z[X]. Also

    log m <= sum_p b_p log p.                        (3)

Negative b_p cause no exception: floor(-1/3)=-1 and the hypothesis then forces
v_p^G(F)>=1. Rounding upward would also preserve integrality, but would add an
unnecessary cost. The integrality of the valuation is exactly what permits
rounding downward.

The primitive part of F is at least as small in absolute value as this mF,
since every integral scalar multiple of the primitive part has a nonzero
integer scalar. Thus the first lemma is the cleanest existence argument;
the floor construction gives an explicit finite normalization when wanted.

## 3. Apply it to the certified superfactorial product

Use the rational coefficients a_c in the checked sparse majorant. Put

    B(N)=product_(j=0)^(N-1) j!,
    b_(K,p)=sum_c a_c v_p(B(cK)),
    m_K=product_p p^floor(b_(K,p)).

All cK are integers because K is a multiple of40. Only primes at most
(max c)K can contribute. The exact inequality (3) gives

    log m_K <= sum_c a_c log B(cK).                  (4)

Let C=148, M=10^7, H_K=floor(K/M), and retain the audited repair

    E_K=D_(H_K)^(3CK) * ((5K)!)^ceil(sqrt K).

The PNT-free normalizer proof establishes, for every prime and sufficiently
large K,

    v_p^G(E_K Delta_K)+b_(K,p) >= 0.

Thus the explicit polynomial

    P_K(X)=E_K m_K Delta_K(X)

has integer coefficients and degree at most h=37K/40. Its value at zeta(5) is
positive, since every scalar is positive and the moment determinant is positive.

Alternatively, choose q=40000 to clear the certificate's rational coefficients
and put T_K=product_c B(cK)^(q a_c). The already audited statement

    E_K^q T_K Delta_K(X)^q belongs to Z[X]

implies exactly what is needed by the primitive-part lemma. One can define
P_K to be the primitive integer part of Delta_K and avoid introducing m_K at
all. These are two legitimate choices of P_K; they need not be equal.

## 4. The same analytic bound, at the original degree

Write F_K=S_K Delta_K as in the full guide. The explicit floor choice satisfies

    log(E_K m_K/S_K)
       <=(A+3C log4/M)K^2+o(K^2).

The checked upper bound is below1.379. There is no additional rounding cost.
For the primitive choice, take the q-th root of the existing bound on
E_K^q T_K Delta_K(zeta(5))^q; the same inequality follows directly from (1).
Combining this with the separately proved analytic exponent below-1.38 gives
positive integer-polynomial values at zeta(5) tending to zero like exp(-cK^2),
with polynomial degree O(K).

If zeta(5)=a/b with positive integer b, then b^h P_K(a/b) is a positive integer.
Its logarithm is at most -cK^2+O(K), so it tends to zero, a contradiction.

This removes only a bookkeeping complication. The local valuation lemmas,
exact profile certificate, small-prime repair and analytic determinant bound
remain substantive ingredients. The PNT-free argument is computer-assisted
and is not claimed to be formally verified by the separate Lean build.
