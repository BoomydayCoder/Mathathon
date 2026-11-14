# Exact p-unit classification for prime-power arity

Status: independently derived complete proof, computationally checked; historical novelty unresolved and no claim of novelty. This is a concrete specialization of the earlier forced-carry framework. It generalizes the familiar odd-Catalan index pattern to every full p^e-ary tree count.

**Theorem.** Let p be prime, e>=1, q=p^e, and

T_q(n)=binom(qn,n)/((q-1)n+1),  n>=0.

Then

v_p(T_q(n)) = s_p((q-1)n)/(p-1) - v_p((q-1)n+1),

where s_p is the base-p digit sum. In particular,

p does not divide T_q(n) iff n=(q^h-1)/(q-1) for an integer h>=0.

Thus the p-unit indices are exactly 0 and the positive base-q repunits 1,11_q,111_q,... . There are exactly H+1 such indices among 0<=n<q^H for H>=1 (h=0,...,H). This is a complete classification, not merely a sufficient family.

**Proof.** The elementary factorial valuation formula is v_p(N!)=(N-s_p(N))/(p-1), proved by counting multiples of p^j. Since multiplication by q=p^e only shifts base-p digits, s_p(qn)=s_p(n). Therefore

v_p binom(qn,n)
= [s_p(n)+s_p((q-1)n)-s_p(qn)]/(p-1)
= s_p((q-1)n)/(p-1).

Division by (q-1)n+1 gives the displayed valuation.

Put a=v_p((q-1)n+1) and write

(q-1)n = Q p^a + (p^a-1),  Q>=0.

This equality includes a=0, when the last term is zero. The trailing a digits are p-1, so

s_p((q-1)n)=s_p(Q)+a(p-1),
v_p(T_q(n))=s_p(Q)/(p-1).

Thus valuation zero is equivalent to Q=0, or (q-1)n+1=p^a. In that case p^e-1 divides p^a-1. Writing a=eh+r with 0<=r<e gives p^a-1 congruent p^r-1 modulo p^e-1. For 0<r<e the latter lies strictly between 0 and p^e-1, a contradiction; hence r=0. Conversely a=eh makes the quotient n=(q^h-1)/(q-1) integral and gives valuation zero. The exceptional modulus p^e-1=1 occurs only at p=2,e=1, where e|a is automatic and the same conclusion holds. QED.

The expression s_p(Q)/(p-1) is an integer: (q-1)n is divisible by p-1, and p^a congruent 1 mod (p-1), so Q is divisible by p-1; digit sums preserve residue modulo p-1. This also checks the positivity/integrality of the tail formula directly.

Examples:

* Binary trees, q=2, p=2: odd Catalan numbers exactly at n=0,1,3,7,15,... .
* Ternary trees, q=3, p=3: 3-units exactly at n=0,1,4,13,40,... .
* Quaternary trees, q=4, p=2: odd counts exactly at n=0,1,5,21,85,... .
* Nine-ary trees, q=9, p=3: 3-units exactly at n=0,1,10,91,... .

The tree interpretation is the full q-ary tree count proved in valuation.md. The present valuation proof remains arithmetic; it does not yet construct an orbit decomposition of trees proving the classification combinatorially. One sharpened project goal is now clear: find a tree operation pairing/grouping every q-ary tree size except the base-q repunit sizes into p-orbits, and account for the exceptional sizes modulo p. The exact classification supplies a target and a check against proposed actions, not the action itself.

Connection to Apéry numbers: both have sparse p-unit index sets described by digits. For T_{p^e}, the allowed indices are base-p^e repunits, imposed by exact removal of all trailing carries. For Apéry A_n, allowed indices have every base-p digit in a finite allowed set, imposed by factorization of the carry-free kernel. See ../problem131/apery-digit-arithmetic.md. These are different digit languages arising from the same factorial/binomial arithmetic.
