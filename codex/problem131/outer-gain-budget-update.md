# Incorporating the checked outer-prime kernel into the ζ(5) normalizer

Status: exact consequence of the independently reviewed local improvement in
`outer-polynomial-correction-kernel.md` and
`outer-correction-adversarial-review.md`, together with the existing complete
assembly audit. This strengthens an alternative derivation retaining Fauzan's
construction. It is not a first-proof or historical-novelty claim.

## Define the changed normalizer explicitly

Let K=40n, N=3n, h=37n, and let tilde-m_(K,M) be the plain-rank normalizer
in `fauzan-arithmetic-simplification.md`. Retain every old local bound except
in the outer-prime band where the new lemma applies. There set

    d_p=K+3N-3p+1>0,
    gamma_new=gamma_plain+d_p.

The checked lemma actually proves the improvement relative to the source's
stronger min(rank,zero-count) estimate, so this weaker plain-rank-plus-d_p
bound is valid as well. For sufficiently large K, every prime in this band
satisfies the source's outer hypotheses. Define

    m_new = tilde-m_(K,M) / product_{p in the eligible band} p^(d_p).

A normalizing multiplier may be a positive rational number. What matters is
that every coefficient of m_new F_K(X) has nonnegative valuation at every
prime; the local lemma proves exactly that. Thus

    Q_new(X)=m_new F_K(X) belongs to Z[X], deg Q_new<=37n.

This is a genuinely changed normalizer. The old one is not silently assigned
a better bound.

## Exact arithmetic and analytic margins

The fixed-M prime averaging already used in the assembly audit gives

    lim (1/K^2) sum_band d_p log p
       = integral_(1/3)^(49/120) (49/40-3y) dy
       =27/3200.

The finite +1 in d_p contributes only o(K^2). At M=10000 the previous
arithmetic upper bound was 2203747/1600000. The new upper bound is therefore

    2203747/1600000 -27/3200
       =2190247/1600000 =1.368904375.

The independently certified analytic upper bound is -2771/2000. Their sum is

    -26553/1600000.

Since K=40n, the final positive integer-polynomial value obeys

    limsup_(n->infinity) n^(-2) log Q_new(zeta(5)) <= -26.553.

In particular, for all sufficiently large n,

    0<Q_new(zeta(5))<exp(-26n^2).

The degree bound is linear in n. If zeta(5)=a/b were rational, then
b^(37n) Q_new(a/b) would be a positive integer tending to zero, a contradiction.
Every other premise is the same one checked in the complete assembly audit:
positive moments, determinant comparison, small-prime control, coefficientwise
integrality, and the energy certificate. This note does not replace those proofs.

A coarser analytic target -137/100 would already suffice with this normalizer:
the remaining exact margin is 1753/1600000, giving eventual exp(-n^2) after
lower-order terms are absorbed. We have not claimed a new shorter certificate
for that weaker target; the existing rational certificate proves the stronger
one and remains available.

## Role in simplification

The outer-kernel proof adds a finite-field jet argument. Its stronger margin
may permit other simplifications, but a larger numerical margin alone does
not establish that the full exposition is shorter. The previous plain-rank
proof and its eventual exp(-13n^2) bound remain valid without this extra lemma.
The new lemma is an optional, independently checked arithmetic strengthening.
