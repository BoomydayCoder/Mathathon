# Original zeta(3) Schur family: final proof status and provenance

Checkpoint: 2026-11-14, 15:50 UTC. This note concerns the original zeta(3)
family only. “Accepted” below means independently checked within this research
workspace; it does not mean external publication or peer review.

## Proved results and their proof locations

Write C_h=5A_h-A_(h-1), D_h=5B_h-B_(h-1), where A,B are the classical
Apéry solutions with (A_0,A_1)=(1,5), (B_0,B_1)=(0,6).

| Result | Proof and independent review |
|---|---|
| Original rank-one root q_h=D_h/C_h for every h>=1 | [Heun–RII proof](rank-one-original-apery-heun-rii-proof.md); [root review](original-schur-apery-root-review.md), [second review](rank-one-heun-rii-second-independent-review.md) |
| Elementary finite-node boundary defect, sufficient for the complete rank-one proof | [Real boundary proof](rank-one-real-boundary-proof.md); [first review](rank-one-real-boundary-independent-review.md), [second review](rank-one-real-boundary-second-independent-review.md) |
| Positive logarithmic inverse-norm identity | [Remainder isometry](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/rank-one-logarithmic-remainder-isometry.md); [review](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/rank-one-remainder-isometry-independent-review.md) |
| Exact positive error tail and q_h-zeta(3) asymptotic | [Tail and rate](rank-one-original-error-tail-and-rate.md); [review](rank-one-original-error-independent-review.md) |
| Fixed integer rank-q Heun quotient; rank-two Pell subsequence and explicit Apéry formula | [Rank-q reduction](zeta3-rank-q-heun-reduction.md); [review](zeta3-rank-q-independent-review.md) |
| All-rank basis from rank-one polynomials and a fractional combination of the two Apéry endpoints h and h-q+1 | [Endpoint basis](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/zeta3-all-rank-apery-endpoint-basis.md); [review](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/zeta3-all-rank-independent-review.md) |

The exact error conclusions include

\[
 q_h-\frac{B_h}{A_h}=\frac{6}{h^3A_hC_h},\qquad
 q_h-\zeta(3)\sim4\pi^3(3-2\sqrt2)(17+12\sqrt2)^{-2h},
\]

and (q_h-zeta(3))/(zeta(3)-B_h/A_h) tends to 3+2sqrt(2). The rate had
previously been open for this optimized family; it is no longer a gap.

The rank-two discriminant is 8h²-12h+5. It is a rational square precisely
when (4h-3)²-2d²=-1, giving exact rational brackets along the corresponding
Pell indices. The explicit formula transfers classical Apéry lcm-cube
arithmetic on those rational indices. The rank-q endpoint formula keeps its
finite spectral coefficients; it is not an unevaluated norm being called a
closed two-term recurrence.

## Primary sources: what is actually identified

The earlier [source dictionary](rank-one-primary-source-context.md) gives the
inspected theorem and parameter mappings for the following classical inputs:
Gautschi's rational Gaussian quadrature, Krattenthaler's rationally modified
moment determinant, Ferguson's fixed-measure Hankel approximation, and
Prévost's all-half Wilson evaluation of Apéry's A_n and associated numerator.
The original measure mu_3/D_(2h) is distinct from Prévost's moving-pole all-half
measure. Their equal scalar quotient now follows from two separate proofs.

The operator construction has a precise classical framework. In
Grünbaum–Vinet–Zhedanov, [*Algebraic Heun operator and band-time limiting*](https://arxiv.org/abs/1711.07862),
[DOI 10.1007/s00220-018-3190-0](https://doi.org/10.1007/s00220-018-3190-0),
§4, equation (4.13), the algebraic Heun operator is

\[
 \mathcal M=\tau_1LZ+\tau_2ZL+\tau_3L+\tau_4Z+\tau_0I.
\]

Our formal operator is its specialization

\[
 L=\mathscr D_M,\quad Z=x^2,\quad
 (\tau_1,\tau_2,\tau_3,\tau_4,\tau_0)=(1,1,-1,-\kappa,0).
\]

Their §5 Lemma 1, equations (5.25) and (5.28), explicitly connect finite
cutoffs with vanishing links of a tridiagonal operator. These general
construction principles are therefore established prior context. They do
not by themselves supply the symmetry needed here: our formal Wilson
operator has a nonzero skew boundary form for the original positive measure.
The project proves that defect directly and cancels it before using the
finite-degree eigenspace. No positive-parameter Wilson orthogonality theorem
is applied to parameters (1,1,1,-2h).

The author PDF of Ismail–Rains–Stanton,
[*Orthogonality of very well-poised series*](https://www-users.cse.umn.edu/~stant001/PAPERS/10W9FINAL.pdf),
[DOI 10.1016/j.indag.2024.09.007](https://doi.org/10.1016/j.indag.2024.09.007),
has a varying-weight orthogonality theorem in §6, Theorem 6.5, and an RII
recurrence in §9, Proposition 9.1. The text after that proposition explicitly
points to Ismail–Masson (1995) for the RII framework. This establishes close
classical context for varying-measure quadratic-factor recurrences. No exact
parameter identification of that family with the original zeta(3) Schur
polynomials is claimed here. The present recurrence existence and coefficients
are proved directly, including the old-node versus new-denominator conversion.

The institutional abstract for Ismail–Rahman (2017),
[*Diophantine properties of orthogonal polynomials and rational functions*](https://stars.library.ucf.edu/scopus2015/5565/),
[DOI 10.1090/proc/12355](https://doi.org/10.1090/proc/12355), describes
factorization of associated Wilson/Askey–Wilson polynomials with an application
to isochronous particle systems. That abstract does not report a zeta
irrationality result. Its title alone is insufficient to attribute the
present family to it. The AMS PDF request returned HTTP 403; no full-text
comparison with that paper has been completed.

## Remaining questions, separated from closed proof obligations

1. **Historical attribution.** We have not established whether the exact
   original adjacent-Apéry family, its boundary repair, or the rank-two Pell
   brackets already occur in another form in the literature. Independent
   reconstruction and unsuccessful exact-fraction searches establish no novelty.
2. **Higher-rank arithmetic.** Every root lies in the field of its finite
   Heun eigenvalue and has the proved two-endpoint formula, but uniform heights,
   denominator bounds for general algebraic ranks, and irreducibility are not
   proved. Distinct Heun eigenvalues need not have distinct Schur roots by any
   theorem established here; no such distinctness claim is made for q>=3.
3. **A fully real proof for all ranks.** Rank one, and separately rank two,
   have finite-difference symmetry proofs. The accepted general-rank reduction
   uses the independently checked all-degree boundary identity from its contour
   proof. Removing that analytic input for every q is an optional simplification,
   not a missing hypothesis in the accepted theorem.
4. **Minimal exposition.** The [design guide](schur-heun-apery-design-guide.md)
   explains why the boundary repair, two cutoffs, RII recurrence, and Pell
   discriminant arise. A publication-ready account would still consolidate the
   scattered inputs and perform a broader primary-literature comparison.

## Reproducibility and source custody

The new Heun primary PDF, extracted text, and URL/SHA256/byte-count record are
`../sources/grunbaum-vinet-zhedanov-algebraic-heun-171107862.*`.
The Ismail–Rains–Stanton files and checksum record are
`../sources/ismail-rains-stanton-very-well-poised-2024.*`.
Earlier source custody is in `../sources/rank-one-primary-context-provenance.json`.

The checks accompany proofs; finite checks are not used as all-index evidence.
In particular, `check_rank_one_real_boundary.py` passed all 364 permitted degree
pairs through M=12. The original rank-one bridge and independent Gram checkers
passed their stated ranges. `check_rank_two_apery_formula.py` passed both roots
and the width norm from original moments for h=2,...,12.
`check_zeta3_all_rank_endpoint_basis.py` passed 110 basis-vector pairing
identities for h=2,...,9 and 2<=q<=min(h,6), including exact basis rank,
lower orthogonality, extreme affine pairings, and Heun invariance.
