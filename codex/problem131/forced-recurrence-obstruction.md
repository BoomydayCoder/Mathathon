# Why taking more derivatives does not automatically prove zeta(5)

Status: new independent derivation in this session, with an all-n proof below.
Historical novelty not established. This is a precise obstruction to one natural
extension, not an obstruction to other irrationality methods.

## The tempting construction

Keep Zudilin's rational function R_n and its short certificate S_n=s_n R_n,
where s_n(t)=4(2n+1)(-2t²+t+(2n+1)²). For m>=1 define

    L_(n,m)=-1/(2m)! sum_(t>=1) R_n^(2m-1)(t).

Write R_n=sum_k[a_(n,k)/(t+k)²+b_(n,k)/(t+k)]. Its double-pole coefficients
sum to A_n, and sum b_(n,k)=0. Differentiation shows

    L_(n,m)=A_n*zeta(2m+1)-B_(n,m),
    B_(n,m)=sum_k[a_(n,k) H_k^(2m+1)
                      +b_(n,k)/(2m) H_k^(2m)].

Moreover 2m D_n^(2m+1) B_(n,m) is an integer, since a_(n,k) are integers
and D_n b_(n,k) are integers. Thus a single odd zeta value is isolated, with
explicit controlled denominators. That achievement alone does not prove
irrationality: we still require sufficiently small nonzero linear forms.

## Exact inhomogeneous recurrence

Let T be the Apéry recurrence operator

    (T u)_n=(n+1)³u_(n+1)-(2n+1)(17n²+17n+5)u_n+n³u_(n-1).

**Theorem.** For n>=1 and m>=1,

    (T L_(_,m))_n=S_n^(2m-1)(1)/(2m)!.

Proof: differentiate the rational certificate T R_n(t)=S_n(t+1)-S_n(t)
2m-1 times, sum from t=1 to M, and multiply by -1/(2m)!. The upper boundary
derivative tends to zero because S_n(t) is a constant plus O(1/t) at infinity.
All derivative series converge absolutely. Only the lower boundary remains.
QED.

For m=1 the right side vanishes: R_n has a double zero at t=1. The same
argument fails for m=2 because a third derivative detects the first asymmetric
coefficient around that zero. This identifies the structural role of exactly
one derivative in the original proof.

## An explicit obstruction for zeta(5)

**Theorem.** For n>=1,

    (T L_(_,2))_n
      =-(2n+1)[3+8n(n+1)(H_(n-1)+H_(n+1))]/[n²(n+1)²].    (1)

In particular this sequence of zeta(5) linear forms cannot be O(theta^n) for
any fixed 0<theta<1.

Proof: at t=1, set x=t-1. The first two nonzero coefficients are

    R_n(1+x)=C_n x²(1+d_n x+O(x²)),
    C_n=1/[n²(n+1)²],
    d_n=-2(H_(n-1)+H_(n+1)).

The constant follows from factorial products; the logarithmic derivative of
R_n(t)/(t-1)² gives d_n. Also

    s_n(1)=16n(n+1)(2n+1),
    s_n'(1)=-12(2n+1).

Hence S_n'''(1)/24=C_n[s_n'(1)+s_n(1)d_n]/4, giving (1).
All terms in its bracket are positive. Its magnitude is at least
3(2n+1)/[n²(n+1)²], which decreases polynomially, not exponentially.
If |L_(n,2)|<=C theta^n eventually, its three-term recurrence combination
would be O(n³ theta^n), contradicting that lower bound. QED.

More precisely H_n=log n+O(1) yields a forcing term asymptotic to
-32 log n/n. No asymptotic estimate is needed for the impossibility proof.

## What this teaches us to change

Increasing derivative order fixes the zeta value but destroys the boundary
cancellation. Increasing zero multiplicity could restore cancellation, but then
changes the pole structure, recurrence, and denominator costs. This is a
concrete design tradeoff rather than a vague claim that higher zeta is harder.

One can subtract a particular solution of the forced recurrence to recover a
homogeneous sequence. But that does not by itself prove irrationality: the
correction must have suitable rational arithmetic and must preserve the desired
zeta coefficient and limit. An arbitrary real particular solution can conceal
exactly the unknown constant. A valid next construction must account for that
correction explicitly.

This also clarifies why a determinant method is a substantial change of
viewpoint: it need not preserve this second-order recurrence or its boundary
conditions. See odd-moment-generalization.md for the analytic construction and
odd-arithmetic-barrier.md for remaining arithmetic questions.
