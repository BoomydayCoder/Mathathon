# Adversarial review of the extra zeta(5) outer-prime factor

**Verdict.** The finite-prime statement in `outer-polynomial-correction-kernel.md` is valid, with the original source local decomposition as an explicit dependency. I found no unresolved gap in the polynomial reduction, full-matrix couplings, old negative-half-weight rows, or parity lift. The improvement is coefficientwise in X. Its prime-averaged saving is27/3200. It is a sharper sufficient arithmetic estimate, not a claim that the determinant's valuation equals this bound or that the entire preprint is thereby independently certified anew.

The research lead came from the parent's read-only consultation of Claude's published summary, as recorded in the construction note. This audit used the actual mathematical construction and source local definitions; no historical-novelty claim is made.

A subsequent simplification is proved in [outer-correction-Hermite-ideal.md](outer-correction-Hermite-ideal.md): an explicit polynomial ideal supplies the entire guaranteed kernel, replacing the target-column count by a degree calculation. It gives the same gain, not a second additive gain. Its exact script checks CRT span membership, independence, and all critical columns.

## 1. Exact statement and baseline

Let K be a positive multiple of40, N=3K/40, h=K-N, and let p satisfy the source's outer hypotheses

    p>=7, 2p<K<3p, p^2>2K, 2N<p, 5N<=2p-2.

Suppose d=K+3N-3p+1>0. In the source's outer CRT basis, write the actual symmetric moment matrix as

    G=A+p^-1 L,

where L is integral and independent of X, and A has the original coefficientwise row bounds v_p^G(A_ij)>=w_i+w_j. The correction acts only on the polynomial quotient of the rational input, not on its simple-pole values.

Let Z be the old zero-weight rows, z=|Z|. Replacing each old zero weight by-1/2 and leaving every negative weight unchanged gives weights w'_i for which

    v_p^G(G_ij)>=w'_i+w'_j.

Indeed the A bounds only weaken, while w'_i+w'_j<=-1 absorbs p^-1 L entrywise. The baseline determinant exponent is2sum w_i-z. In this range the source rank bound r=K+4N-2p+2 exceeds z, so this equals the original refined exponent gamma_out. The proof of the additional gain below does not assume a rank drop over Q_p.

## 2. The Bernoulli correction is an actual finite-field derivative sum

For every p-integral polynomial Q(t),

    p mu(Q) = (1/24)sum_(a in F_p)
                [D_v^3(v^5Q(-v^2))]_(v=a) (mod p).       (1)

For a polynomial P(v), the identity pL(P)=sum_a P(a) modulo p follows termwise from von Staudt--Clausen and the finite-field power sum. For the constant term, both sides are0 modulo p; for degree1 both sides are also0. Applying it to P=D_v^3(v^5Q(-v^2)) proves(1), with the positive sign and factor1/24 shown. There is no Bernoulli B_1 convention ambiguity in this step.

In particular, if Q has a zero of order at least4 at t=-a^2 with a nonzero modulo p, its contribution at v=+/-a to(1) vanishes: substitution t=-v^2 is a local coordinate there. At a=0 the factor v^5 kills the third derivative for every polynomial Q.

For comparison, the exact derivative is

    60v^2 Q(-v^2)-150v^4 Q'(-v^2)
                 +72v^6 Q''(-v^2)-8v^8 Q'''(-v^2).

This shows explicitly why arbitrary unit factors cannot restore a contribution when the polynomial has order at least4.

## 3. Polynomial division is safe despite collided denominator nodes

Write D=D_K/D_N and W=D_N^5. A source CRT row is g_(a,i)=P_a q_(a,i), P_a=D/Q_a. After reduction modulo p,

    Q_a=(t+a^2)^(ell_a-delta_a),
    q_(a,i)=(t+a^2)^i,

including the rows with the special bad-pole factor. Their actual CRT unit factors remain in P_a; they are not being replaced by1.

For a removed class a<=N, W has order5 at its node. Here its denominator multiplicity ell_a-1 is at most5. At every other node g_(a,i) contains the full denominator multiplicity. Thus W g_(a,i) is divisible by D **as a polynomial over F_p**, and

    H_(a,i)=bar(W)bar(g_(a,i))/bar(D)

is an honest polynomial. Since D is monic over Z_p, quotient and remainder in polynomial division reduce correctly even though bar(D) has repeated roots. Therefore the polynomial quotient of W g_c g_j/D reduces to H_c bar(g_j).

This is not a substitution into a rational function at an unresolved pole. The multiple denominator nodes are cancelled by a polynomial identity before(1) is applied. That distinction is essential and is satisfied here.

## 4. Audit every possible coupling from a selected row

Select C consisting of all old zero-weight rows in removed classes:

    ell=5: local orders i=2,3;
    ell=6: local orders i=2,3,4.

For c=(a,i) in C, at its own node H_c has order6+i-ell_a. Thus its minimum order is3 in a removed size-five class and2 in a removed size-six class. At a different removed node it has order at least5. At an unremoved node its guaranteed order is0.

These facts give the complete coupling check:

- At a removed node other than a, H_c already has order5, so every column is harmless there.
- At a, a column from another class contains the full local denominator multiplicity ell_a-1>=4, so the product has order at least4.
- At a, another selected same-class row has order at least2. The product order is at least5 when ell=5 and at least4 when ell=6. Hence L_(C,C)=0 modulo p.
- The only old weight-1/2 rows in this band are local order1 in removed size-five classes. At their own selected node the product order is3+1=4. Hence L_(C,T)=0 modulo p for this entire old weight-1/2 set T.
- At an unremoved node, a column from any other class has the full multiplicity ell>=4. Only a column originating at that node can contribute. Its zero-weight local orders are ell-2,ell-1. Thus both columns can matter for ell=4, only the first can matter for ell=5, and neither can matter for ell=6.
- The zero node contributes nothing because of v^5. A zero-class column contains full multiplicity at all ordinary nodes, so it produces no additional exception.

Every unknown CRT unit and its derivatives are harmless in these vanishings: multiplication by a unit cannot lower a zero's order. At potentially nonzero target columns no unit is discarded; the proof only bounds their number, so it does not need a guessed explicit kernel vector.

There is one especially important distinction. A selected row can have a nonzero polynomial correction against a removed size-five row of local order0, whose weight is-3/2. This does **not** invalidate the proof: such deeper negative-half-weight columns vanish after the baseline scaling, as checked in section6. The old weight-1/2 rows do not vanish automatically, which is why they needed the separate order-four argument above.

## 5. Exact counts, including the finite +1

Put v=K-2p and m=(p-1)/2. An ordinary class a has size

    ell_a=4+1_(a<=v)+1_(a>=p-v).

The condition d>0 implies v>N. For example, 3p<K+3N+1 gives

    v-N=(K-N)-2p>13K/120-2/3>0,

using K>=40. Hence no removed class has size4. The number of removed size-six classes is

    u=(N+v-p+1)_+=(K+N-3p+1)_+.

Thus |C|=2N+u. The potentially nonzero unremoved target columns number

    b=2*(number of unremoved size-four classes)
         +(number of unremoved size-five classes)
      =(3p-K-N-1)_+.

For v<=m this is2(m-v)+(v-N); for v>=m+1 it is(p-v-N-1)_+. These cases are exhaustive since p is odd. Therefore

    |C|-b=K+3N-3p+1=d.                              (2)

Also z=|Z|=p+u, including the one zero-weight row in the zero class. The source rank bound obeys

    r-z=d+N+1-u>=d+1>0,

since u<=N. This confirms the identification with the original exponent in section1.

The reduced matrix L_(C,Z) has at most b nonzero columns. Its left kernel over F_p consequently has dimension at least d. This is an ordinary constant-coefficient kernel, not a kernel requiring rational functions of X: L itself comes solely from the polynomial moments and is independent of X. A basis can be lifted to a unimodular change of coordinates on the C rows by ordinary Gaussian elimination over F_p and unit lifts to Z_p.

## 6. The full scaled matrix, including all negative weights

Put pi=sqrt(p) and B_ij=p^(-w'_i-w'_j)G_ij. This matrix is integral over the ramified valuation ring. The selected C rows have w'=-1/2.

Against an old zero row, B_Cj=pA_Cj+L_Cj, and A_Cj is integral. Against an old weight-1/2 row, the same formula holds, and A_Cj is again integral: its old lower bound is-1/2, but all coefficients lie in Q_p, so their valuations are integers and must be at least0. Therefore in these two cases the reduction of B_Cj is exactly the correction entry checked above.

Against an old weight -k-1/2 with k>=1, scaling gives p^(k+1)A_Cj+p^k L_Cj. The old A bound rounds up to-k, so both terms are divisible by p. This includes the potentially nonzero correction against the old-3/2 rows identified in section4.

Against an old integer negative weight -k, k>=1, both terms have valuation at least1/2 after scaling. These columns lie in the opposite parity sector and need no additional cancellation.

Thus each of the d lifted combinations annihilates the reduction of the **entire same-parity principal block**, not just its selected self block. All opposite-parity entries have the required automatic half valuation. Apply the change simultaneously to rows and columns; symmetry is preserved. Selected-selected entries are divisible by p, selected-other same-parity entries by p, and selected-opposite-parity entries by pi. Dividing each selected row and column by pi leaves an integral matrix and removes a determinant factor p^d.

Consequently

    v_p^G(det G)>=2sum_i w'_i+d=gamma_out+d.           (3)

This argument only uses a kernel **modulo p** followed by an integral lift. It does not claim that the polynomial correction's rank drops over Q_p, and it does not need such a claim. All operations and valuation bounds are coefficientwise in X. No invertibility assumption on a complementary block is used.

## 7. Asymptotic budget and independent checks

With N=3K/40,

    d=p[(49/40)(K/p)-3]+1.

The fixed-band prime-number-theorem conversion gives the saving

    integral_(120/49)^3 [(49/40)x-3]/x^3 dx=27/3200.

The finite +1 and interval endpoints affect only o(K^2) terms. This is relative to the source's refined outer bound. Relative to a plain-rank bound there is at least this much gain, but one must avoid double-counting any separate recovery of the original zero-weight refinement.

I independently checked the residue identity using exact Bernoulli fractions and literal finite-field derivative sums at primes7,11,29,31,47,97 across a range of monomial degrees, and checked the integral exactly. I reran the supplied actual-CRT correction script: all seven cases passed, including all forbidden columns and the old weight-1/2 rows. Those finite checks support, but do not replace, the full-matrix proof in sections3–6.

The valid conclusion is a real extra divisor in the coefficientwise determinant bound. It does not identify the exact determinant valuation, give a new standalone irrationality theorem, or resolve historical novelty. The other source prime regimes and the analytic period estimate remain separate dependencies.
