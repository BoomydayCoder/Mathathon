# Odd zeta values: symmetry, failed extensions, and a changed research baseline

Updated 2026-11-13 Los Angeles time. Status: the lemmas below have complete independently derived proofs; their novelty is not claimed. Source reports are explicitly separated from local verification.

## 1. Current-source correction

The earlier project prompt called irrationality of ζ(5) a major open problem. That is no longer a safe description of the current literature. The following primary sources were retrieved in this session:

* Aabir Fauzan, *ζ(5) is irrational*, dated 17 September 2026, [Zenodo record](https://zenodo.org/records/22826419), DOI 10.5281/zenodo.22826419. Local PDF/text: `codex/sources/fauzan-zeta5.*`.
* Dan Romik's [Lean formalization](https://github.com/danromik/zeta5-irrationality), pinned commit `fc2fefa77076e22e7e80b4c101ded5a98d53e14d` (committer date 2026-09-25T09:08:54Z). Four pinned source files are saved under `codex/sources/romik-*`.
* Fischler–Lai–Sprang–Zudilin, [*Many more odd zeta values are irrational*](https://arxiv.org/abs/2611.16156), submitted 12 November 2026. PDF/text saved as `codex/sources/fischler-lai-sprang-zudilin-2026.*`.

Fauzan's Theorem 2.1 supplies integer polynomials of degree 37n whose nonzero positive values at ζ(5) are below exp(-139n²/5). Section 1.1 deduces irrationality by clearing a rational argument's denominator. The construction uses normalized Hankel determinants, positivity, real decay, and local arithmetic estimates. This is a primary-source theorem claim, not a proof independently audited here.

The pinned formalization's `cert/Check.lean` spells out the conclusion directly with Mathlib's `Irrational` and the series sum of (v+1)^(-5). `Zeta5/Axioms.lean` declares the standard prime number theorem θ(x)/x→1. Its README reports no `sorry`, the usual Lean logical axioms plus that PNT axiom, and successful builds. **We inspected these source statements; we did not run Lean or independently reproduce the asserted axiom output.** Neither `lean` nor `lake` was installed. A README report is not a local verification result.

The November follow-up explicitly uses Fauzan's determinant method and states that at least 2.78 sqrt(s log log s/log s) values among ζ(3),ζ(5),...,ζ(s) are irrational for sufficiently large odd s. This is substantial evidence that the preprint is being used in current expert research. The follow-up was retrieved, not fully checked here.

The project should therefore seek a **simpler explanation or proof than the currently reported ζ(5) construction**, and separately investigate all higher odd values. Do not claim a first ζ(5) proof, or call the existing result false merely because we have not checked it.

## 2. A general symmetry lemma (proved)

Let R(t) be a rational function over Q, with poles only at 0,-1,...,-n, pole order at most a, and R(t)=O(t^-2). Write its unique proper partial fractions as

    R(t)=sum_{k=0}^n sum_{j=1}^a c[j,k]/(t+k)^j.

Assume R(-n-t)=-R(t). Then

    sum_{t=1}^infinity R(t)
      = sum_{odd j, 3<=j<=a} C[j] ζ(j) - B,
    C[j]=sum_k c[j,k],
    B=sum_k sum_j c[j,k] H_k^(j) in Q.

Here H_0^(j)=0 and H_k^(j)=sum_{v=1}^k v^(-j), including j=1.

**Proof.** Substituting -n-t into the partial fractions and using uniqueness gives
c[j,n-k]=(-1)^(j+1)c[j,k]. Thus C[j]=0 for even j. Decay O(t^-2) also gives C[1]=0 by comparing the coefficient of 1/t at infinity. Sum over t=1,...,T before taking limits. For j>=2, the shifted sums tend to ζ(j)-H_k^(j). For j=1, the common divergent harmonic term cancels because C[1]=0, and H_(T+k)-H_T→0. The remaining finite expression is exactly the claimed identity. QED.

**Why useful.** Reflection is a projection onto odd pole orders. It removes all even zeta values at once, but it does not distinguish ζ(3) from ζ(5). This is the algebraic content behind the well-poised symmetry relevant here; the full hypergeometric terminology is not needed to prove this lemma.

## 3. Exact mixed-value obstruction in the simplest extension (proved)

Increasing the denominator power of Ball's n=1 summand leads naturally to

    R(t)=(2t+1)(t-1)(t+2) / [t^6(t+1)^6].

It is antisymmetric under t→-1-t. Its partial fractions are

    -2/t^6 +9/t^5 -21/t^4 +33/t^3 -33/t^2
    +2/(t+1)^6 +9/(t+1)^5 +21/(t+1)^4
    +33/(t+1)^3 +33/(t+1)^2.

Consequently

    sum_{t>=1} R(t)=18 ζ(5)+66 ζ(3)-98.

All equalities can be checked by clearing the common denominator. This proves a concrete limitation of the naive change from fourth-order to sixth-order poles: it creates a nonzero ζ(3) coefficient. A small nonzero integer linear form in 1,ζ(3),ζ(5) does not establish irrationality of ζ(5), because ζ(3) is already irrational.

This example is not a no-go theorem for all constructions. One could try eliminating ζ(3) between several forms, but then must independently recheck coefficient growth, denominators, residual decay, and nonvanishing. Cancellation of unwanted values can destroy the analytic gain.

## 4. A different extension really does isolate each odd value (proved), but loses the easy decay argument

Use Zudilin's

    R_n(t)=[product_{j=1}^n(t-j)/product_{j=0}^n(t+j)]²
          =sum_k [a_k/(t+k)²+b_k/(t+k)].

For m>=1 define

    L[n,m]= -1/(2m-1)! sum_{t>=1} R_n^(2m-1)(t).

Termwise differentiation of the finite partial fractions and absolute convergence give

    L[n,m]=2m A_n ζ(2m+1)-B[n,m],
    A_n=sum_k a_k,
    B[n,m]=sum_k [2m a_k H_k^(2m+1)+b_k H_k^(2m)].

The unwanted ζ(2m) vanishes because sum_k b_k=0. If D_n=lcm(1,...,n), then

    A_n is integral and D_n^(2m+1) B[n,m] is integral.

For completeness, write the unsquared rational function as sum_k d_k/(t+k), where d_k are signed integers given by binomial coefficients. Squaring it gives a_k=d_k² and
b_k=2 d_k sum_{l!=k} d_l/(l-k), hence D_n b_k is integral. Multiplication by D_n^(2m+1) clears every harmonic denominator in B[n,m]. This proves the arithmetic assertion without denominators up to 2n.

Why is this not already an irrationality proof for all odd values? The first-derivative case benefits from double zeros: R_n'(1),...,R_n'(n) vanish. For m>=2 the higher derivatives need not vanish. For example

    L[1,2]=20ζ(5)-20,
    -R_n'''(1)/6 = 2(H_(n-1)+H_(n+1))/(n²(n+1)²), n>=1.

The latter identity follows by writing R_n(t)=(t-1)²h(t), evaluating h(1)=1/[n²(n+1)²], and logarithmically differentiating h. This surviving term is only polynomially small. **It does not by itself prove a lower bound on the whole sum**, since other terms may cancel. It does prove that the root-based removal of the initial summation segment no longer works. A new global estimate and a nonvanishing proof are required; the old positivity argument cannot simply be copied.

This yields a useful research distinction: symmetry may leave mixed odd values, while differentiation can isolate a single odd value at the price of disrupting the smallness mechanism. Stronger root multiplicities are a concrete next experiment, but they also change pole order, arithmetic cost, and the resulting zeta combination.

## 5. A polynomial irrationality criterion (proved independently)

Let ξ be real. Suppose P_n in Z[X] are nonzero polynomials, with degrees d_n>=1, and

    0 < |P_n(ξ)| <= exp(-a_n),    a_n/d_n → infinity.

Then ξ is irrational. If ξ=u/v with v>0, the number v^d_n P_n(u/v) is a nonzero integer. Its absolute value is at most exp(d_n log v-a_n), which is eventually below 1. Contradiction.

No coefficient-height estimate is needed for this conclusion. Heights become relevant for quantitative approximation to nearby rational arguments, because the small value at ξ alone does not prevent P_n from having a rational root nearby.

The classical linear-form criterion is the degree-one case. The theorem permits degrees growing like n if values decay like exp(-c n²). This makes determinant constructions conceptually natural rather than an arbitrary escalation of complexity.

## 6. Why Hankel determinants connect approximation, positivity, and algebra (proved)

Let μ be a positive measure on R with finite moments through degree 2h-2, and suppose no nonzero polynomial of degree <h vanishes μ-almost everywhere. Let

    G[i,j]=integral x^(i+j) dμ(x),    0<=i,j<h.

Then det G>0. Indeed, for every nonzero real vector v,

    v^T G v = integral (sum_i v_i x^i)² dμ(x)>0.

Thus G is positive definite. Moreover,

    det G = (1/h!) integral ... integral
              product_{i<j}(x_i-x_j)² product_i dμ(x_i).

Expand the two Vandermonde determinants; integrate termwise (finite sums, with absolute integrability following from the moment assumptions); each permutation pairing contributes the Gram determinant, and symmetry gives the factor h!. This is also the squared volume of 1,x,...,x^(h-1) in L²(μ), and via Gram–Schmidt the product of squared norms of the monic orthogonal polynomials.

If the moments are affine expressions u_j ξ+v_j with rational coefficients, det G is the value at ξ of a rational polynomial of degree at most h. The Gram argument supplies nonvanishing without needing to understand its coefficients. The Vandermonde formula supplies a route to decay. Clearing denominators supplies integrality, but it can erase the entire analytic gain. These three independent requirements explain the architecture of determinant irrationality proofs.

This gives an unexpected but exact link: rational approximation → moment matrices → orthogonal polynomials → Vandermonde interaction. It also identifies the crucial remaining task for any proposed generalization: a denominator estimate on the same n² scale as the determinant's real decay.

## 7. Classical family result and bounded next steps

The retrieved primary Rivoal paper [math/0008051](https://arxiv.org/abs/math/0008051), Theorem 1, states: for every ε>0 and all sufficiently large n,

    dim_Q span{1,ζ(3),ζ(5),...,ζ(2n+1)}
       >= (1-ε) log(n)/(1+log 2).

It gives information about the family, not every specified value. The PDF/text are saved in `codex/sources/rivoal-0008051.*`. Its proof and modern strengthenings were not fully rederived here.

For the team's next stage:

1. Verify the pinned formalization in a suitable Lean environment, including `cert/Check.lean` and kernel replay; record the precise commit and axiom output. This is an audit task, not an invitation to assume a README is certification.
2. Explain the independently proved polynomial criterion and Gram determinant identity before opening the large local-arithmetic calculations in Fauzan.
3. Find a smaller illustrative determinant whose moment entries involve one constant and whose arithmetic can be checked by hand. It need not prove ζ(5) irrational to teach the mechanism faithfully.
4. Study stronger root multiplicities in the derivative-isolation family; keep a ledger of all extra zeta coefficients and actual denominator exponents.
5. Treat a proof for every odd value as a research ambition. No such proof, and no simplification of Fauzan sufficient to replace his estimates, is supplied here.
