# Pairing colliding negative squares removes the upper prime band

Status: complete proof, accepted independent review in rank-one-prime-band-independent-review.md, and independent exact checkers. This explains the prime-support cancellation observed in the balanced trial experiment. The proper-rational-function hypothesis is essential and is stated explicitly. No asymptotic denominator or irrationality claim.

Let s=2m+1>=3, K>=m+1, and D_K(t)=product_{j=1}^K(t+j^2). For an integer polynomial F with deg F<K, write

mu_s(F/D_K)=a_F zeta(s)+b_F.

**Theorem.** For every prime p>K, both a_F and b_F are p-integral. In particular, for a **monic degree-d** trial p_trial=P/T with **K=2d+m+1**, its degree-drop rational approximant can have a denominator prime above K only if that prime divides T; its exponent there is at most2v_p(T). The exact zeta coefficient is then (-1)^m, so no further division by an arbitrary affine coefficient is involved.

The theorem applies to every proper integer numerator F, not just a square or a balanced polynomial. The rank-one choice K=2d+m+1 satisfies all its hypotheses.

## 1. The rational pole values are equal at paired residues

For 1<=j<=K put

alpha_j=j^{s-1},
beta_j=-j^{s-1}H_j^{(s)}-1/(s-1)+1/(2j),

so mu_s((t+j^2)^{-1})=alpha_j zeta(s)+beta_j. If p>K, the harmonic denominators and j are p-units. Also p does not divide s-1=2m: p>K>m and p is odd since K>=2. Thus beta_j is p-integral.

For p>K but p<=2K, a collision between nodes -j^2 and -k^2 modulo p occurs exactly when j+k=p. The nodes within1,...,K are otherwise distinct modulo p; each collision is a pair, not a larger cluster. Write k=p-j.

Since s is odd, pairing a with p-a gives H_{p-1}^{(s)}=0 modulo p. Therefore

H_{p-j}^{(s)}=-sum_{a=1}^{j-1}(p-a)^{-s}
=H_{j-1}^{(s)}=H_j^{(s)}-j^{-s} modulo p.

Since s-1 is even, alpha_{p-j}=alpha_j modulo p. Substituting the harmonic reflection into beta yields

beta_{p-j}
=-j^{s-1}(H_j^{(s)}-j^{-s})-1/(s-1)-1/(2j)
=beta_j modulo p.

The correction +1/(2j) is crucial; omitting it would spoil the reflection.

## 2. The barycentric p-residues are opposite

Set z_j=-j^2 and c_j=F(z_j)/D_K'(z_j). Properness gives F/D_K=sum_j c_j/(t+j^2). Since

D_K'(z_j)=(-1)^{j-1}(K-j)!(K+j)!/(2j^2),

a coefficient c_j has at worst a simple p denominator: K-j<p, K+j<2p, and j<p. It has such a possible denominator exactly in the colliding range j>=p-K.

For a colliding pair j,k=p-j, factor

D_K'(z_j)=(z_j-z_k)B_j,
D_K'(z_k)=-(z_j-z_k)B_k,

where B_j and B_k are products over all other nodes. Those products are p-units, and B_j=B_k modulo p because z_j=z_k modulo p. Also F(z_j)=F(z_k) modulo p and z_j-z_k=p(p-2j). Consequently

p c_j=-p c_k modulo p.

Because alpha and beta have equal residues at the pair, both

c_j alpha_j+c_k alpha_k,
c_j beta_j+c_k beta_k

are p-integral: multiply by p, use the cancellation, and recall there was at most one p in the denominator. Noncolliding coefficients are already p-integral. Summing proves the theorem for K<p<=2K. For p>2K there are no factorial denominators involving p, so the conclusion is immediate.

In fact a_F is an integer globally: its barycentric sum is (-1)^m times the leading coefficient of the remainder of t^m F(t) modulo the monic integer polynomial D_K. This stronger observation is optional; the pair argument treats both affine coefficients uniformly.

## 3. Why arbitrary improper numerators are excluded

Without deg F<K the polynomial quotient contributes Bernoulli moments and can reintroduce primes above K. An explicit counterexample to an unrestricted 'any integer P' assertion is

s=3, K=2, P=t^2,
P^2/D_2=t^2-5t+21+(-85t-84)/D_2.

The quotient moment is59/12. The proper part is -341 zeta(3)+405, so

mu_3(t^4/D_2)=4919/12-341 zeta(3).

Its rational coefficient has a3-denominator, with3>K and3<=2K. The proper-part pairing works exactly as proved, but cannot remove a denominator from the polynomial quotient. The balanced rank-one application is proper and is unaffected by this counterexample.

## 4. Refine the exact fixed-divisor normalizer

Use the notation of balanced-rank-one-trial-experiment.md: p_trial=P/T is monic of degree d, P is primitive integer, G=gcd_{j=0}^{2d}P(-j^2), Lambda_K=lcm(1,...,K), K=2d+m+1. The previously proved clearing integer was

C_grid=T^2*m*(2K)!*Lambda_K^{2m+1}/G^2.

Let Pi_K=product_{K<p<=2K, p prime}p. Then the sharper clearing integer is

**C_band=C_grid/Pi_K.**

It is integral and C_band*q_bal is an integer. Thus the entire upper primorial band can be removed, even when a band prime divides T.

Here is the precise arithmetic justification. The primitive integer polynomial f(j)=P(-j^2) has degree2d and coefficient content one. Its Newton expansion has every divided-difference numerator divisible by G; multiplying by(2d)! shows G divides(2d)! (compare coefficient contents). Hence G has no prime factor above2d<K. Each prime in Pi_K occurs exactly once in(2K)!, so C_grid has band-prime valuation1+2v_p(T). Dividing Pi_K leaves valuation2v_p(T), sufficient by the theorem applied to the integer numerator P^2 and subsequent division by T^2. At primes outside the band no normalizer factor is removed. This proves both integrality and clearing.

The same argument proves the exact support statement: no new denominator prime above K comes from the proper pole basis. Any remaining such prime is inherited from the polynomial coefficient denominator T, to at most twice its exponent. In all18 computed balanced cases T itself is supported on primes<=K, so the observed support of the primitive B is now explained by a theorem rather than an empirical cancellation guess. No all-index prime-support theorem for the balanced coefficient denominator T is asserted.

This removes a proved large factor from a sufficient normalizer; it does not establish that the resulting bound is minimal or arithmetically effective for irrationality.

## 5. What the fixed-divisor data do and do not suggest

The exact finite ratios G/(n!)^{m+1} for n=1,...,6 are

m=1: 1,1/2,1,1/2,1/6,1/12;
m=2: 1,1/2,2,1/3,3/2,1/4;
m=3: 1,1/4,4/9,1/18,36/5,1/10.

They justify probing a factorial normalization near this scale, but show no simple uniform divisibility relation: ratios may exceed one or have different small-prime denominators. No exponent-growth or asymptotic statement is inferred. The exact Newton fixed divisors and prime valuations remain available in the experiment JSON for testing any proposed formula.

For m=1 the existing Wilson(1,1,1,0) identification in one-pole-schur-structure.md supplies the all-degree classical hypergeometric formula (equivalently permute parameters to(0,1,1,1)). It proves the integer-grid binomial sum recorded in the experiment note for every n, so the earlier m=1 finite recognition can be upgraded using that established theorem. The m>1 evaluation pattern remains conjectural.

## Reproduce

check_rank_one_prime_band.py independently tests the harmonic reflection and pair cancellation for all monomial proper numerators in a bounded range, includes the improper counterexample, and recomputes the refined normalizers/support bounds for the18 exact experiment rows. The main experiment generator now records Pi_K and C_band as well, so rerunning it preserves the added data.

Independent-review follow-up: rank-one-prime-band-independent-review.md accepts the proper-numerator theorem and both normalizer arguments, including band primes dividing T. Its separate checker tests39 arbitrary rational trials, including27 with genuinely inherited above-K primes, and records a proper-affine-root counterexample illustrating why slope±1 is needed for the root-support conclusion.
