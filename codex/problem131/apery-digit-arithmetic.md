# Apéry numbers, carries, digit products, and algebraic series modulo p

Status: independent complete elementary proofs; not a historical novelty claim. This bridges the two project strands through base-p arithmetic. It does not purport to settle any requirement of an unverified official problem statement. No Claude notes used.

## 1. A general digit-factorization theorem

For positive integers a,b define the integer polynomial

B_n^{a,b}(T)=sum_{k=0}^n binom(n,k)^a binom(n+k,k)^b T^k.

**Theorem.** For every prime p, n>=0, and 0<=r<p,

B_{pn+r}^{a,b}(T) = B_n^{a,b}(T^p) B_r^{a,b}(T) in F_p[T].

Consequently, for a fixed integer t, the sequence B_n^{a,b}(t) has the Lucas property

B_{pn+r}^{a,b}(t) congruent B_n^{a,b}(t) B_r^{a,b}(t) (mod p).

This includes the zeta(3) Apéry numbers A_n=B_n^{2,2}(1), the zeta(2) binomial-sum sequence B_n^{2,1}(1), and central Delannoy numbers B_n^{1,1}(1). The theorem concerns these explicit sums; any recurrence identifications should be checked separately.

### Elementary Lucas lemma

In F_p[X], (1+X)^{pu+v}=(1+X^p)^u(1+X)^v for 0<=v<p, because each middle coefficient of (1+X)^p is divisible by p. Comparing coefficients of X^{pj+s}, 0<=s<p, gives

binom(pu+v,pj+s) congruent binom(u,j)binom(v,s) (mod p),

where a binomial coefficient is zero when its lower index is outside [0,upper]. This proves the only external-looking lemma needed below.

### Proof of the theorem, including the carry case

Write k=pj+s with 0<=s<p and extend the sum by zero terms to 0<=j<=n, 0<=s<p. Lucas gives

binom(pn+r,pj+s) congruent binom(n,j)binom(r,s).

If r+s<p, Lucas applied to the second binomial yields

binom(p(n+j)+(r+s),pj+s) congruent binom(n+j,j)binom(r+s,s).

The product therefore factors into its j-part and s-part, with T^{pj+s}=(T^p)^j T^s.

If r+s>=p, the upper number in the second binomial instead has low digit r+s-p. This is strictly smaller than s, since r<p, so that binomial is zero modulo p. The corresponding small-digit factor binom(r+s,s) is also zero modulo p by the same lemma (its upper number is p+(r+s-p)). Thus this case contributes zero on both sides. This is exactly where a carry must be handled: replacing n+k by digitwise n+k without this check would be an invalid proof.

Terms with s>r are already killed by the first binomial. Summing the factorization now separates j and s and proves the polynomial identity. For integer t, t^p=t in F_p proves the numeric assertion. QED.

## 1b. Reflection of the digit polynomial

**Reflection theorem.** For 0<=r<p,

B_{p-1-r}^{a,b}(T) = B_r^{b,a}((-1)^{a+b}T) in F_p[T].

For 0<=k<p, division by k! is legal modulo p. Writing each binomial as a product of k consecutive factors gives

binom(p-1-r,k)=(-1)^k binom(r+k,k),
binom(p-1-r+k,k)=(-1)^k binom(r,k) (mod p).

These identities hold also when one side is zero. Extend both polynomial sums to k=0,...,p-1 using zero binomial coefficients, raise the identities to powers a,b, and multiply: their roles swap and the sign becomes (-1)^{(a+b)k}. Summation proves the claim. QED.

For Apéry numbers a=b=2, this yields A_{p-1-r}=A_r mod p. In particular A_{p-1}=A_0=1 and P_p(X) is reciprocal of degree p-1:

X^{p-1}P_p(1/X)=P_p(X).

Thus the forbidden digits for p-unit indices occur in reflected pairs r,p-1-r, with the possible central digit (p-1)/2 unpaired for odd p. This explains the symmetric initial-residue tables without invoking the Apéry recurrence. It is a directly proved structural constraint, not a claim of novelty.

## 2. Apéry consequences

Let n=sum_i n_i p^i. Iterating the theorem gives

A_n congruent product_i A_{n_i} (mod p).

**Unit-digit criterion.** Let D_p={r:0<=r<p and A_r is nonzero modulo p}. Then p does not divide A_n iff every base-p digit of n belongs to D_p. In particular a nonzero binomial summand need not imply A_n is a unit: cancellation among the surviving summands remains possible. The criterion tests the small digit *sums*, not just the carries of one summand.

This is an exact finite-state description: store the product of the digit residues in F_p (or only whether it has become zero). Unlike Fuss–Catalan valuations, no multiplication carry state is required after summing the Apéry kernel.

Examples proved by computing only the indicated initial p terms:

* p=2: (A_0,A_1)=(1,5)=(1,1) mod 2, so every A_n is odd.
* p=3: initial residues (1,2,1) give A_n=(-1)^{sum digits}=(-1)^n mod 3, since 3 is odd.
* p=5: initial residues (1,0,3,0,1). Thus A_n is a 5-unit iff every base-5 digit belongs to {0,2,4}. Among 0<=n<5^h exactly 3^h indices are 5-units. Their natural density is zero: if 5^{h-1}<=N<5^h, the proportion below N is at most 3^h/5^{h-1}=5(3/5)^h.

More generally the exact number of p-unit indices below p^h is |D_p|^h. This is an elementary digit-fractal connection, with no need for a theorem about automatic sequences.

## 3. Algebraic generating function after reduction modulo p

Let F(X)=sum_{n>=0} A_n X^n in F_p[[X]], and let P_p(X)=sum_{r=0}^{p-1} A_r X^r. Splitting indices pn+r gives

F(X)=P_p(X) F(X^p).

The Frobenius identity for formal series gives F(X^p)=F(X)^p. Since F(0)=P_p(0)=1, division is legitimate, and

**P_p(X) F(X)^{p-1}=1.**

Thus F is algebraic over F_p(X), of degree at most p-1. This conclusion requires only the explicit finite binomial sum, even though over characteristic zero the recurrence and differential equation carry much richer information. This is a useful contrast of arithmetic regimes, not a claim that characteristic-zero F is algebraic.

Iteration also gives the well-defined formal product

F(X)=product_{j>=0} P_p(X^{p^j}).

Indeed the ignored remainder is F(X^{p^h})=1+O(X^{p^h}), so every coefficient stabilizes. A solution with constant term 1 is unique, by the digit functional equation, or because p-1 is invertible and the coefficient of each successive unknown in F^{p-1} is p-1.

For p=2 this gives F=1/(1+X). For p=3 it again gives F=1/(1+X), agreeing with the residue examples. For p=5, P_5=1+3X^2+X^4=(1-X^2)^2, hence F^4=(1-X^2)^{-2}; these equations and constant term imply F^2=(1-X^2)^{-1} (the alternate square root has constant term -1).

All of this extends verbatim to a numeric specialization B_n^{a,b}(t), replacing the digit polynomial accordingly.

## 4. Bridge to the Fuss–Catalan project

For C_m(n)=binom((m+1)n,n)/(mn+1), let alpha=v_p(mn+1). The proof in ../problem137/valuation.md gives the valuation as all carries in n+mn except its alpha forced initial carries. There, division removes some carries. In the Apéry kernel, a carry in n+k makes binom(n+k,k) vanish modulo p. There, summation kills all cross-digit carry terms and the remaining terms separate into digit sums. Both are consequences of the same base-p factorial/binomial arithmetic, but the mechanism after that arithmetic differs: quotient cancellation versus sum factorization.

A project-sized extension is to classify the small digit polynomials for the family B^{a,b}(t), then compare which specializations have a transparent recurrence or integral interpretation. The theorem above supplies a complete proof framework; conjectures about additional congruences modulo p^2 or p^3 are not implied and need separate work. Likewise digit criteria alone do not construct the sought tree action for Fuss–Catalan exact valuations.

A particularly concrete secondary theorem is now proved in ../problem137/prime-power-arity.md: for full q-ary trees with q=p^e, the count binom(qn,n)/((q-1)n+1) is a p-unit exactly when n=(q^h-1)/(q-1), a base-q repunit (including n=0). Indeed its valuation is s_p((q-1)n)/(p-1)-v_p((q-1)n+1); deleting the trailing p-1 digits leaves valuation equal to the digit sum of the remaining tail divided by p-1. It vanishes exactly when that tail is empty. This yields a striking comparison of digit languages: Apéry units permit arbitrary strings over an allowed digit alphabet, whereas prime-power-arity Fuss–Catalan units permit only repunits. Neither statement asserts a new historical result.

## 5. Audit of the Fuss–Catalan finite-state transducer

The code ../problem137/check_valuations.py maintains a multiplication carry u, an addition carry c, and a Boolean t, initially (0,0,true). Before processing digit j, u is the carry into digit j in multiplication m*n, c is the carry into digit j in n+mn, and t says all previously produced digits of mn equal p-1. For input digit d, the updates are

e=(m*d+u) mod p,
u'=floor((m*d+u)/p),
c'=floor((d+e+c)/p),
t'=t and (e=p-1).

These are precisely the elementary multiplication/addition identities, proving the invariant by induction from j=0. If 0<=u<m and 0<=d<p, then m*d+u<=mp-1, so 0<=u'<m; also 0<=c'<=1. Thus at most 4m states suffice for fixed m,p.

The output weight is c'-1_{t'}. If alpha=v_p(mn+1), then t' is true exactly on j=0,...,alpha-1. The forced-carry proof shows c'=1 at each such position. Therefore each reachable output is nonnegative, and the total output is total carries minus alpha, exactly v_p(C_m(n)). This proves the transducer for every integer n>=0, rather than inferring it from finite checks.

After all input digits are consumed, append zero digits. Each step replaces u by floor(u/p), so eventually u=0. Then a remaining c=1 is cleared on the next zero digit since 1<p, and t becomes false once e=0. Hence flushing terminates. Leading zero digits do not change the answer after the terminal state (0,0,false), which has a zero-output loop. For n=0, the first zero step immediately reaches that state and outputs zero. Reading ordinary least-significant-first digit strings with a finite terminal flush defines a finite automaton for p-unit indices by permitting only zero-output paths. Reversing digit languages also preserves regularity, if a most-significant-first statement is desired; that last standard closure fact is not needed for the computational algorithm.

## References and limits of source verification

The complete proofs above do not rely on an unread modern draft. Primary historical sources:

* Édouard Lucas, *Théorie des fonctions numériques simplement périodiques*, American Journal of Mathematics 1 (1878), pp. 184–240 and 289–321. Publisher-archive first installment: https://www.jstor.org/stable/2369308 . A published translation by Sidney Kravitz, edited Douglas Lind (Fibonacci Association, 1969), is hosted at https://www.fq.math.ca/simply-periodic.html and https://www.fq.math.ca/Books/Complete/simply-periodic.pdf . The scan was retrieved locally. Its title/introduction pages identify this as a translation of the first part (1878, pp. 184–240). Section XXI, translated pp. 60–62 (PDF pages 64–66), was OCR checked: it gives the Pascal-triangle congruence, the two-digit binomial factorization, and iteration to all digits. OCR excerpts are retained under ../sources/lucas-ocr/; OCR equations are imperfect, so the proof above is supplied independently rather than copied from OCR. The downloaded published translation is ../sources/lucas-simply-periodic-translation.pdf.
* E. E. Kummer, *Über die Ergänzungssätze zu den allgemeinen Reciprocitätsgesetzen*, J. reine angew. Math. 44 (1852), 93–146, https://doi.org/10.1515/crll.1852.44.93 . Metadata and publisher link verified through the national-library record https://ndlsearch.ndl.go.jp/books/R100000136-I1361699994682421120 ; full original text has not yet been read in this session. Our elementary factorial proof supplies the carry theorem actually used.

No claim that the Apéry Lucas congruence, the general binomial-family version, or the algebraic-series consequence is new. A novelty assessment would require targeted existing literature beyond these primary historical anchors.
