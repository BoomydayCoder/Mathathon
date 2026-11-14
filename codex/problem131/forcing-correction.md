# Boundary forcing: what a rational correction can and cannot repair

Status: independently verified forcing asymptotic and proved obstruction to finite harmonic corrections. No historical novelty claim. This note audits `forced-recurrence-obstruction.md` and develops the correction question posed there.

Write

    (T u)_n=(n+1)³u_{n+1}-(2n+1)(17n²+17n+5)u_n+n³u_{n-1}.

The higher derivative forms are L_(n,m)=-(2m)!^(-1)sum_{t>=1}R_n^(2m-1)(t).

## 1. Independent audit of the all-m forcing asymptotic

The expansion in the companion note is correct:

    R_n(1+x)=x²/[n²(n+1)²] * exp(sum_{j>=1} ell_(n,j)x^j),
    ell_(n,j)=(2/j)[(-1)^j H_(n+1)^(j)-H_(n-1)^(j)].

After removing the zero, numerator factors are 1-x/l for l=1,...,n-1 and denominator factors are 1+x/l for l=1,...,n+1. Twice their logarithms give the stated signs and factor 2/j. The Taylor series converge for |x|<1, uniformly sufficient for taking any fixed finite coefficient.

Thus ell_(n,1)=-4log n+O(1), while each fixed ell_(n,j), j>=2, stays bounded. For fixed k the exponential coefficient is

    [x^k]exp(sum ell_j x^j)=(-4log n)^k/k!+O((1+log n)^(k-1)).

The multiplier has coefficients

    s_n(1+x)=16n(n+1)(2n+1)-12(2n+1)x-8(2n+1)x².

For k=2m-3, the leading term of [x^(2m-1)]S_n is therefore

    (32/n)*(-4log n)^k/k!.

Since S_n^(2m-1)(1)/(2m)! equals that coefficient divided by 2m, the final result is

    (T L_(_,m))_n
      ~ -16*4^(2m-3)/[m(2m-3)!] * (log n)^(2m-3)/n.

The constant and sign in the companion theorem are correct. Terms arising from the linear and quadratic multiplier coefficients are O(n^-3) times lower powers of log n and cannot alter the leading term. For m=2 this is -32log n/n, consistent with the exact formula below.

## 2. Simplify the zeta(5) forcing before searching for a correction

Let H_n=H_n^(1). Since H_(n-1)+H_(n+1)=2H_n-1/[n(n+1)], the exact forcing becomes

    E_n:=(T L_(_,2))_n
      =-16(2n+1)H_n/[n(n+1)]+5(2n+1)/[n²(n+1)²].       (1)

A natural hope is to find a particular solution P_n=a(n)H_n+b(n), where a,b are rational functions of n with rational coefficients. Such a correction would have understandable arithmetic and could potentially leave a homogeneous, exponentially small remainder.

That hope fails in this precise class.

## 3. A rational-function obstruction

**Lemma.** The equation

    T a(n)=-16(2n+1)/[n(n+1)]                           (2)

has no rational-function solution a(n). Also T a=0 has no nonzero rational-function solution.

**Proof by extreme poles.** Group the finitely many poles of a into classes modulo integer translation. In a class not containing integers, choose the pole z with smallest translate index. At n=z-1, only the term (n+1)³a(n+1) can have a pole. Its coefficient z³ is nonzero, so the pole cannot cancel; the right side has no pole there. Thus every pole must be an integer.

Suppose integer poles exist, with smallest m and largest M. At n=m-1 only the forward term can be singular, unless its coefficient vanishes. Because the right side has poles only at 0,-1, we must have m=0 (possible coefficient zero) or m-1 in {0,-1}. Thus m belongs to {0,1}.

At n=M+1 only the backward term can be singular. Its coefficient vanishes only for M=-1; otherwise the right side can have a pole there only if M+1 is 0 or -1. Thus M belongs to {-1,-2}. This contradicts m<=M. The argument allows coefficient zeros without assuming they completely cancel a pole, so it remains valid for poles of any order.

Hence a has no finite poles and is a polynomial. Then T a is a polynomial, which cannot equal (2)'s nonpolynomial right side.

For the homogeneous equation the same argument allows only smallest pole m=0 and largest pole M=-1, again impossible. A polynomial a of degree d>=0 with leading coefficient c has T a of degree d+3 with leading coefficient -32c, so no nonzero polynomial solves the homogeneous equation. QED.

**Theorem.** There is no particular solution of T P=E of the form

    P_n=sum_{j=0}^d a_j(n) H_n^j,

with finitely many rational functions a_j(n).

**Proof.** First, H_n satisfies no nonzero polynomial relation over rational functions of n when evaluated for all sufficiently large positive integers. To see this, expand the finitely many rational coefficients at infinity, select their largest power of n, and within that group select the highest nonzero power of H_n. Since H_n=log n+gamma+o(1), the resulting term dominates: lower powers of log n cannot cancel it, and lower powers of n multiplied by any fixed log power are negligible.

We may therefore compare powers of H_n after using H_(n+1)=H_n+1/(n+1) and H_(n-1)=H_n-1/n. If d>=2, the coefficient of H_n^d gives T a_d=0, impossible by the lemma. If d=1, the coefficient of H_n gives (2), also impossible. If d=0, no H_n term can arise, whereas (1) has a nonzero H_n coefficient. QED.

This includes any finite linear combination of shifted H_(n+j) with rational-function coefficients, and any finite products of such shifts: each reduces to a polynomial in H_n with rational-function coefficients.

## 4. Stronger conclusion: these finite corrections cannot even leave exponential decay

For a polynomial P_n in H_n with rational-function coefficients, E_n-(T P)_n is again in the same class. The theorem shows it is nonzero. The dominance argument in its proof shows any such nonzero expression has a leading asymptotic c n^d(log n)^j, with c nonzero and integers d,j (j>=0). It cannot decay exponentially.

If L_(n,2)-P_n were O(theta^n) with theta<1, applying T would give E_n-(T P)_n=O(n³theta^n), a contradiction.

**Therefore no fixed finite harmonic-rational correction in this class can turn the original zeta(5) linear forms into exponentially small forms.** This is stronger than merely failing to guess a particular solution.

The restriction is essential: it does not cover growing-length sums, generalized harmonic sums, hypergeometric expressions, or determinant constructions.

## 5. Local corrections improve powers of n but not the exponential rate

The contribution at the first zero t=1 is

    P_n^(1)=-R_n'''(1)/24
      =[H_n-1/(2n(n+1))]/[n²(n+1)²].                  (3)

It captures the leading log n/n^4 term. Direct substitution gives

    E_n-(T P^(1))_n = (32-32H_n)/n^5
                       +(80H_n-80)/n^6
                       +(282-232H_n)/n^7
                       +O((1+H_n)/n^8).

So it removes several leading powers in the forcing but leaves a nonzero polynomial-logarithmic remainder. This exact asymptotic expansion was checked by symbolic rational algebra.

More generally, for fixed positive integer j and n>=j,

    -R_n'''(j)/24
      =[(j-1)!^4 (n-j)!²/(2(n+j)!²)]
          [H_(n-j)+H_(n+j)-2H_(j-1)].                (4)

Proof: write R_n=P_n², differentiate three times at a simple zero of P_n, and use

    P_n'(j)=(-1)^(n-j)(j-1)!²(n-j)!/(n+j)!,
    P_n''(j)/(2P_n'(j))=2H_(j-1)-H_(n-j)-H_(n+j).

For fixed j, the factorial ratio is a rational function of n and the harmonic shifts equal H_n plus rational functions. Thus subtracting the contributions from any fixed number of zeros lies within the obstruction theorem. It can improve algebraic accuracy but cannot produce exponential decay.

Subtracting all n zero contributions is a different operation: the correction has growing length and escapes this theorem. The truncation investigation in `truncated-derivative.md` studies that option. Its arithmetic cost requires a separate argument and is not certified merely by the present analytic identities.

## 6. A rational particular sequence always exists, but its arithmetic is another question

The absence of a rational-function/harmonic-polynomial closed form does not mean no rational sequence solves T P=E. Starting with rational P_0,P_1 and applying the forced recurrence gives rational terms.

There is an exact variation-of-parameters formula. Let A_n be the positive Apéry solution and define

    W_n=A_nP_(n+1)-A_(n+1)P_n,
    V_n=(n+1)³ W_n.

Then T P=E is equivalent to

    V_n-V_(n-1)=A_n E_n.

Consequently

    P_n=A_n [P_0/A_0
          +sum_{j=0}^{n-1} (W_0+sum_{i=1}^j A_i E_i)
                             /((j+1)³ A_j A_(j+1))].   (5)

Every finite term is rational if the initial data are rational. Formula (5) is proved by telescoping first V_n and then P_(n+1)/A_(n+1)-P_n/A_n.

However, this formula introduces denominators involving Apéry numbers as well as powers of consecutive integers. Cancellations may exist, but formula (5) alone does not control them by a usable exponential bound. Setting initial data to cancel the growing component typically involves infinite sums; allowing arbitrary real constants there could simply hide the unknown zeta value. Any proposed irrationality proof must prove the correction's arithmetic independently.

## Verification

`codex/computations/check_forcing_correction.py` checks the exact rewriting (1), the residual expansion from (3), and finite coefficient expansions for the all-m forcing formula. The no-rational-solution theorem is a pole argument, not a finite search. No conclusion about the complete Fauzan preprint or a general impossibility of zeta(5) irrationality proofs follows.

## 7. Independent audit of the truncated-tail asymptotic

Reviewed the completed `truncated-derivative.md`, including its exact moving-boundary recurrence and its new sharp asymptotic

    E_n ~ (pi/2)16^(-n)log n/n³.

No error found in the reviewed argument. In particular:

- The moving boundaries contribute +(n+1)³R_(n+1)'''(n+1) and -n³R_(n-1)'''(n) to 24g_n. Their signs agree with the negative normalization -1/24.
- The root contribution divided by r_n=R_n(n+1) is exactly 3(n+1)h_n-3/2: the root derivative ratio is P_(n+1)'(n+1)/P_n(n+1)=1/(2n+2).
- The logarithmic-derivative shifts d1,d2,d3 stated there follow by removing the last terms from the ordinary and higher harmonic sums. The factor n³R_(n-1)(n)/r_n=4n(2n+1)² is correct.
- After the cancellation, the forcing has leading constant -287pi/32 times 16^(-n)log n. Applying T to f_n=16^(-n)log n/n³ gives leading constant -287/16. Their ratio is pi/2, with the correct positive sign.
- The maximum-principle squeeze is valid. The added epsilon*2^n handles infinity because E_n is polynomially bounded. The correction M*eta^n with alpha<eta<1/16 handles the finite starting boundary without changing the normalized limit. The proof does not assume a sign for E_n in advance.
- The prime-band denominator argument and the mixed-lcm upper bound were checked algebraically. The warning that one prime family does not prove a global exponential denominator lower bound is essential and correctly retained.

Independent 180-digit numerical checks using exact rational coefficients give

    n       E_n / [16^(-n)log n/n³]
    2       0.153451123541758...
    5       1.125663610940199...
    10      1.471060294264214...
    20      1.507413690975363...
    40      1.520534948088725...

The predicted limit is pi/2=1.570796326794896.... Convergence is slow because h_n=log n+gamma-log2+O(1/n), producing a visible 1/log n correction. These finite values are consistency checks; the recurrence squeeze is the proof.

This tail result and the finite-correction obstruction are compatible. Subtracting all n root contributions is a growing-length operation and escapes the rational-H_n-polynomial class. It repairs the analytic side precisely, while exposing new prime denominators.
