# Why this rational square? A classification, an arithmetic optimum, and a deformation that fails

Status: independently derived design theorems and an all-n deformation analysis. Historical novelty has not been assessed. These results classify a precisely stated ansatz, not all irrationality proofs. They supplement rather than repeat the general pole-denominator theorem in `denominator-design.md`.

## 1. The square is forced by a saturated interpolation problem

Fix n>=1. Consider rational functions R(t) satisfying:

1. Their only possible poles are0,-1,...,-n, each of order at most two.
2. R(t)=O(t^-2) at infinity.
3. R(j)=R'(j)=0 for j=1,...,n.

**Classification theorem.** Every such function has the form

    R(t)=C [product_(j=1)^n(t-j)/product_(j=0)^n(t+j)]².       (1)

In particular, imposing lim_(t→infinity)t²R(t)=1 gives exactly Zudilin's R_n. Squaring has not been assumed; it follows from the constraints.

**Proof.** Multiply R by D(t)², where D(t)=product_(j=0)^n(t+j). The resulting numerator is a polynomial of degree at most2n, because D² has degree2n+2 and R decays at least quadratically. The n double zeros force divisibility by product_(j=1)^n(t-j)², already of degree2n. Hence the numerator is a constant multiple of that product. QED.

This also proves that every allowed pole is genuinely double for a nonzero solution: none of the prescribed positive zeros can cancel a nonpositive pole. The solution space is one-dimensional, and positivity of the normalization C makes R nonnegative on the real line away from poles.

**Minimal pole budget.** With poles of order at most two, n distinct prescribed double zeros, and decay O(t^-2), at least n+1 distinct poles are necessary. If there were M poles, the numerator degree would be at most2M-2, while the zero conditions require degree at least2n. Thus M>=n+1. The construction uses exactly that budget.

These are conditional minimality statements. Requiring only R'(j)=0, without R(j)=0, is a different interpolation problem and is not covered. Likewise a proof could use higher poles, other contours, or a different summation mechanism.

## 2. Integer spacing minimizes the universal pole-difference denominator

For n+1 distinct integers b_0,...,b_n define

    L(b)=lcm{|b_i-b_j|:i<j}.

**Theorem.** D_n=lcm(1,...,n) divides L(b). The consecutive set0,1,...,n attains L(b)=D_n.

**Proof.** For every prime power p^e<=n, there are more than p^e integers in the set. Two are congruent modulo p^e, so one nonzero difference is divisible by p^e. Thus L(b) contains every prime power required by D_n. Consecutive integers have exactly the difference set1,...,n, proving attainment. QED.

This is the relevant arithmetic quantity because, when

    P(t)=sum_k c_k/(t+b_k),  c_k integral,

squaring gives simple-pole coefficients

    b_k'=2c_k sum_(l!=k)c_l/(b_l-b_k),

whose denominators are cleared by L(b). Therefore the consecutive lattice minimizes this *universal differences-based multiplier in divisibility order*. It does not prove optimality of the actual reduced denominators for every specially chosen residue vector; cancellations may improve those.

Equality in the theorem does not force consecutiveness. For example the five poles with b-set{0,2,3,4,6} have L=12=D_4. This example also exposes the other arithmetic costs. With zeros1,2,3,4 and monic numerator, the corresponding simple-pole residues are

    1/6, -45/2, 280/3, -105, 35.

They need scaling by6 to become integral, unlike the standard consecutive construction. Moreover the shifted harmonic sums reach6 rather than4, so D_6=60 enters a naive constant-term bound. Optimizing pairwise differences alone does not optimize the whole proof.

Among n+1 distinct nonnegative integer pole positions with minimum0, the smallest possible span is n, uniquely attained by consecutive positions. Thus the standard choice simultaneously minimizes span and the universal difference multiplier, while producing integer residues without scaling.

## 3. Why one derivative kills exactly the unwanted zeta value

For any rational R with these grid double poles and O(t^-2), write

    R(t)=sum_k[a_k/(t+k)²+b_k/(t+k)].

Decay at infinity gives sum b_k=0. Hence

    -sum_(t>=1)R'(t)
       =2(sum a_k)ζ(3)
          -sum_k[2a_kH_k^(3)+b_kH_k^(2)].              (2)

There are only two possible zeta values after one derivative: ζ(3) from the double poles and ζ(2) from the simple poles. One scalar relation, sum b_k=0, eliminates the latter. Thus the lowest nontrivial pole multiplicity is exactly enough for one cancellation.

For the classified function, R'(1),...,R'(n) all vanish, since the zeros are double. The same multiplicity also makes the lower boundary S_n'(1) vanish when the short certificate S_n=s_nR_n is differentiated and summed. These are separate benefits: deleting the first n summands and eliminating the recurrence's boundary term.

With the same pole budget and quadratic decay, requiring zeros of order at least r+1 at those n points is impossible for any r>=2: it would require numerator degree at least(r+1)n>2n. Thus the *multiplicity-based* way of preserving all initial zeros is compatible with one derivative and incompatible with higher derivatives in this minimal double-pole ansatz. This does not rule out accidental derivative cancellations without high-multiplicity zeros, or other rational families.

The broader tradeoff theorem in `denominator-design.md` explains why higher poles can restore multiplicity while introducing extra zeta coefficients and larger denominator exponents.

## 4. Which reflection is actually present?

The original R_n does **not** satisfy R_n(-n-t)=R_n(t). For n=1, reflection changes its numerator(t-1)² to(t+2)² while leaving the denominator t²(t+1)² unchanged.

In fact, imposing that reflection in addition to the classification hypotheses would force the zero function for n>=1. Every positive double zero j would require another double zero -n-j. These2n distinct zeros would give numerator degree at least4n, exceeding the allowed2n.

The actual symmetry is reciprocal. Put

    P_n(t)=1/t product_(j=1)^n (t-j)/(t+j),
    B_n(t)=tP_n(t).

Then

    B_n(-t)=1/B_n(t),
    R_n(t)R_n(-t)=t^-4.                              (3)

On the imaginary axis, every factor(iy-j)/(iy+j) has modulus one. In the open right half-plane it has modulus strictly less than one, because

    |z+j|²-|z-j|²=4j Re z>0.

Thus B_n is a finite rational inner function for a half-plane, and

    |P_n(z)|<=1/|z|  when Re z>0.

This is a concrete connection between the auxiliary rational function and elementary complex interpolation. It explains the balanced positive-root/negative-pole product and supplies useful size bounds. It does not by itself give the delicate exponentially small derivative sum.

A genuinely additive reflection R(-n-t)=R(t) is a possible different family: its double-pole coefficients satisfy a_(n-k)=a_k and its simple-pole coefficients satisfy b_(n-k)=-b_k, so sum b=0 automatically. Such a family cannot simultaneously carry all the original positive double zeros with the same pole budget. Symmetry is therefore an optional alternative constraint, not the explanation for ζ(2) cancellation in Zudilin's actual function.

## 5. A natural integer dilation preserves arithmetic but breaks exponential decay

Fix an integer d>=2 and deform the roots and poles together:

    R_(n,d)(t)=[product_(j=1)^n(t-dj)/product_(j=0)^n(t+dj)]²
             =d^-2 R_n(t/d),
    F_(n,d)=-sum_(t=1)^infinity R_(n,d)'(t).

The residues of the unsquared rational function remain the same integer binomial coefficients. If R_n has coefficients a_k,b_k, then

    R_(n,d)(t)=sum_k[a_k/(t+dk)²+(b_k/d)/(t+dk)].

Thus F_(n,d) still isolates ζ(3), with the same integral coefficient2A_n. The constant coefficient has denominator dividing D_(dn)³: dD_n divides D_(dn), and every shifted harmonic sum ends at dk<=dn. The reciprocal symmetry and right-half-plane bound also survive.

Nevertheless the exponential approximation mechanism fails.

**Theorem.** For each fixed integer d>=2,

    F_(n,d) ~ (4/d³)[Γ(1/d)/Γ(1-1/d)]²
                         n^(-4/d)log n.              (4)

In particular this sequence is eventually positive but decreases only polynomially. Thus preserving integer residues and algebraic zeta cancellation is not sufficient: the summation lattice must also align with the zeros and the telescoping boundary.

### 5.1 An exact forced recurrence identifies the defect

Let T be the usual Apéry operator and S_n=s_nR_n its short certificate. Scaling gives

    T R_(n,d)(t)=S_(n,d)(t+d)-S_(n,d)(t),
    S_(n,d)(t)=d^-2S_n(t/d).

Differentiate and sum over positive integers. The shift is now d, so d lower boundary terms remain:

    (T F_(_,d))_n=d^-3 sum_(r=1)^d S_n'(r/d)
                =d^-3 sum_(r=1)^(d-1) S_n'(r/d).     (5)

The last term at r=d vanishes, but the others lie strictly between0 and1, where R_n has no zero. This formula is exact, not an asymptotic prediction.

### 5.2 The first unmatched boundary dominates

For each fixed x∈(0,1), gamma products give

    R_n(x)=[Γ(x)/Γ(1-x)]²
                [Γ(n+1-x)/Γ(n+1+x)]²
           ~C(x)n^(-4x),
    C(x)=[Γ(x)/Γ(1-x)]²>0.

Its logarithmic derivative satisfies

    R_n'(x)/R_n(x)=-4log n+O_x(1).

The latter estimate also follows directly from the finite product: subtract2/j from1/(j-x)+1/(j+x), leaving a convergent series. The existence of a positive limit for n^(4x)R_n(x) likewise follows from the logarithmic product and H_n-logn→γ; gamma notation supplies its familiar closed value.

Since s_n(x)~32n³ and s_n'(x)=O_x(n),

    S_n'(x)~-128C(x)n^(3-4x)log n.

Among the finitely many x=r/d in(5), the smallest x=1/d dominates. Therefore

    (T F_(_,d))_n~-(128/d³)C(1/d)n^(3-4/d)log n.      (6)

Already this proves that F_(n,d) cannot be exponentially small: a polynomial-coefficient recurrence applied to an exponentially small sequence would remain exponentially small, contradicting(6).

### 5.3 A comparison argument proves the full asymptotic

For integer t>0 not equal to a zero, |R_(n,d)(t)|<=t^-2 and its logarithmic derivative has absolute value O(n+1), because all nonzero distances to integer roots and poles have absolute value at least one. At the roots R' is zero. Summing therefore gives F_(n,d)=O(n+1).

Set

    b_n=(4/d³)C(1/d)n^(-4/d)log n.

For any fixed power times logn, neighboring ratios tend to1, so

    (T b)_n~-32n³b_n,

which matches(6). Hence the residual g_n-(Tb)_n is o(n³b_n), where g_n is the right side of(5).

Choose any θ strictly between17-12sqrt2 and1. For sufficiently large n, T(θ^n)<0 and T(2^n)<0. For arbitrary ε>0, the positive barrier

    εb_n+Mθ^n+η2^n

has T-value at most minus the absolute residual after n is large enough; M can handle the fixed left boundary. The polynomial bound on F ensures that η2^n dominates at infinity. The same discrete maximum principle used in `truncated-derivative.md` bounds both signs of F-b by this barrier. Let η decrease to zero, divide by b_n and let n→infinity, then let ε decrease to zero. Since θ^n/b_n→0, the ratio F_(n,d)/b_n tends to1, proving(4).

The d=1 construction is exceptional in this deformation family: there are no unmatched lower boundary terms at all, so its recurrence is homogeneous and exponential decay becomes possible.

## 6. A concise design principle, with its limits

Within this ansatz, the main ingredients fit together as follows:

- Minimal double-pole budget plus initial double zeros forces the square.
- Consecutive integer poles minimize the universal difference denominator and the pole span, while factorial ratios make residues integral.
- Quadratic decay cancels the one unwanted zeta value after differentiation.
- Unit spacing aligns the discrete summation boundary with the zeros. Integer dilation preserves the first three properties but provably loses exponential decay.
- The actual reciprocal reflection explains a half-plane modulus property; additive reflection is a different and incompatible constraint when all original zeros are retained.

This is a structural explanation of why the construction is unusually well matched to ζ(3). It does not prove that no alternative elementary construction could handle ζ(5), nor that the chosen ansatz is optimal among all rational approximations.

## Reproducible checks

`codex/computations/check_zudilin_design_classification.py` verifies the nonconsecutive pole example, reciprocal reflection, scaled partial fractions and forced recurrence at exact finite instances, and high-precision asymptotic ratios for d=2,3,4. The classification, arithmetic minimum, and deformation asymptotic have general proofs above; finite checks are supplementary.
