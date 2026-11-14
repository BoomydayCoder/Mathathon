# Adversarial audit of the extra inner-prime factor for zeta(7)

**Verdict.** No gap found in the claimed full-matrix improvement. The extra factor is coefficientwise in X and survives both paired residue sources and all cross-class entries. The parity argument yields **one whole p-power per selected class**, not merely one half-power. The saving1/25200 is a valid improvement of a sufficient arithmetic envelope; it does not establish the actual normalizer's asymptotic rate or prove zeta(7) irrational.

I read `inner-leading-block-cancellation.md`, the parent's `inner-leading-block-independent-review.md`, the local arithmetic reconstruction, and the fixed-odd-index distribution/weight extension. The verification below reconstructs the relevant entire-matrix argument rather than inferring it from a singular numerical self block.

## 1. Fix the exact finite-K allocation

Take s=7,r=6,K=40n,N=3n,h=37n,alpha=3/40,H=23/20, fixed M>=40, L_0=4M+10, and K>=200M^2. Let p be prime with

    3<K/p<70/23.

Then this is within the original inner range, p>200M, p^2>5K, floor(K/p)=3, and floor(N/p)=0. Put m=(p-1)/2 and S=HK-L_0. The CRT allocation is mT+E=S with0<=E<m.

There is a direct uniform proof of T=6. Namely

    S-6m=p(HK/p-3)-L_0+3>86M-7>0,
    S-7m=p(HK/p-7/2)-L_0+7/2<0.

Thus E=S-6m and the number of nonextra classes is exactly

    d=m-E=7(p-1)/2-HK+L_0.                         (1)

There are precisely K-3p classes of size ell_K=7, with representatives a<=K-3p; all other ordinary classes have size6. Moreover K-3p<N, so all size-seven classes lie in a<=N. Finally

    E-N=(H-alpha)K-3p-L_0+3>N-L_0+3>0.

Therefore assign extras to all a<=N, including all size-seven classes, and assign any remaining extras within a>N. Every nonextra class then has a>N and ell_K=6. This respects the required decreasing-ell_K allocation.

Moving extras among tied size-six classes does not change the old bound: their additional row has weight0 at s=7. The minimum ordinary cap in the zero-block weights also stays0 because d>0. Hence the comparison really is with the same old row-weight sum, not a secretly weakened reallocation.

## 2. Record the full matrix's weight sectors

For an ordinary class let delta=1_(a<=N), epsilon in{0,1} mark an extra, and Z=6+epsilon. Its dimension and weights are

    L_a=6-3delta+epsilon,
    w_(a,i)=i+3delta-(ell_K+6)/2, 0<=i<L_a.

For ell_K=6 all these weights are integers and the largest is epsilon-1, hence at most0. For ell_K=7 the class necessarily has an extra, its weights are strict half-integers, and its largest weight is-1/2.

The zero-block weights simplify to

    w_(0,i)=min(2i-5/2,0).

Thus its first two weights are-5/2,-1/2 and every later weight is0. In particular **every external integer-sector weight is nonpositive**. This fact makes the external-column analysis completely explicit.

The old coefficientwise bounds v_p^G G_uv>=w_u+w_v hold for the actual symmetric moment matrix in the CRT basis. The same s=7 distribution proof applies: p>7, its factor is p^-6, Y=p^7X+C with C in p^7 Z_p, and the zero-source degree is at most7+4L_0<p+1. Ordinary leading degrees are at most14. Thus the baseline local integrality hypotheses are valid under the displayed finite bounds, not merely in an unspecified eventual range.

## 3. The CRT leading units really can be removed

Fix a selected nonextra class a>N. Its dimension is6 and its weights are w_i=i-6. Write its rows as

    E_(a,i)(t)=F_a(t)(t+a^2)^i,
    F_a(t)=product_(c!=a)(t+c^2)^(L_c).

At v=a+pu, t=-v^2,

    t+a^2=-2ap u-p^2u^2,
    F_a(-v^2)=F_a(-a^2)+p*(polynomial in u)+...

and F_a(-a^2) is a p-unit, independent of i. Thus after extracting p^i the row's leading term is

    F_a(-a^2)(-2a)^i u^i.

The only index-dependent leading unit is(-2a)^i. It is independent of X. Removing these units after the row-weight scaling is legitimate. It is not necessary, and would be inaccurate, to claim that this weighted operation itself is unimodular in the unscaled original CRT basis.

The near-pole interval is exactly[-3,2], while I_N is empty: a>N>K-3p and a<p/2 give the necessary floor inequalities. Every near factor contributes the same leading coefficient-2a, so the common denominator is

    U(u)=(u+3)(u+2)(u+1)u(u-1)(u-2).

All other leading factors are units constant in u and common to the relevant selected-row indices.

The paired source p-a+pu is not a source of destructive or differently weighted interference. The full pulled-back entry v^7 R(-v^2) is odd, and

    p-a+pu=-(a+p(-1-u)),
    tau_(7,Y)(g(-1-u))=-tau_(7,Y)(g(u)).

The two minus signs cancel. The paired contributions are equal and add, introducing the p-unit2. This applies to the complete pulled-back bilinear entry, before one truncates its local expansion.

## 4. Which columns can contribute at the selected source?

There are exactly three cases after the old weights are factored.

1. A column in the selected class has leading monomial u^j for0<=j<=5. These are the six self columns.
2. A column in another CRT block vanishes at t=-a^2 to order exactly L_a=6. Its local half-weight at a is0. If its assigned integer weight is negative, the normalized contribution is already divisible by p, so there is nothing to cancel.
3. An external integer column of assigned weight0 can contribute a leading term. Its monomial is a unit multiple of u^6, with the multiplier independent of the selected row index. This includes the zero-block columns of weight0.

Strict-half-integer external columns are controlled by parity in section7 and do not need to be annihilated modulo p as integer-sector columns.

Thus checking columns u^0,...,u^6 is sufficient for the **whole integer-sector block** at this source. No leading u^7 border column is present. This degree limit is essential: the annihilator below would not kill the next moment, since tau_7(u^7/[u(u+1)])=1/6.

## 5. The local expansion reduction is coefficientwise valid

For the selected row indices used below, i<=4. Combining them with a potentially relevant self/border column gives numerator degree at most i+6<=10. Each first-order correction in p increases that degree by at most one:

- (a+pu)^7 has first correction proportional to u.
- A far factor or its reciprocal has a unit constant term and first correction affine in u.
- A near factor splits into its exact linear pole factor and a unit factor with first correction affine in u.
- Powers of CRT factors have the same property, regardless of how large their exponents are.

Consequently the first two numerator coefficients U_0,U_1 in the local expansion have degrees at most10 and11 respectively. Dividing U_1 by the degree-six U leaves a polynomial quotient of degree at most5. Under tau_7, degrees below5 vanish and the degree-five moment is1/6, which is p-integral. The residues have unit denominator differences and poles only at-3,...,2, whose harmonic values are p-integral. Thus tau_(7,Y)(U_1/U) is integral.

This checks the otherwise dangerous pU_1 term. The full functional can lose one p-power on arbitrary polynomial degrees; integrality of the whole expansion would not alone justify dropping pU_1 modulo p. Here that additional check is explicit. Terms p^jU_j with j>=2 vanish modulo p even under the general one-power loss.

Since Y=p^7X+C has all coefficients divisible by p, the leading reduction is tau_(7,0), with **constant coefficients independent of X**. Every assertion is a Gauss-valuation assertion about all polynomial coefficients, not an identity obtained only after evaluating X.

## 6. The vector kills all other sources as well

Set

    q(u)=U(u)/(u(u+1))=u^4+2u^3-7u^2-8u+12.

For0<=j<=6, polynomial division gives

    tau_(7,0)(q(u)u^j/U(u))
      =tau_(7,0)(u^j/[u(u+1)])=0.

The quotient has degree at most4, below the fifth derivative, and the remaining simple pole values at0 and-1 are both zero. After the diagonal-unit correction from section3, its coefficient vector therefore kills every potentially leading integer-sector column at the paired sources a,p-a.

Now take another ordinary source c. A selected original row of index i<=4 vanishes there to order L_c. Its local half-weight exceeds its assigned weight by

    Z_c-ell_K(c)/2+3-i.

For a size-six class Z_c is6 or7, so this is at least2. For a size-seven class Z_c=7, it is at least5/2. The other column's local half-weight is at least its assigned weight by the baseline proof. Therefore every such source contributes a p-divisible normalized entry against an integer-sector column, before any cancellation is used.

At the zero source, the selected row has surplus

    (2L_0-5/2)-(i-6)=2L_0+7/2-i>=2L_0-1/2.

It likewise cannot restore a nonzero leading entry. This completes the cross-class check: the vector is in the left kernel of the reduction of the **entire** integer-sector principal block, not only of its six-by-six self block.

The coefficient of u^4 is1. After dividing by its diagonal unit it still has a unit pivot, so the vector lifts to a unimodular row operation over Z_p. The supports for different selected classes are disjoint, yielding d independent constant kernel vectors and a simultaneous unimodular completion.

## 7. The parity lift gives exactly the asserted determinant lower bound

Let pi=sqrt(p), and scale the CRT matrix symmetrically by p^-w_i, obtaining B. All entries are integral over the ramified ring. Integer-sector and half-sector principal entries lie in Z_p[X], while entries between the two sectors lie in pi Z_p[X].

Apply the lifted kernel operations, and their transposes, **inside the integer sector**. Each selected row then has entries in p Z_p[X] against every integer-sector column, and in pi Z_p[X] against every other column. Symmetry gives the corresponding column statement. Selected-selected entries are in p Z_p[X].

Divide each selected row and column by pi. All entries remain integral: selected-selected entries are divided by p, selected-other-integer entries by pi, and selected-half entries by pi. Hence the determinant before this final scaling has valuation at least d. Restoring the old diagonal factors gives

    v_p^G(det G)>=2sum_i w_i+d=gamma_old+d.            (2)

This step uses symmetry and the integer parity of the selected supports. A generic null vector over the ramified residue field would not alone justify the same conclusion. The reviewed argument meets both requirements.

## 8. Budget conversion and checks

From(1), with fixed M,

    d=p(7/2-(23/20)(K/p))+O_M(1).

The usual fixed-interval prime-number-theorem conversion therefore lowers the sufficient arithmetic-cost envelope by

    integral_3^(70/23) (7/2-23x/20)/x^3 dx
       =(6*(23/20)-7)^2/252=1/25200.

The finite reserve term contributes o(K^2); interval endpoints affect only lower-order terms. This is a genuine improvement of a proved local lower bound for the determinant. Neither the old nor new envelope is asserted to equal the actual minimal arithmetic cost.

I independently reran `check_inner_kernel_allocation.py`: all755 finite-K prime patterns passed, including exact counts and the rational integral. I also reran `check_inner_leading_blocks.py`: the six-by-seven kernel and all listed exact ranks passed. The full-matrix conclusion rests on sections3–7 above, not on those finite rank checks.

The revised scalar budget remains far short of a zeta(7) irrationality proof. The valid discovery is the mechanism: low-degree polynomial moments and the two vanishing pole values cause an actual coefficientwise determinant cancellation that the generic row-weight bound does not capture.
