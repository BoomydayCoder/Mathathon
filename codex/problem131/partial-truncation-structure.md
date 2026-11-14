# Partial truncation: exact shape, reciprocal saddles, and two crossings

Status: exact cutoff identities and nonvanishing of the harmonic bracket proved here. The contour identity was independently checked against the integral agent's derivation. The crossing asymptotics below combine that identity with the already proved Apéry asymptotic. No historical novelty claim; no irrationality theorem for zeta(5) is asserted.

Use R_n(t)=[(t-1)...(t-n)/(t(t+1)...(t+n))]². Let

    E_(n,k)=-(1/24)sum_{t=k+1}^infinity R_n'''(t), 0<=k<=n,
    a_(n,j)=[(j-1)!²(n-j)!/(n+j)!]²,
    h_(n,j)=H_(n-j)+H_(n+j)-2H_(j-1),
    c_(n,j)=a_(n,j)h_(n,j)/2=-R_n'''(j)/24.

Then

    E_(n,k)-E_(n,k-1)=-c_(n,k).                         (1)

This note distinguishes the **minimum value** of E as a function of k from the **smallest absolute value**. The first lies near the saddle; the second lies near one of two zero crossings a growing distance away.

## 1. Exact shape and unique minimum

For 1<=j<n,

    h_(n,j+1)-h_(n,j)
      =-1/(n-j)+1/(n+j+1)-2/j<0.                      (2)

Thus the bracket strictly decreases. It starts positive. For n>=3 its last value H_(2n)-2H_(n-1) is negative: it is negative at n=3 and decreases with n. Hence E strictly decreases until the bracket changes sign and strictly increases thereafter, provided no bracket is zero.

**Theorem: h_(n,j) is never zero.** If j<=n/2, both H_(n-j) and H_(n+j) exceed H_(j-1), so the bracket is positive. If j>n/2, let 2^v be the largest power of two not exceeding n+j. Then

    2^v>(n+j)/2>=max(n-j,j-1).

The harmonic number H_(n+j) has exact 2-adic valuation -v: only the summand 1/2^v has that denominator exponent. Both remaining harmonic terms have strictly larger 2-adic valuation, so their addition cannot cancel it. In fact v_2(h_(n,j))=-v in this range. Therefore the bracket is nonzero. QED.

Consequently E_(n,k) has a unique discrete minimum. For n>=3 its index is

    k_min=max{j:h_(n,j)>0},

strictly between 0 and n. For n=1,2 all brackets are positive, and the minimum is at k=n. There can be at most two sign crossings. The location of the minimum depends only on rational harmonic numbers, not on the value of zeta(5).

This exact shape makes a full scan unnecessary: binary search on the strictly decreasing harmonic bracket locates the minimum; the two monotone branches then locate any sign crossings.

## 2. Gamma interpolation and an almost explicit cutoff

For real 0<t<n+1 define

    U_n(t)=Gamma(t)² Gamma(n+1-t)/Gamma(n+1+t),
    a_n(t)=U_n(t)²,
    h_n(t)=psi(n+1-t)+psi(n+1+t)-2psi(t).

Here psi is the logarithmic derivative of Gamma. At integer j these agree with a_(n,j),h_(n,j). Differentiation gives

    a_n'(t)=-2a_n(t)h_n(t),
    c_(n,j)=-a_n'(j)/4.                               (3)

Moreover

    h_n'(t)=-psi_1(n+1-t)+psi_1(n+1+t)-2psi_1(t)<0,

because the trigamma function psi_1 is positive and strictly decreasing. These facts follow directly from its convergent series sum_{r>=0}(x+r)^(-2). Thus log a_n is strictly convex and has one minimum, at a point t_* where h_n(t_*)=0.

The exact shift ratio is particularly simple:

    a_n(t+1)/a_n(t)=t^4/[(n-t)²(n+t+1)²].              (4)

Set

    t_0=[sqrt(1+8n(n+1))-1]/4.

Then a_n(t_0+1)=a_n(t_0). Strict log convexity places the continuous minimum strictly between t_0 and t_0+1. Since h_(n,j) never vanishes, for n>=3

    k_min=floor(t_*),
    k_min belongs to {floor(t_0),ceil(t_0)}.           (5)

When t_0 is an integer the two candidates coincide. One exact harmonic comparison therefore chooses the minimum from at most two explicit candidates.

The standard digamma expansion, or harmonic asymptotics at nearby integers, gives

    t_*=n/sqrt2+(1+sqrt2)/4+O(1/n).                   (6)

Indeed, for t=xn,

    h_n(t)=log[(1-x²)/x²]
             +[1/(1-x²)+1/x]/n+O(n^-2).

The leading zero is x=1/sqrt2, its derivative is -4sqrt2, and the next term there is 2+sqrt2. This yields (6).

## 3. Exact reciprocal link to the Apéry binomial saddle

If

    v_(n,j)=binom(n,j)²binom(n+j,j)²,

then factorial cancellation gives

    a_(n,j)=1/[j^4 v_(n,j)].                          (7)

Thus the tiny correction weights near their minimum are reciprocal to the large binomial summands in the Apéry sequence, up to the simple factor j^4. The same saddle n/sqrt2 appears for a concrete algebraic reason, rather than as a numerical coincidence.

There is also an exact weighted convexity identity in the cutoff direction:

    [E_(n,k+1)-E_(n,k)]/a_(n,k+1)
      -[E_(n,k)-E_(n,k-1)]/a_(n,k)
      =-[h_(n,k+1)-h_(n,k)]/2>0.                     (8)

This is a second-order difference equation in k with explicit rational forcing. It supplies shape control, although it alone cannot determine the vertical offset involving zeta(5).

## 4. Independent audit of the contour identity

Let c=k+1/2 and let F_n=-sum_{t>=1}R_n'(t), the positive Apéry linear form for zeta(3). The integral agent derived

    E_(n,k)=-zeta(2)F_n
       +(pi/8) integral_{-infinity}^infinity sech²(pi y)U_n(c+iy)² dy.   (9)

Here is an independent check of all signs and constants.

Reflection of Gamma gives

    R_n(z)=sin²(pi z)U_n(z)²/pi².

For K(z)=pi cot(pi z), close the upward vertical line Re z=c to the right. The contour is clockwise, so the sum of residues is the negative of the vertical integral. R_n has no poles in this half-plane. Its derivatives decay rationally, and on half-integer closing boundaries the cotangent kernel is bounded; horizontal and right boundaries vanish as they recede.

For the first derivative, all omitted terms at 1,...,k vanish because they are derivatives at double zeros. Hence

    F_n=(1/(2pi i)) integral K(z)R_n'(z)dz
        =-(1/(2pi i)) integral K'(z)R_n(z)dz
        =(1/(2pi)) integral U_n(c+iy)²dy,

since K'=-pi²csc²(pi z). Integration-by-parts boundary terms vanish because R_n=O(z^-2) and the kernel derivatives are bounded or exponentially small on the vertical line.

For the third derivative,

    E_(n,k)=-(1/(48pi i)) integral K'''(z)R_n(z)dz.

Using

    K'''=-2pi^4 csc²(pi z)[3cot²(pi z)+1]

and cot(pi(c+iy))=-i tanh(pi y), this becomes

    (pi/24) integral [1-3tanh²(pi y)]U_n(c+iy)²dy.

Finally 1-3tanh²=-2+3sech² and the first-derivative integral give (9), with -pi²F_n/6=-zeta(2)F_n. Thus the contour formula's constants, orientation, and omitted-zero terms all check.

The beta-factorization

    U_n(z)=B(z,z)B(n+1-z,2z)

is valid on this line because 0<c<n+1. Taking absolute values in the two beta integrals yields

    |U_n(c+iy)|<=U_n(c).

Consequently the correction in (9) has absolute value at most a_n(c)/4, since integral sech²(pi y)dy=2/pi. The integral is real by conjugation symmetry.

## 5. Why the minimum becomes negative and has the Apéry decay base

Put rho=(sqrt2-1)^4. The Apéry asymptotics already proved in `apery-asymptotics.md` give

    F_n ~ C_F rho^n n^(-3/2),
    C_F=2^(3/4)pi^(3/2)/(1+sqrt2)².

Uniform Stirling estimates near the minimum give

    a_n(c) ~ C_a rho^n n^-2,
    C_a=8pi²/(1+sqrt2)²,

whenever c=n/sqrt2+O(1). Thus at k_min, the correction integral in (9) is O(rho^n n^-2), smaller than F_n by a factor n^-1/2. Therefore

    E_(n,k_min)/F_n -> -zeta(2),
    E_(n,k_min) ~ -zeta(2)C_F rho^n n^(-3/2).           (10)

The minimum is strictly negative eventually, and its exponential rate is exactly rho. Full truncation at k=n has the slower decay base 1/16. This quantifies the benefit of stopping near the reciprocal saddle.

The conclusion concerns the minimum value, not the minimum absolute value. Near a zero crossing the latter can be smaller by cancellation.

## 6. Two crossings occur on the sqrt(n log n) scale

The contour correction admits a useful uniform refinement. If c_n/n->1/sqrt2, then for every fixed y,

    U_n(c_n+iy)/U_n(c_n)->1.

This follows from the digamma expansion: the first logarithmic derivative tends to zero and higher derivatives are O(1/n) locally. The beta bound gives domination by sech²(pi y), so dominated convergence in (9) yields

    E_(n,k)=-zeta(2)F_n+(1+o(1))a_n(c)/4              (11)

uniformly when |c-n/sqrt2|<=C sqrt(n log n), for fixed C. Uniformity follows by applying the same dominated-convergence argument to every sequence in the window.

In that window, Taylor expansion of log a_n at t_* and uniform Stirling estimates give

    a_n(c)=(1+o(1))C_a rho^n n^-2
                   exp[4sqrt2(c-t_*)²/n].            (12)

The cubic Taylor remainder is O((log n)^(3/2)/sqrt n), hence tends to zero. The shift of t_* from n/sqrt2 is O(1) and is already recorded in (6).

Balancing the positive term in (11) with zeta(2)F_n shows that the crossing points satisfy

    (c-t_*)²
      = n/(4sqrt2) [ (1/2)log n+log K+o(1) ],
    K=4zeta(2)C_F/C_a=pi^(3/2)/(6*2^(1/4)).           (13)

In particular their distance from the minimum is asymptotic to

    sqrt[n log n/(8sqrt2)].                          (14)

To make the sign conclusion rigorous, choose the squared distance in (13) with an added fixed epsilon>0 or -epsilon. The ratio of the positive correction to zeta(2)F_n tends to exp(epsilon) or exp(-epsilon). Thus outer points are positive and inner points negative. The exact monotonic shape gives one crossing on each branch. Adjacent integer cutoffs change the bracket in (13) by o(1), so the same asymptotic applies to the discrete crossing locations, allowing an exact zero at one endpoint.

Since c=k+1/2, a useful numerical prediction is

    k approximately n/sqrt2+(sqrt2-1)/4
       +/- sqrt{ n/(4sqrt2)[(1/2)log n+log(pi^(3/2)/(6*2^(1/4)))] }.

It predicts approximately 12.14 on the left for n=20, 25.03 on the left for n=40, 61.91 on the right for n=80, and 105.2 on the left for n=160. These correspond closely to the parent's observed smallest-absolute-error cuts. Formula (13), rather than those finite numerical matches, is the asymptotic theorem.

**A sign-controlled small remainder.** Choose the first strictly negative cutoff on the left branch. Its predecessor is nonnegative. Equation (1), together with (11)–(13), then bounds its absolute error by one step:

    0<-E_(n,k)<=O(rho^n sqrt(log n)/n²).               (15)

Indeed a_n(k)=O(F_n) at the crossing and h_n(k)=O(sqrt(log n/n)). This proves an upper bound and nonvanishing for a sign-selected sequence. It does not establish a matching lower bound for the globally smallest absolute error: unusually close integer crossings could lead to further cancellation.

The cutoff selection in (15) is defined through the sign of the actual linear form; no claim is made that a simple fixed rounding rule realizes it for every n. The explicit near-saddle choice in (10) already gives a nonzero sequence with the optimal proved base rho and does not require such sign selection.

## 7. Arithmetic remains separate

Partial truncation changes the rational coefficient while leaving A_n as the coefficient of zeta(5). The improved decay does not by itself control the new factorial denominators in the correction. The identities above therefore establish analytic structure, not irrationality.

The exact reciprocal identity (7) explains the tradeoff especially clearly: the small correction comes from reciprocals of large binomial coefficients. Those reciprocals help analysis while potentially harming arithmetic. The denominator investigation must account for cancellations in their sum, not simply ignore them.

## Verification

`check_partial_truncation_structure.py` checks the strict harmonic ordering, nonvanishing and 2-adic assertion, reciprocal identity, discrete cutoff candidates, and high-precision finite partial-tail shape. Gamma/contour and asymptotic arguments are proved above rather than inferred from numerical scans.
