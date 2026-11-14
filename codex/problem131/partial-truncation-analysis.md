# Partial truncation recovers the Apéry exponential rate

Status: independently derived construction and proof in this session. A separate odd-zeta agent audited Sections 1–7, including contour orientation, convergence, beta strip, Stirling estimates and gamma-phase/Fourier constants; no substantive error was found (see `odd-zeta-independent-review.md`, Section 8). Historical novelty unresolved. This note proves analytic improvement but does not claim that its rational coefficient has adequate denominators for an irrationality proof. Arithmetic is treated separately in `partial-truncation-arithmetic.md`.

## 1. Definitions and the main result

Let

    R_n(t)=[∏_{j=1}^n(t-j)/∏_{j=0}^n(t+j)]²,
    F_n=-Σ_{t=1}∞R_n'(t),
    E_(n,k)=-(1/24)Σ_{t=k+1}∞R_n'''(t),  0≤k≤n.

The F_n normalization is Zudilin's: F_n=2A_nζ(3)-2B_n, and F_n equals the positive Beukers triple integral. The partial higher-derivative tail has the form

    E_(n,k)=A_nζ(5)-B_(n,2)-C_(n,k),
    C_(n,k)=(1/2)Σ_{j=1}^k P_n'(j)²
                    [H_(n-j)+H_(n+j)-2H_(j-1)].

Here A_n=Σ binom(n,j)²binom(n+j,j)². Put ρ=17-12sqrt2=(sqrt2-1)^4.

**Theorem.** For any integer sequence k_n=n/sqrt2+O(1),

    E_(n,k_n)/F_n = -ζ(2)+O(n^-1/2).                  (1)

Consequently these partial tails are eventually negative and

    lim |E_(n,k_n)|^(1/n)=ρ.                          (2)

Thus cutting near the middle saddle genuinely improves the fully truncated rate 1/16 to the Apéry rate ρ≈0.029437. Positivity is not preserved: the optimal central tail has the opposite sign to the full tail. Nonvanishing nevertheless follows from (1).

The proof rests on an exact identity and inequality valid for every n,k.

## 2. A gamma-function envelope

For 0<Re z<n+1 define

    U_n(z)=Γ(z)² Γ(n+1-z)/Γ(n+1+z).

The reflection identity Γ(w)Γ(1-w)=π/sin(πw) gives

    R_n(z)=sin²(πz) U_n(z)²/π².                      (3)

To check the integer shift sign before squaring, write

    P_n(z)=Γ(z)²/[Γ(z-n)Γ(z+n+1)]
          =(-1)^n sin(πz)U_n(z)/π.

The sign disappears in R_n. At integer 1≤j≤n, U_n(j)=|P_n'(j)|. Thus U_n² interpolates the reciprocal-binomial envelope of the initial-zero corrections.

There is also a positive-integral factorization:

    U_n(z)=B(z,z) B(n+1-z,2z).                       (4)

Indeed the Γ(2z) factors cancel. For c∈(0,n+1), both beta integrals converge absolutely. Taking absolute values inside each integral gives the useful uniform inequality

    |U_n(c+iy)|≤U_n(c)  for all real y.               (5)

No asymptotic approximation is involved. In particular log U_n(c) is convex in c, either by Hölder's inequality applied to the positive beta-product integral or by regarding (4) as one positive Mellin-type integral.

## 3. Sum derivatives using a vertical contour

Fix k∈{0,...,n} and c=k+1/2. The rational function R_n has poles only at nonpositive integers. For r=1 or r=3, close the line Re z=c to the right, with the vertical line traversed upwards. This is clockwise orientation. The poles of πcot(πz) inside are exactly the integers k+1,k+2,..., each with residue one. Consequently

    Σ_{t=k+1}∞R_n^(r)(t)
       =-1/(2πi)∫_(c-i∞)^(c+i∞) πcot(πz)R_n^(r)(z)dz. (6)

For a direct justification, first take rectangles with right edge at a half-integer tending to infinity and horizontal edges at imaginary parts tending to ±∞. On these edges cot is uniformly bounded; R_n^(r)=O(|z|^(-r-2)) for fixed n. The other three edges tend to zero, and the vertical integral is absolutely convergent. This avoids an unproved contour-at-infinity assumption.

Since R_n'(j)=0 at j=1,...,n, the sum for r=1 is independent of k and equals -F_n. Let K(z)=πcot(πz). Integration by parts along the vertical line has no endpoint contributions, because derivatives of R_n decay and K and its fixed derivatives are bounded on that line. Thus

    F_n=1/(2πi)∫ K R_n' dz
       =-1/(2πi)∫ K'R_n dz
       =1/(2πi)∫ U_n(z)² dz
       =1/(2π)∫_R U_n(c+iy)²dy,                     (7)

using K'=-π²csc²(πz) and (3). The integral on the last line is real, although its integrand need not be pointwise real.

Similarly,

    E_(n,k)=1/(48πi)∫ K R_n''' dz
           =-1/(48πi)∫ K'''R_n dz.

The exact derivative identity is

    K'''=-2π^4 csc²(πz)[3cot²(πz)+1].

On Re z=k+1/2, cot(πz)=-i tanh(πy). Substituting this, (3), and dz=i dy gives

    E_(n,k)=π/24∫_R [1-3tanh²(πy)]U_n(c+iy)²dy.

Since 1-3tanh²=-2+3sech², equation (7) proves the central exact identity

    E_(n,k) = -ζ(2)F_n
         + (π/8)∫_R sech²(πy)U_n(c+iy)²dy.           (8)

The coefficient is ζ(2)=π²/6, not 1/6; the factor 1/(2π) in (7) is crucial. This normalization was independently checked against exact rational/high-precision evaluation of both linear forms.

Combining (5), (8), and ∫_R sech²(πy)dy=2/π gives

    |E_(n,k)+ζ(2)F_n|≤U_n(k+1/2)²/4.                (9)

This is an exact inequality for every permitted cut.

## 4. Bound the envelope at the saddle

For integer j between 1 and n,

    U_n(j)²=1/[j^4 binom(n,j)²binom(n+j,j)²].

Set a=j/n. Stirling's formula, uniformly for a in a compact subset of (0,1), gives

    U_n(j)²~ [4π²(1-a)/(a²(1+a))] n^-2 exp(-2n Φ(a)),
    Φ(a)=(1+a)log(1+a)-(1-a)log(1-a)-2a log a.         (10)

The derivative is

    Φ'(a)=log[(1-a²)/a²],

so the unique maximum is at a_0=1/sqrt2, with

    Φ(a_0)=2log(1+sqrt2),   exp(-2Φ(a_0))=ρ.

If j=n/sqrt2+O(1), then n(Φ(j/n)-Φ(a_0))=O(1/n), since Φ'(a_0)=0. Thus U_n(j)²=O(n^-2ρ^n). For c=k+1/2 with k=n/sqrt2+O(1), log convexity bounds U_n(c) between the geometric interpolation of U_n(k),U_n(k+1), in particular by their maximum. Therefore

    U_n(c)²=O(n^-2ρ^n).                              (11)

Only this upper bound is required for the main theorem.

## 5. A simple lower bound for the positive Beukers integral

The companion integral proof gives

    F_n=∫_[0,1]^3 h(x,y,z)^n/[1-z+xyz] dxdy dz,
    h=x(1-x)y(1-y)z(1-z)/(1-z+xyz).

Its unique maximum is ρ at the interior point

    x=y=sqrt2-1,  z=1/sqrt2.

The function is smooth near that point, and its gradient vanishes there. Taylor's theorem therefore gives h≥ρ-C||v||² in a fixed small neighborhood, where v is displacement from the maximizing point. Integrate over a ball of radius d/sqrt n, with fixed sufficiently small d. On that ball h^n≥cρ^n for a constant c>0, while its volume is a constant times n^-3/2. The positive denominator is bounded above and below there. Hence

    F_n≥c_1 n^-3/2ρ^n.                              (12)

The already-proved upper bound F_n≤2ζ(3)ρ^n, together with (12), shows F_n^(1/n)→ρ. Combining (9), (11), and (12) gives (1), including its O(n^-1/2) error, and then (2).

This avoids a detailed three-dimensional Laplace expansion; a small ball near the maximum is enough.

## 6. More precise local correction and general fixed cut proportions

Suppose k_n/n→a∈(0,1), put c_n=k_n+1/2, and set

    b=log[(1-a²)/a²].

For every fixed real y, the standard gamma ratio asymptotic Γ(x+iy)/Γ(x)~x^(iy) as x→+∞ gives

    U_n(c_n+iy)/U_n(c_n)→exp(-iby).

The uniform bound (5) permits dominated convergence in the sech² integral. Its Fourier transform is

    ∫_R sech²(πy)exp(-2iby)dy=2b/[π sinh b],

with limiting value 2/π at b=0. For example, set u=exp(2πy) and evaluate the resulting beta integral; the identity Γ(1+iv)Γ(1-iv)=πv/sinh(πv) completes the calculation.

Thus (8) yields the sharper asymptotic

    E_(n,k_n)+ζ(2)F_n
       ~ U_n(c_n)²· b/[4sinh b],                    (13)

where b/sinh b is interpreted as 1 at b=0. The coefficient is always positive.

For k_n=an+O(1) with fixed a≠1/sqrt2, the envelope has exponential rate exp(-2Φ(a)), strictly larger than ρ, so F_n is exponentially smaller than the correction (13). Therefore

    E_(n,k_n)>0 eventually,
    lim E_(n,k_n)^(1/n)=exp(-2Φ(a)).                  (14)

In contrast, for a bounded-distance central cut, (1) shows eventual negativity. There must therefore be sign transitions between the central negative region and either fixed-proportion side. Locating them requires comparison of the envelope size with F_n, not just the sign of individual root terms. An independent companion study addresses the discrete minimizer and the transition scale.

At central cuts, (13) specializes to

    E_(n,k_n)+ζ(2)F_n~U_n(c_n)²/4.

Numerical ratios to this correction are approximately .9560,.9770,.9881,.9939 for n=10,20,40,80 with k=floor(n/sqrt2).

## 7. What has and has not improved

The main analytic improvement is rigorous: partial truncation near the saddle achieves the same exponential base ρ as the original ζ(3) approximation, whereas complete truncation has base 1/16. The exact identity (8) explains why a ζ(2) multiple of the ζ(3) linear form governs the central ζ(5) approximation error. This is a connection between higher derivatives, contour kernels, beta integrals and the original Beukers form, not just a matching numerical rate.

However, the rational correction still contains factorial denominators; its arithmetic must be assessed independently. A very small raw error does not imply a small integer after clearing the reduced denominator. The separate arithmetic study includes exact examples demonstrating this distinction. Neither (1) nor (8) proves ζ(5) irrational.

The contour proof invokes standard complex analysis, gamma reflection, beta integrals and Stirling/gamma-ratio asymptotics. It is an advanced structural explanation. It does not replace the elementary Beukers proof of ζ(3) with fewer prerequisites.

## Reproducibility and review

`codex/computations/check_partial_truncation.py` evaluates the exact rational coefficients before high-precision cancellation. It scans fixed cut proportions, checks the predicted central ratio and the envelope inequality, and tests the contour formula numerically for small n. These are supporting checks; the residue calculation, beta bound and local maximum argument provide the proof.

Independent review has checked contour orientation, the π² normalization, convergence of vertical integrals, and the exact strip needed for both beta factors. Section 8 below was subsequently added from the parent’s proposal and independently checked algebraically by the integral agent.

## 8. Uniform structural extension to all fixed odd derivative orders

The parent agent proposed the following extension; the integral agent independently checked its coefficient recurrence and contour normalization. Let m≥1 be fixed, q=2m-1, and define

    E_(n,k,m)=-(1/(2m)!)Σ_{t=k+1}∞R_n^(2m-1)(t).

Define integer coefficients c_(m,j), 1≤j≤m, by

    (d/dθ)^(2m-1)cot θ=Σ_{j=1}^m c_(m,j)csc^(2j)θ.

They satisfy

    c_(1,1)=-1,
    c_(m+1,j)=-4j²c_(m,j)+2(j-1)(2j-1)c_(m,j-1),

with missing indices interpreted as zero. To prove this recurrence, differentiate twice and use

    (d/dθ)² csc^(2j)θ
       =2j(2j+1)csc^(2j+2)θ-4j²csc^(2j)θ.

In particular,

    c_(m,1)=(-1)^m4^(m-1).

Repeated integration by parts in (6), using q odd, gives the exact identity

    E_(n,k,m)=a_m F_n
       -π^(2m-2)/[2π(2m)!]
          Σ_{j=2}^m c_(m,j)∫_R sech^(2j-2)(πy)U_n(c+iy)²dy,
    a_m=(-1)^(m+1)(2π)^(2m-2)/(2m)!.                 (15)

For m=1 the sum is empty and E_(n,k,1)=F_n/2, as required by its normalization. For m=2 it is exactly (8). For every fixed m≥2, (5) bounds the remaining integrals by constants times U_n(c)², since all powers of sech in the sum are integrable. Hence

    |E_(n,k,m)-a_m F_n|≤C_m U_n(c)².                 (16)

At any bounded-distance saddle cut k=n/sqrt2+O(1), (11)–(12) imply

    E_(n,k,m)/F_n=a_m+O_m(n^-1/2).                   (17)

Thus every fixed odd-zeta partial-tail family has the same optimal-cut exponential rate ρ, and its eventual sign is (-1)^(m+1). The first constants are

    a_1=1/2,
    a_2=-π²/6=-ζ(2),
    a_3=π^4/45=2ζ(4),
    a_4=-π^6/630=-(3/2)ζ(6).

This is a uniform relation between higher odd-zeta approximation errors and the ζ(3) linear form, with even-zeta multipliers. It still gives no denominator improvement by itself. Each new higher derivative requires its own arithmetic analysis or the general denominator bound from `higher-truncated-tail.md` adapted to the partial cut.
