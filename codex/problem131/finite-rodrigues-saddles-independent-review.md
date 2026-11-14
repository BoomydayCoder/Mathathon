# Independent exact audit of the finite higher-Rodrigues saddle certificates

**Verdict.** The finite certificates and their use in Section 4 of `eulerian-rodrigues-period-decomposition.md` pass independent review. For each \(m=2,3,4,5\), they prove the existence of every required strip saddle, the correct logarithm branches, and a unique dominant nonreal conjugate pair with modulus strictly below one. Combined with the independently reviewed conditional component theorem, this proves the displayed full-form oscillatory asymptotic, including the ζ(7), ζ(9), and ζ(11) cases. No mathematical correction was needed.

The independent verifier `check_finite_saddles_independent.py` uses only integer and rational arithmetic. It reads the saved rational centers but does **not** call a numerical root finder, use floating-point arithmetic, or trust the saved Rouché margins, modulus intervals, or branch labels without checking them.

## 1. Rouché inequalities genuinely isolate the roots

For either polynomial
\[
 P_\varepsilon(x)=(m-x)(1+x)^m-\varepsilon x^{m+1},
 \qquad\varepsilon\in\{1,-1\},
\]
the certificate compares \(P(c+z)\) with \(P'(c)z\) on \(|z|=r\). Its rational inequality is sufficient because
\[
 |P(c)|+\sum_{k\ge2}\frac{|P^{(k)}(c)|r^k}{k!}
 \le |P(c)|_1+\sum_{k\ge2}\frac{|P^{(k)}(c)|_1r^k}{k!},
\]
whereas \(|P'(c)|\ge\max(|\Re P'(c)|,|\Im P'(c)|)\). Thus the remainder is strictly smaller than the linear term everywhere on the circle. Rouché gives exactly one zero, counted with multiplicity, inside. In particular that zero is simple.

The generator checks disjoint disks in number equal to each polynomial's degree, so its approximate-root procedure is merely a way to propose rational centers. A wrong or duplicated proposal would fail an exact check. The use of a floating-point estimate to propose the branch integer in that script is similarly harmless: the strict rational branch inequality is what accepts it.

My independent verifier instead checks the selected component centers directly, using complex numbers represented as pairs of Fractions and Taylor coefficients \(\sum_{i\ge k}p_i\binom ikc^{i-k}\). It checks the root in every upper-half or real component disk, as well as disk disjointness. There are \(m\) such disks for each multiplicity, accounting by conjugation for \(2m-1\) component roots. Their counts in the two polynomial branches leave exactly one plus-branch root. The already proved positive-axis monotonicity identifies that remaining root as the unique positive real Laurent critical point. No component is omitted.

A real-centered disk has its unique root fixed by conjugation, hence real. Every nonreal selected disk has its imaginary center strictly larger than its radius. The disk bounds also prove \(0<\Re t<1\) after \(t=-x\); this is a prerequisite for the contour theorem, not an inference from the root's decimal approximation.

## 2. The logarithm branches are certified with rational intervals

For \(t=a+ib\) in an upper-half disk, the three required arguments are
\[
 A=\arctan(b/a),\quad B=\arctan(b/(1-a)),\quad
 D=\arctan(b/(a+m)).
\]
Thus
\(\Im\psi_m'=(2m+2)A+2mB-2D\).
All denominators are positive on the certified disks, so monotonicity of arctangent validates the endpoint ratios used in the interval calculation.

The alternating arctangent series has the stated remainder direction: an even number of terms ends with a negative term, giving a lower bound, and the next positive term bounds the error above. Reciprocal and quarter-turn reductions keep the small argument within \([0,1/2]\). Machin's formula gives a rigorous π interval from the same series. For clarity, the underlying identity follows from
\[
 \tan(4\arctan(1/5))=120/119,\qquad
 \frac{120/119-1/239}{1+(120/119)/239}=1,
\]
and the principal-angle range. My verifier uses 60 terms instead of the generator's 50 and independently recomputes all bounds.

The algebraic equation already implies \(\Re\psi_m'=0\) and \(\Im\psi_m'\in2\pi\mathbb Z\). A strict enclosure between \((2k-1)\pi\) and \((2k+1)\pi\) therefore identifies the value exactly as \(2k\pi\). The rational comparisons against lower and upper π bounds have the safe directions in the source script. The independently checked branch labels run through \(j=1-m,\ldots,0\); conjugation gives the positive labels. The sign relation \(\varepsilon=(-1)^{j+m+1}\) holds for every certified disk.

For the real disk, positivity of the ratio in the real saddle equation and its branch sign imply the full real-logarithm equation with \(j=0\). No argument interval is needed there.

## 3. Dominance and nonreal critical values

At a certified upper root,
\[
 |\eta|=\left[\frac{m^2((1-a)^2+b^2)}{(a+m)^2+b^2}\right]^m.
\]
The certificate uses the lower numerator with the upper denominator, and conversely, so its rational intervals are conservative despite correlations between \(a\) and \(b\). I recomputed these intervals from the disk coordinates. Each upper bound is below one, and the lower bound for \(j=1-m\) exceeds every other component's upper bound. Thus the full decomposition has exactly one dominant conjugate pair.

Nonreality is certified separately, rather than inferred from a nonreal saddle. Since
\[
 -\arg\eta/\pi\equiv2m(B+D)/\pi\pmod2,
\]
a strict interval between successive integers proves that \(\eta\) is not real. These inequalities passed for every noncentral component, including the dominant one.

The four decimal enclosures printed in the main note are also safe: the independent verifier parses each decimal endpoint as an exact rational number and checks that it strictly contains the recomputed rational modulus interval. No decimal rounding is used to prove dominance.

## 4. Passage to the full pure form

The earlier independent review established the conditional component formula
\[
 T_{n,m,j}=\mathcal C_m(t_j)n^{-m-1/2}\eta_j^n(1+O(n^{-1}))
\]
for each full saddle in the strip, and the exact Eulerian decomposition. The new certificates establish every missing finite hypothesis for \(m=2,3,4,5\). There are no poles to cross when moving the vertical line to \(\Re z=n\Re t_j\), since the entire move stays in \(0<\Re z<n+1\).

The outer Eulerian coefficients are both 1. The two corresponding components are complex conjugates, so their combined leading coefficient is exactly
\(2/(2m-1)!\), not \(1/(2m-1)!\). Every other component has a strict exponential modulus gap. With fixed \(m\), those finitely many terms can be absorbed into the relative \(O(1/n)\) error on the dominant scale. Hence
\[
 \mathscr L_{n,m}=
 \frac{2n^{-m-1/2}}{(2m-1)!}
 \Re\{\mathcal C_m(t_\star)\eta_\star^n\}
 +O(n^{-m-3/2}|\eta_\star|^n).
\]
The constant \(\mathcal C_m(t_\star)\) is nonzero: its prefactor has no zero in the strip, and \(\Re\psi_m''>0\). The coefficients of the real-part term are respectively \(1/3,1/60,1/2520,1/181440\) for \(m=2,3,4,5\), agreeing with the prior ζ(5) normalization.

Because \(\eta_\star\) is nonreal, its normalized phase rotation either has a dense orbit, or a finite orbit of at least three points. In the finite case the real parts of \(\mathcal C\) times that orbit sum to zero and cannot all vanish. Both cases supply infinitely many positive and negative leading terms bounded away from zero after division by \(|\mathcal C\eta_\star^n|\). The error is smaller on those subsequences. This proves both sign changes infinitely often and
\[
 \limsup_{n\to\infty}|\mathscr L_{n,m}|^{1/n}=|\eta_\star|<1.
\]
No assertion that every individual term is nonzero is needed or follows from this argument.

## 5. Scope, reproducibility, and the arithmetic limitation

Run `python codex/computations/check_finite_saddles_independent.py`. It verifies 14 real/upper component disks across \(m=2,3,4,5\), their conjugate accounting, all logarithm branch conditions, nonreal values, strict dominant gaps, and the safe printed modulus intervals. It passed with no floating-point operations. The saved centers come from `rodrigues_finite_saddles_certificate.json`; their original numerical provenance is irrelevant once the exact inequalities have been verified.

The new statements are analytic results for these four fixed multiplicities. They do not prove a corresponding result at \(m=6\) or for all odd zeta values. They also do not remove the denominator cost: the established sufficient normalizer \(D_{mn}^{2m+1}\) remains much larger than these analytic decay rates. The note correctly treats that as a gap in available bounds, not a lower-bound theorem for primitive denominators or an impossibility theorem for irrationality proofs.
