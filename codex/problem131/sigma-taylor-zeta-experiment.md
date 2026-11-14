# Parameter Taylor coefficients isolate higher odd zeta values

Status: exact finite experiment and proved local denominator lemma. No all-index convergence or sign theorem for the Taylor coefficients is claimed. The positive deformed measure and its rational approximants are established inputs from [the deformation note](logarithm-zeta3-wilson-deformation.md). The finite experiment is accepted in [the independent review](sigma-taylor-independent-review.md), which reconstructed all twenty jets and rational certificates separately. No historical novelty claim is made; the parent investigation identified the generic diagonal recurrence in Askey (1984), equation 22, so the deformation should not be described as historically new.

Let \(q_n(\sigma)=b_n(\sigma)/A_n(\sigma)\) use exactly the associated numerator defined in the deformation note. The target has the convergent parameter expansion
\[
 Z_\sigma=\sum_{h\ge1}\frac1{h(h^2-\sigma)}
 =\sum_{r\ge0}\zeta(2r+3)\sigma^r,\qquad |\sigma|<1.
\]
We compute \(q_{n,r}=[\sigma^r]q_n(\sigma)\) for \(r=1,2\), \(1\le n\le10\). These are exact rational jets, rather than numerical derivatives. The jet construction truncates each factor \((h^2-\sigma)^{-1}\) after degree 2, constructs the full polynomial \(\Phi_n(t)\), checks its evaluation at \(-(n+1/2)^2\), and checks every coefficient of the polynomial difference quotient defining \(b_n\). Moments are the explicit Bernoulli-polynomial parameter jets in [the generic recurrence note](generic-sigma-diagonal-recurrence.md).

## Finite results

At index 1 the whole rational function is
\[
 q_1(\sigma)=\frac1{1-\sigma}+\frac1{5-\sigma}.
\]
Consequently \(q_{1,1}=26/25>\zeta(5)\), whereas \(q_{1,2}=126/125<\zeta(7)\). Positivity of the original parameter-dependent error therefore cannot justify a common one-sided sign for its Taylor coefficients.

Every error sign and primitive-error assertion below is certified by an exact Euler–Maclaurin rational interval for the target: \(N=80,L=60\), with the Bernoulli-polynomial coefficient-sum remainder bound documented in the base-Gaussian experiment. Decimal values are summaries of those certified rational intervals. Put \(B=\operatorname{den}(q_{n,r})\); the primitive error is \(B|q_{n,r}-\zeta(2r+3)|\).

| n | target | denominator digits | side | absolute error | primitive error |
|---:|---|---:|---|---:|---:|
| 1 | ζ(5) | 2 | above | 0.00307224485663 | 0.0768061214158 |
| 1 | ζ(7) | 3 | below | 0.000349277381923 | 0.0436596727404 |
| 2 | ζ(5) | 5 | above | 4.59614204948e-6 | 0.391885455707 |
| 2 | ζ(7) | 8 | below | 2.7992717117e-6 | 69.6937141422 |
| 3 | ζ(5) | 10 | above | 4.89903925954e-9 | 39.7715823572 |
| 3 | ζ(7) | 14 | below | 3.99452541643e-9 | 187436.850282 |
| 4 | ζ(5) | 15 | above | 4.73726617301e-12 | 641.886470652 |
| 4 | ζ(7) | 20 | below | 4.39995223469e-12 | 314793299.351 |
| 5 | ζ(5) | 21 | above | 4.39370870402e-15 | 1145857.52984 |
| 5 | ζ(7) | 29 | below | 4.39565864484e-15 | 3.75551815552e+14 |
| 6 | ζ(5) | 22 | above | 3.98756203765e-18 | 9520.62089034 |
| 6 | ζ(7) | 32 | below | 4.18601730022e-18 | 1.54432160644e+14 |
| 7 | ζ(5) | 30 | above | 3.57315857068e-21 | 2657236519.19 |
| 7 | ζ(7) | 44 | below | 3.87966984757e-21 | 4.24829428603e+22 |
| 8 | ζ(5) | 33 | above | 3.17595664989e-24 | 2372282049.2 |
| 8 | ζ(7) | 48 | below | 3.53555239892e-24 | 8.71430699001e+23 |
| 9 | ζ(5) | 41 | above | 2.80747140939e-27 | 3.13372584075e+13 |
| 9 | ζ(7) | 57 | below | 3.18599625843e-27 | 1.68011080311e+30 |
| 10 | ζ(5) | 43 | above | 2.47212748863e-30 | 2.3430521003e+13 |
| 10 | ζ(7) | 62 | below | 2.84856928973e-30 | 3.71678627571e+31 |

Exactly three of these twenty primitive errors are below 1: both index-1 values and the index-2 ζ(5) value. The other seventeen exceed 1. The tested ζ(5) coefficients are all above their target, and the tested ζ(7) coefficients all below; these are twenty individual finite sign certificates, not an all-index statement. There are substantial arithmetic dips (notably index 6), so no monotonic primitive-growth assertion follows from the table.

The parameter derivative reduces the local denominator power relative to differentiation in the sampling variable. At degree 10 the central base-measure sample had primitive errors about \(1.68\cdot10^{26}\) for ζ(5), \(1.83\cdot10^{60}\) for ζ(7); the present values are about \(2.34\cdot10^{13}\) and \(3.72\cdot10^{31}\). The separately optimized sampling locations of other families have not been optimized here: this construction is the diagonal parameter coefficient. The observed gains do not establish an irrationality argument or its impossibility.

## Exact denominator mechanism

Here is a general local lemma, independent of any coprimality guess. Let
\[
 A(x)=A_0+A_1x+\cdots,\qquad b(x)=b_0+b_1x+\cdots,
\]
with all coefficients through order \(r\) in \(\mathbb Z_p\). Suppose \(v_p(A_0)>0\) and \(A_1,b_0\) are units. In the coefficient of \(x^r\) in \(b/A\), the unique term with denominator \(A_0^{r+1}\) is
\[
 \frac{(-1)^r b_0 A_1^r}{A_0^{r+1}}.
\]
Every other term uses at most \(r\) powers of \(A_0\) in its denominator and has integral numerator. Thus its valuation is strictly greater, and
\[
 \boxed{v_p([x^r](b/A))=-(r+1)v_p(A_0).}
\]
This proves a denominator power \(r+1\), compared with \(2r+1\) for the base-measure sampling-variable derivative targeting the same ζ(2r+3).

For each experiment, construct the largest divisor \(G\) of the integer \(A_n(0)\) coprime to (i) all denominators of the actual computed \(A,b\) jets, (ii) the actual numerators of \(A_1,b_0\), and (iii) \(n!\). Repeated gcd removal constructs this divisor without factorization. The checker verifies that the complete part of the reduced denominator supported on primes dividing \(G\) is **exactly** \(G^{r+1}\), for all twenty rows. At degree 10, \(G=A_{10}(0)=13657436403073\). Thus noncancellation is established on that support, rather than inferred from the formal denominator. Prime support outside \(G\) is not classified here.

## Reproduction and remaining questions

Run `python codex/computations/explore_sigma_taylor_zeta.py`; full rational approximants, signed rational error intervals, primitive certificates and good denominator factors are stored in `codex/computations/sigma_taylor_zeta_experiment.json`. The script reuses only the independently documented exact zeta enclosure and integer gcd-removal helpers from the base experiment.

Useful next theory would be an analytic coefficient-error formula that explains the opposite first/second signs and proves all-index convergence, together with arithmetic control of the actual good factors. A positive real error for each parameter alone is insufficient for either conclusion about its Taylor coefficients. No more finite sampling is needed to establish the present conclusions.
