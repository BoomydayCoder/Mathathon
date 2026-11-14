# A short adjacent-balance telescope for the deformed Apéry diagonal

Status: the finite adjacent-balance identity is proved below. A self-contained Legendre-operator derivation of the grid relation then proves the stated diagonal recurrence. The associated numerator now has an all-index proof in [the Pearson symmetry note](deformed-associated-pearson-symmetry.md), accepted in its [independent audit](deformed-associated-pearson-independent-review.md); its proof is separate from the denominator argument here. No historical novelty claim is made. The certificate was supplied by the parent researcher and independently verified here.

**Exact classical antecedent.** The full denominator family and recurrence
are the specialization a=d=0, b=-tau, c=tau, e=tau, f=-tau of Askey and
Wilson, *A recurrence relation generalizing those of Apéry* (1984),
equations (22)–(23), printed pp.271–273. Their F_n is exactly our A_n;
their three coefficients are (alpha_n/2,-beta_n/2,gamma_n/2).
The [primary-source dictionary](koecher-askey-wilson-primary-dictionary.md)
records the full symbolic substitution, original-page inspection and PDF hash.
The finite adjacent-balance proof below is independently reconstructed.
This exact denominator match does not by itself prove the associated
numerator recurrence.


Write
\[
 C_{nk}=\binom nk\binom{n+k}k,\quad
 w_k=\frac{(k!)^2}{\prod_{j=1}^k(j^2-\sigma)},\quad
 K(n,m)=\sum_{k=0}^{\min(n,m)}C_{nk}C_{mk}w_k,
 \quad A_n=K(n,n).
\]
All statements initially hold as rational-function identities in an indeterminate \(\sigma\). They therefore hold at every specialization at which their displayed denominators are defined. In particular they apply throughout \(0\le\sigma\le1/4\).

## Exact finite telescope

For every \(n\ge1\),
\[
 \boxed{A_n+A_{n-1}=c_nK(n,n-1),\qquad
 c_n=\frac{2(3n^2-\sigma)}{n^2-\sigma}.} \tag{1}
\]
Put \(T_k=C_{nk}^2w_k\). For \(0\le k\le n\),
\[
 \frac{C_{n-1,k}}{C_{nk}}=\frac{n-k}{n+k},\qquad
 \frac{T_{k+1}}{T_k}
 =\frac{(n-k)^2(n+k+1)^2}{(k+1)^2((k+1)^2-\sigma)}.
\]
Define
\[
 R_k=-\frac{4k^2(k^2-\sigma)}{(n^2-\sigma)(n+k)^2}.
\]
Then direct cancellation gives
\[
 T_{k+1}R_{k+1}-T_kR_k
 =T_k\left[1+\left(\frac{n-k}{n+k}\right)^2
       -c_n\frac{n-k}{n+k}\right]. \tag{2}
\]
Indeed the bracket equals
\[
 \frac{4[(2n^2-\sigma)k^2-n^4]}
      {(n^2-\sigma)(n+k)^2},
\]
and the needed polynomial identity is only
\[
 (n^2-k^2)^2-k^2(k^2-\sigma)
 =n^4-(2n^2-\sigma)k^2.
\]
Sum (2) from \(k=0\) to \(n\). Both endpoints vanish: \(R_0=0\), and the terminating summand gives \(T_{n+1}=0\). This proves (1), including the last summand where \(C_{n-1,n}=0\).

## Why the diagonal recurrence has these coefficients

The ordinary-grid relation is
\[
 u_nK(n+1,m)-(2n+1)[n(n+1)+2m(m+1)+1-\sigma]K(n,m)
 +v_nK(n-1,m)=0,
\]
where \(u_n=(n+1)((n+1)^2-\sigma)\) and \(v_n=n(n^2-\sigma)\). Here is a self-contained derivation. Set
\[
 Q_k(t)=\frac{\prod_{r=0}^{k-1}(t+(r+1/2)^2)}
 {\prod_{r=1}^k(r^2-\sigma)},\qquad U(u^k)=Q_k(t).
\]
This is an invertible degree-preserving linear map over \(\mathbb Q(\sigma)\). The Newton products directly give
\[
 tQ_k=((k+1)^2-\sigma)Q_{k+1}-(k+1/2)^2Q_k. \tag{4}
\]
Define \(J=-\frac d{du}u(1-u)\frac d{du}\), and write \(u\) also for multiplication by \(u\). The operator
\[
 T_\sigma=\tfrac12(uJ+Ju)-\sigma u+\tfrac14
\]
satisfies
\[
 Ju^k=k(k+1)u^k-k^2u^{k-1},\qquad
 T_\sigma u^k=((k+1)^2-\sigma)u^{k+1}-(k+1/2)^2u^k.
\]
Thus (4) proves \(UT_\sigma=tU\), without requiring an integral representation. For shifted Legendre polynomials
\[
 L_n(u)=P_n(1-2u)=\sum_{k=0}^n(-1)^kC_{nk}u^k,
\]
the differential equation and three-term relation give
\[
 JL_n=n(n+1)L_n,\qquad
 uL_n=\tfrac12L_n-\frac{n+1}{2(2n+1)}L_{n+1}
                    -\frac n{2(2n+1)}L_{n-1}.
\]
Apply \(\tfrac12(uJ+Ju)-\sigma u+\tfrac14\) and then \(U\). Since \(UL_n=\Phi_n\), this yields
\[
 t\Phi_n=-\frac{u_n}{2(2n+1)}\Phi_{n+1}
 +\left[\frac{n(n+1)}2+\frac14-\frac\sigma2\right]\Phi_n
 -\frac{v_n}{2(2n+1)}\Phi_{n-1}. \tag{5}
\]
Finally \(\Phi_n(-(m+1/2)^2)=K(n,m)\): the Newton product at this point is \((-1)^k(k!)^2C_{mk}\), and vanishes for \(k>m\). Substituting this value of \(t\) into (5) proves the grid relation. This derivation includes \(n=0\), with the previous-index term defined to be zero.

Set \(m=n\), use symmetry of \(K\), and substitute (1) at \(n,n+1\). This proves
\[
 \boxed{\alpha_nA_{n+1}-\beta_nA_n+\gamma_nA_{n-1}=0,}\tag{3}
\]
\[
 \alpha_n=\frac{(n+1)((n+1)^2-\sigma)^2}{2(3(n+1)^2-\sigma)},\qquad
 \gamma_n=\frac{n(n^2-\sigma)^2}{2(3n^2-\sigma)},
\]
\[
 \beta_n=(2n+1)[3n(n+1)+1-\sigma]-\alpha_n-\gamma_n.
\]
Thus the diagonal coefficients arise by dividing the two grid-edge coefficients by their adjacent-balance factors, rather than by fitting diagonal values. At \(\sigma=0\), multiplying (3) by 6 gives exactly
\[
 (n+1)^3A_{n+1}-(2n+1)(17n^2+17n+5)A_n+n^3A_{n-1}=0.
\]
At \(\sigma=1/4\), (1) and (3) remain regular and apply to the logarithmic endpoint of the positive-measure deformation. This observation alone supplies no irrationality claim for that endpoint.

## Associated numerator: evidence and remaining proof obligation

Use the precise numerator from [the positive deformation](logarithm-zeta3-wilson-deformation.md), with \(\sigma=\tau^2\):
\[
 \Phi_n(t)=\sum_{k=0}^n(-1)^kC_{nk}
 \frac{\prod_{r=0}^{k-1}(t+(r+1/2)^2)}{\prod_{r=1}^k(r^2-\sigma)},
\]
\[
 b_n=A_n\sum_{r=1}^n\frac1{r(r^2-\sigma)}+
 \frac12\nu_\tau\!\left(\frac{A_n-\Phi_n(t)}{t+(n+1/2)^2}\right).
\]
Polynomial moments can be computed without adjoining \(\tau\):
\[
 \nu_\tau(t^j)=\frac{(-1)^j}{j+1}
 \sum_{h=1}^{j+1}\binom{2j+2}{2h}
 B_{2j+2-2h}(1/2)\sigma^{h-1}.
\]
Exact generic calculations give
\[
 b_0=0,\quad b_1=\frac{2(3-\sigma)}{(1-\sigma)^2},\quad
 b_2=-\frac{3(\sigma^3-32\sigma^2+247\sigma-468)}
 {(\sigma-4)^2(\sigma-1)^2}.
\]
The residual \(\alpha_nb_{n+1}-\beta_nb_n+\gamma_nb_{n-1}\) vanishes identically for \(n=1,2,3\), using independently constructed polynomial moments and difference quotients. This is finite evidence only. A proof should establish an adjacent balance for the corresponding **mixed linear forms**, not infer it from the denominator balance. The numerator cannot satisfy the homogeneous relation at \(n=0\), since \(b_0=0\) and \(b_1\ne0\); the intended range is \(n\ge1\).

Reproduce with `python codex/computations/check_generic_sigma_diagonal_recurrence.py`. It checks the telescope symbolically with free \(n,k,\sigma\), six generic diagonal balances and recurrences, twenty grid relations, and the three numerator residuals above.
