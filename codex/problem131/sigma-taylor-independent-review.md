# Independent review: uniform parameter errors and exact Taylor arithmetic

**Verdict:** the analytic argument in
[sigma-taylor-error-asymptotic.md](sigma-taylor-error-asymptotic.md) is accepted
without correction, using its stated common recurrence. The all-index
associated-numerator recurrence is a separate input from
[deformed-associated-pearson-symmetry.md](deformed-associated-pearson-symmetry.md);
its [independent audit](deformed-associated-pearson-independent-review.md) is now complete and accepted. This review does not
replace that audit with finite recurrence checks.

The local denominator lemma and all twenty saved finite certificates in
[sigma-taylor-zeta-experiment.md](sigma-taylor-zeta-experiment.md) are also
accepted. These conclusions give no arithmetic asymptotic or new odd-zeta
irrationality theorem. The denominator family and its recurrence have the
exact classical Askey–Wilson antecedent documented in the source dictionary;
historical priority of the coefficient asymptotic is not asserted.

## 1. The Casoratian has constant one in this normalization

Write \(q_n=b_n/A_n\), \(E_n=Z-q_n\), and
\[
 \alpha_n(\sigma)=
 \frac{(n+1)((n+1)^2-\sigma)^2}{2(3(n+1)^2-\sigma)}.
\]
The common recurrence has backward coefficient \(\alpha_{n-1}\),
not a differently normalized coefficient. Consequently, for
\(W_n=b_nA_{n-1}-b_{n-1}A_n\),
\[
 \alpha_nW_{n+1}=\alpha_{n-1}W_n.
\]
The initial data give
\(W_1=b_1=2(3-\sigma)/(1-\sigma)^2=1/\alpha_0\).
Thus \(W_n=1/\alpha_{n-1}\), and
\[
 q_k-q_{k-1}=h_k:=
 \frac1{\alpha_{k-1}A_kA_{k-1}}.
\]
At \(\sigma=0\), \(\alpha_{k-1}=k^3/6\), recovering the ordinary
Apéry Casoratian. Both the factor six and the index shift are correct.

For real \(0\le\sigma\le1/4\), the previously proved Padé estimate
\(0<E_n\le2/((2n+1)^3A_n^2)\) implies convergence to \(Z\).
It follows on that real interval that \(E_n=\sum_{k>n}h_k\).
The extension to complex parameters needs the following separate argument.

## 2. A genuinely common zero-free closed disk

Put \(w_k(\sigma)=\prod_{j=1}^k(1-\sigma/j^2)^{-1}\) and
\(G(\sigma)=\prod_{j\ge1}(1-\sigma/j^2)^{-1}\). Product convergence
is uniform on compact disks of radius below one, since \(\sum j^{-2}\)
converges. Euler's product gives \(G(1/4)=\pi/2\).

Each \(w_k\) has nonnegative Taylor coefficients and constant term one.
For \(|\sigma|\le1/4\), therefore,
\[
 |w_k(\sigma)-1|\le w_k(1/4)-1\le\pi/2-1=:c_0<1.
\]
The quotient \(A_n(\sigma)/A_n(0)\) is a convex combination of these
\(w_k\), with weights \(C_{nk}^2/A_n(0)\). Hence
\[
 \left|A_n(\sigma)/A_n(0)-1\right|\le c_0,
 \qquad |A_n(\sigma)|\ge(1-c_0)A_n(0)>0
\]
for **every** \(n\), on the same closed disk. This is stronger than
eventual nonvanishing and suffices for all subsequent divisions.
The finite quotient definition of \(b_n\) has parameter poles only among
positive integer squares, so it is holomorphic near this disk. The recurrence
also gives this directly from its rational coefficients and initial data.
Thus \(q_n\) and \(E_n\) are holomorphic inside and continuous on the disk.

For completeness, a uniform majorant of the tail does not require a precise
Apéry asymptotic. Its \(k=n\) summand gives
\[
 A_n(0)\ge\binom{2n}{n}^2\ge\frac{16^n}{(2n+1)^2}.
\]
Moreover
\[
 \left|\frac{\alpha_{k-1}(0)}{\alpha_{k-1}(\sigma)}\right|
 =\left|\frac{1-\sigma/(3k^2)}{(1-\sigma/k^2)^2}\right|
 \le\frac{52}{27}.
\]
It follows that \(|h_k(\sigma)|\le
\frac{52}{27}(1-c_0)^{-2}h_k(0)\). The preceding exponential lower
bound on \(A_n(0)\) makes \(\sum h_k(0)\) converge. The Weierstrass
test proves uniform convergence of the tails on the closed disk. Their
holomorphic sum agrees with \(E_n\) on a real interval; the identity theorem
and boundary continuity extend the exact tail identity to the entire disk.

## 3. Finite-head escape proves the full uniform multiplier

For a fixed \(M\), \(\sum_{k\le M}C_{nk}^2\) grows polynomially in
\(n\), whereas \(A_n(0)\) has the exponential lower bound above.
Thus the convex weight of this fixed head tends to zero. On the remaining
indices, \(w_k\) is uniformly close to \(G\). Choosing \(M\) first
and \(n\) second rigorously gives
\[
 A_n(\sigma)/A_n(0)\longrightarrow G(\sigma)
 \quad\text{uniformly on }|\sigma|\le1/4.
\]
The same nonvanishing bound permits taking inverses. Together with
\(\alpha_{k-1}(0)/\alpha_{k-1}(\sigma)\to1\) uniformly, this proves
\(h_k(\sigma)/h_k(0)\to G(\sigma)^{-2}\) uniformly.

Finally the normalized tail is itself a convex combination, because
\(h_k(0)>0\):
\[
 \frac{E_n(\sigma)}{E_n(0)}
 =\sum_{k>n}\frac{h_k(0)}{\sum_{j>n}h_j(0)}
               \frac{h_k(\sigma)}{h_k(0)}.
\]
Its uniform error is bounded by the supremum of the summand errors over
\(k>n\). This tends to zero. There is no unproved interchange of an
asymptotic expansion with a tail, and no saddle-point estimate is needed.

## 4. Coefficient extraction and eventual signs

For a fixed circle \(|\sigma|=\rho<1/4\), Cauchy's formula bounds
the error of the coefficient of order \(r\) by \(\rho^{-r}\) times
the uniform error on that circle. Hence
\[
 \frac{\zeta(2r+3)-[\sigma^r]q_n}{E_n(0)}
 \longrightarrow
 [\sigma^r]\left(\frac{\sin(\pi\sqrt\sigma)}
                         {\pi\sqrt\sigma}\right)^2
 =\frac{(-1)^r2^{2r+1}\pi^{2r}}{(2r+2)!}.
\]
The apparent square root is harmless: the quotient squared is an entire
function of \(\sigma\). All the displayed coefficients are nonzero.
Since \(E_n(0)>0\), the error has eventual sign \((-1)^r\) for each
fixed \(r\). In particular the first coefficient is eventually above
\(\zeta(5)\) and the second below \(\zeta(7)\). This does not prove
those signs for every finite \(n\), nor a threshold uniform in \(r\),
nor monotonicity of the higher coefficient approximants.

## 5. The local denominator lemma is exact under its unit hypotheses

Let the jets of \(A,b\) through order \(r\) lie in \(\mathbb Z_p\),
let \(e=v_p(A_0)>0\), and suppose \(A_1,b_0\) are units. Formally expand
\[
 \frac bA=\frac b{A_0}\sum_{j=0}^r(-1)^j
       \left(\frac{A-A_0}{A_0}\right)^j+O(x^{r+1}).
\]
To contribute to \(x^r\) with denominator power \(r+1\), the unique
possibility is \(j=r\), the constant coefficient \(b_0\), and \(r\)
copies of \(A_1x\). Its valuation is \(-(r+1)e\). Every other term
has integral numerator and denominator power at most \(r\), so valuation
at least \(-re\), strictly larger. Therefore
\[
 v_p([x^r](b/A))=-(r+1)e.
\]
This includes \(r=0\), for which there are no competing terms.
Neither a coprimality guess nor generic noncancellation is needed.
The hypotheses are essential: the assertion does not classify primes
where a jet is nonintegral or either designated numerator is not a unit.

Removing repeatedly from \(A_n(0)\) its gcd with the recorded bad factor
leaves exactly its largest divisor \(G\) supported away from those bad
primes. Every retained prime satisfies the lemma's hypotheses. Consequently
the part of the reduced denominator supported on \(G\) is precisely
\(G^{r+1}\). Removing \(n!\) as well is an additional conservative
restriction, not an assumption about all other primes. There is no theorem
here that \(G=A_n(0)\) for every \(n\), or that \(G\) has any particular
asymptotic growth. At \(n=10\), the saved exact value
\(G=A_{10}(0)=13657436403073\) is verified.

## 6. Independent reconstruction of every finite certificate

[check_sigma_taylor_independent.py](../computations/check_sigma_taylor_independent.py)
uses only Python's standard library and imports none of the author's code.
It reconstructs the numerator by an independently arranged Newton quotient:
with \(d_j(\sigma)=\prod_{a=1}^j(a^2-\sigma)\),
\[
 N_n(\sigma)=\sum_{k=0}^n C_{nk}^2w_k(\sigma)
    \sum_{m=1}^k(-1)^{m-1}
       \frac{d_{m-1}(\sigma)(n-m)!}{m(n+m)!},
 \qquad b_n=A_nT_n+N_n/2.
\]
This follows directly by the finite divided-product identity and
\(\nu_\sigma(D_j)=d_j/(j+1)\). It avoids both generating \(b_n\)
by the recurrence and reconstructing the full polynomial in \(t\).
Exact jets through order four check 14 constructed pairs, 24 recurrence
identities, and 12 Casoratians. These finite checks corroborate, and do not
prove, the all-index recurrence input.

For the saved twenty rows, the checker reproduces every rational approximant,
good factor, forced denominator component, digit count, error endpoint, and
primitive-error endpoint. Its Bernoulli numbers are generated independently
by the defining rational recurrence. The Euler–Maclaurin enclosure uses
\(N=80,L=60\). The remainder bound is valid because
\[
 \sup_{0\le x\le1}|B_{2L}(x)|
 \le\sum_{j=0}^{2L}\binom{2L}{j}|B_j|,
\]
and integration of the absolute \(2L\)-th derivative of \(x^{-s}\)
gives the radius
\[
 \frac{\sum_j\binom{2L}{j}|B_j|}{(2L)!}
     \frac{(s)_{2L-1}}{N^{s+2L-1}}.
\]
Thus no numerical zeta evaluation certifies the saved signs. All twenty
exact enclosures pass: precisely \((n,r)=(1,1),(1,2),(2,1)\) have
primitive errors below one; all other seventeen exceed one. The checker
also verifies 240 local prime-valuation examples through order four.

The author's asymptotic checker was rerun successfully. Its high-precision
real and complex samples are useful corroboration only; the uniform limit
and eventual signs rest on §§1–4 above. The finite arithmetic table exhibits
substantial dips and supplies no monotonic primitive-growth assertion or
asymptotic obstruction. No further finite sampling is needed to establish
these twenty certificates.
