# One interpolation-and-balance recipe for the zeta(2) and zeta(3) recurrences

Status: elementary general identities proved below, with exact checks. Independent
review passed; see `zeta-balance-independent-review.md`. These are classical sequences and recurrence formulas; no
historical novelty is claimed. The new explanatory organization is not an account
of the authors' historical thought process.

The previous Newton-kernel argument explains17and5 in Apéry's zeta(3) recurrence.
A closely related kernel explains11and3 in the zeta(2) recurrence. Crucially, the
same finite balance calculation predicts the different irrational saddle locations.

## 1. Keep the interpolation basis, change its partner

Write
\[
\phi_k(n)=\binom nk\binom{n+k}k
=\frac{\prod_{j=0}^{k-1}[n(n+1)-j(j+1)]}{(k!)^2},
\qquad \psi_k(m)=\binom mk.
\]
The two elementary multiplication rules are
\[
n(n+1)\phi_k(n)=k(k+1)\phi_k(n)+(k+1)^2\phi_{k+1}(n),
\]
\[
(2m+1)\psi_k(m)=(2k+1)\psi_k(m)+2(k+1)\psi_{k+1}(m).
\tag{1}
\]
Define the mixed kernel
\[
H(n,m)=\sum_{k=0}^{\min(n,m)}\phi_k(n)\psi_k(m).
\tag{2}
\]
Unlike the zeta(3) kernel K(n,m), this one is not symmetric. Its diagonal is
\[
C_n=H(n,n)=\sum_{k=0}^n\binom nk^2\binom{n+k}k,
\quad C_0=1,C_1=3,C_2=19,C_3=147.
\]

To transpose the second multiplication rule in(1), seek a nearest-neighbor
operator on n with
\[
\mathcal B_n\phi_k=(2k+1)\phi_k+2k\phi_{k-1}.
\]
Requiring this for k=0,1 fixes its two coefficients:
\[
\boxed{\quad
\mathcal B_nf=
\frac{(n+1)^2f(n+1)-n^2f(n-1)}{2n+1}.
\quad}
\tag{3}
\]
The ratio identities for phi verify the displayed action for every k; polynomial
continuation handles the apparent zero denominators. The backward term is zero
at n=0. Transposing(1) inside the finite sum now proves
\[
(n+1)^2 H(n+1,m)=(2n+1)(2m+1)H(n,m)+n^2H(n-1,m).
\tag{4}
\]
The square powers arise from pairing the quadratic-grid basis with an ordinary
binomial basis. Pairing it with another quadratic-grid basis instead produced
the cubic operator in the zeta(3) argument.

For completeness, the other variable also has a simple equation:
\[
(m+1)^2H(n,m+1)+m^2H(n,m-1)
=[2m^2+2m+1+n(n+1)]H(n,m).
\tag{5}
\]
Indeed the operator
\((m+1)^2\Delta_++m^2\Delta_-\) sends
\(\psi_k\) to \(k(k+1)\psi_k+k^2\psi_{k-1}\); transpose the first rule in(1).
Equations(4)--(5) are two-variable statements. A diagonal recurrence needs an
additional identity, not the invalid substitution H(n+1,n)=C_(n+1).

## 2. A finite balance identity closes the diagonal

Fix an integer N>=1 and put
\[
w_k=\binom Nk^2\binom{N+k-1}k,\qquad0\le k\le N.
\]
Its consecutive ratio gives
\[
w_k(N-k)^2(N+k)=w_{k+1}(k+1)^3.
\tag{6}
\]
The last birth factor vanishes at k=N and the first death factor at k=0.
Sum(6) and shift the finite sum by one. Because
\((N-k)^2(N+k)-k^3=N(N^2-Nk-k^2)\), one gets
\[
\boxed{\quad\sum_{k=0}^N(k^2+Nk)w_k=N^2\sum_{k=0}^Nw_k.\quad}
\tag{7}
\]
This can be read as detailed balance for a finite birth-and-death chain, but it
is simply an exact finite summation. No probability theorem is being invoked.

Now let U=H(N,N-1), V=H(N-1,N). Relative to the common weights w_k, their
summands and those of C_N are
\[
C_N:\ 1+\frac kN,\qquad
U:\ 1-\frac{k^2}{N^2},\qquad
V:\ 1-\frac kN.
\]
Therefore(7) gives
\[
C_N=2U+V.
\tag{8}
\]
A second identity is termwise. With common base
\(\phi_k(N-1)\psi_k(N)\), the ratios for H(N,N-1) and C_(N-1) add to2.
Equivalently, direct binomial substitution gives
\[
C_{N-1}=2V-U.
\tag{9}
\]
All identities include the endpoints: factors that appear to divide by N-k in
a ratio argument are canceled in the displayed polynomial weights first.
Solving(8)--(9),
\[
H(N,N-1)=\frac{2C_N-C_{N-1}}5,\qquad
H(N-1,N)=\frac{C_N+2C_{N-1}}5.
\tag{10}
\]
The integer5 is thus the determinant of a concrete adjacent-value relation.

## 3. The coefficients11and3 are now forced

Set m=n in(4) and insert(10) at N=n+1 and N=n:
\[
(n+1)^2(2C_{n+1}-C_n)
=5(2n+1)^2C_n+n^2(C_n+2C_{n-1}).
\]
Collecting terms proves
\[
\boxed{\quad
(n+1)^2C_{n+1}=(11n^2+11n+3)C_n+n^2C_{n-1}.
\quad}
\tag{11}
\]
Specifically,
\[
5(2n+1)^2+(n+1)^2+n^2=22n^2+22n+6.
\]
There is no guessed coefficient here: the multiplication rule, its transposed
operator, and the finite balance identity determine it.

Compare zeta(3): the symmetric kernel has adjacent identity
\(A_n+A_{n-1}=6K(n,n-1)\). Substitution in its cubic eigenvalue equation yields
\[
(n+1)^3A_{n+1}=(2n+1)(17n(n+1)+5)A_n-n^3A_{n-1},
\]
where17=6*3-1 and5=6*1-1. The sign difference between the backward terms is
also inherited from the two operators, rather than chosen to fit a desired error.

Identifying the companion limits with zeta(2) or zeta(3), controlling their
rational denominators and proving nonvanishing are separate steps. This note
explains the recurrences for the integer denominator sequences; it does not
silently replace those remaining irrationality obligations.

## 4. The golden-ratio and square-root-of-two scales have the same origin

For the zeta(2) edge weights, the leading balance polynomial in k/N is
\[
1-x-x^2.
\]
Its positive root is \((\sqrt5-1)/2\). For the zeta(3) squared edge weights
\(w_k=[\binom Nk\binom{N+k-1}k]^2\), detailed balance instead gives exactly
\[
\sum k^2w_k=\frac{N^2}2\sum w_k,
\]
with leading polynomial1-2x² and positive root1/sqrt2.

The mode statement can be made without claiming that a moment identity proves
concentration. The consecutive ratios are strictly decreasing in k: the logarithmic
derivative of \((N-k)^r(N+k)^s/(k+1)^{r+s}\) is
\(-r/(N-k)+s/(N+k)-(r+s)/(k+1)<0\) when s<=r, which covers the three cases below.
The unique crossing of the ratio through1 therefore locates the mode. Scaling
k/N in(6), or its squared counterpart, shows these mode locations approach the
stated roots, with O(1) unscaled error by the nonzero derivative at the crossing.
A Gaussian concentration theorem would require an additional argument and is not
inferred merely from(7).

The recurrence growth rates are correspondingly
\[
\frac{11+5\sqrt5}{2}=\left(\frac{1+\sqrt5}2\right)^5
\quad\text{and}\quad17+12\sqrt2=(1+\sqrt2)^4.
\]
A sharp asymptotic constant needs more work, but the quadratic fields themselves
are already visible in elementary balance equations.

## 5. Why these low-order examples are unusually simple

Consider the general positive integer weight family
\[
w_k^{r,s}=\binom Nk^r\binom{N+k-1}k^s,\qquad r,s\ge1.
\]
The same finite shift gives
\[
\sum_k w_k^{r,s}\big[(N-k)^r(N+k)^s-k^{r+s}\big]=0.
\tag{12}
\]
For positive integer r,s, the polynomial in brackets has degree at most2 in k
**exactly for**
\[
(r,s)=(1,1),(2,1),(2,2).
\tag{13}
\]
Proof: if r is odd, the leading coefficient is-2, so degree=r+s; only(1,1)
works. If r is even, the leading powers cancel. The next coefficient is
(s-r)N. If s!=r, degree=r+s-1, so only(2,1) works. If s=r, expansion of
\((N^2-k^2)^r-k^{2r}\), with even r, has degree2r-2 and nonzero coefficient;
only r=2 works.

These three cases connect to the familiar denominator sequences
\[
\sum_k\binom Nk\binom{N+k}k=P_N(3),\quad
\sum_k\binom Nk^2\binom{N+k}k=C_N,\quad
\sum_k\binom Nk^2\binom{N+k}k^2=A_N.
\]
The first is the Legendre family for log2, developed in recurrence-design.md;
the other two are the zeta(2) and zeta(3) families. This is a classification only
of the degree of the single balance polynomial(12). It does not classify all
short recurrences, all useful binomial sums, or all irrationality proofs.
Higher powers can still have useful structure, but this particular two-moment
closure no longer follows automatically.

## 6. The balancing identity is an actual insertion bijection

The weights in(12) count concrete finite objects. An object at level k consists
of r subsets of size k of an N-element set, together with s weak compositions
of k into N parts. Represent each composition by k stars and N-1 bars.

To move upwards, insert one new element into each subset and one new star into
each stars-and-bars word. There are N-k choices for each subset and N+k gaps
in each word. Mark every newly inserted element and star. In the reverse
direction, an object at level k+1 has k+1 choices of a marked element in each
subset and k+1 choices of a marked star in each word; erase all the marks and
their selected items. These operations are inverse, including when several
insertions yield the same unmarked word: the marked star retains the chosen gap.
Thus counting both ends of the same edges gives exactly
\[
w_k^{r,s}(N-k)^r(N+k)^s=w_{k+1}^{r,s}(k+1)^{r+s}.
\]
Summing cancels the edge counts and proves(12). The birth-and-death terminology
is therefore optional; the underlying mechanism is a reversible insertion.

Why use the edge weights with N+k-1, rather than the N+k in the Apéry sums?
They provide a common counting basis for neighboring values of the kernels.
Changing from weak compositions into N parts to N+1 parts multiplies the count
by (N+k)/N. This is precisely the small factor that appears in the adjacent-value
identities above and in the zeta(3) closure. The common weight is chosen to make
those neighboring formulas comparable, not introduced independently of them.

Reproduce the arbitrary-index basis identities, finite balance, adjacent closure,
and recurrence with computations/check_zeta_balance_unification.py.
