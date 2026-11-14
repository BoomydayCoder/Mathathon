# Actual scalar-normalization obstruction for every fixed single depletion

Status: proved local coefficient-denominator bounds and their consequence
for the canonical determinant, subject to the separately audited monic
Jensen lower bound cited below. This develops the parent's proposed
normalized Euler obstruction and its extension to all fixed integer
depletions. It concerns this determinant family, not irrationality proofs
in general. No historical novelty is asserted.

Fix odd \(s\ge7\), a fixed integer \(L\ge2\), and
\[
K=40n,\quad N=3n,\quad h=K-N,\quad m=s-1,\quad C=1-L^m.
\]
Use the normalized positive density
\[
\sigma_{L,s}(y)=w_s(y)-L^{-1}w_s(Ly)>0.
\]
It deletes the frequencies divisible by \(L\) from the original positive
exponential sum and retains the real decay exponent \(2\pi\).
This is the normalized, integer-node model, not the unscaled alternating
integer-node model.

## 1. Moments, poles and monic basis

The polynomial moments satisfy
\[
\Sigma_{L,s}(t^e)=(1-L^{-2e-2})\mu_s(t^e)
\in\mathbb Z[1/L].
\]
For example, if \(\kappa_k\) is the \(k\)-th logarithmic derivative at
zero of
\((1+e^z+\cdots+e^{(L-1)z})/L\), then
\[
\Sigma_{L,s}(t^e)
=(-1)^eL^{-2e-2}\kappa_{2e+2}
                    \binom{2e+s}{s-1}.
\]
The cumulant recurrence proves \(\kappa_k\in\mathbb Z[1/L]\).
Thus every polynomial moment is integral at \(p\nmid L\), with no
additional condition \(p>s\).

The exact pole value is
\[
\Sigma_{L,s,X}((t+j^2)^{-1})=Cj^mX+c_j,
\qquad
c_j=j^m(L^mH_{Lj}^{(s)}-H_j^{(s)})
                         +\frac{1-L^{-1}}{2j}.          \tag{1}
\]
For \(L=2\), the last term is \(1/(4j)\).

Put \(D=D_K/D_N\), \(W=D_N^5\), and let \(\Delta(X)\) be the
source-like determinant of \(\Sigma_{L,s,X}(Wt^{i+j}/D)\),
\(0\le i,j<h\). Its leading coefficient is
\[
\ell_K=(-1)^{h(h-1)/2}C^h
\prod_{j=N+1}^Kj^mD_N(-j^2)^5\ne0.
\]
Set \(F_K=\Delta/\ell_K\), so \(F_K\) is monic.
Use the integer basis \(f_j=D/(t+j^2)\), \(N<j\le K\).
Its change-of-basis determinant need not be a \(p\)-adic unit: it
cancels between \(\Delta\) and its leading coefficient.

Writing
\[
R_j=W(-j^2)D'(-j^2),\qquad C_j=Cj^mR_j,
\]
the coefficient matrix of \(X\) is diagonal. Off-diagonal integrands
are integer polynomials, and diagonal integrands have the form
\[
Wf_j^2/D=Q_{jj}+R_j/(t+j^2),\qquad Q_{jj}\in\mathbb Z[t].
\]
With \(P_{jk}=\Sigma_{L,s}(Q_{jk})\), the constant matrix is
\(B=P+\operatorname{diag}(R_jc_j)\), and exactly
\[
\boxed{F_K(X)=\det(XI+\mathcal C^{-1}B),\quad
\mathcal C=\operatorname{diag}(C_j).}                   \tag{2}
\]

## 2. Fixed prime exceptions and the diagonal pole residues

For \(1\le k\le2L-1\), put
\[
u_k=L^mH_k^{(s)}-H_{\lfloor k/L\rfloor}^{(s)},\qquad
u_*=L^mH_L^{(s)}-\frac{1+L^{-1}}2.
\]
All these rational numbers are strictly positive:
\[
u_k\ge L^m-\zeta(s)>0,\qquad
u_*\ge L^m-\frac{1+L^{-1}}2>0.
\]
Exclude the finite set of primes dividing \(2LC\), or the numerator
or denominator of any \(u_k,u_*\).
Take \(K/2<p\le LK\) outside this set, with \(p^2>2K\).
These conditions hold throughout the band eventually, except for the
fixed exclusions.

For an ordinary node \(p\nmid j\), separating multiples of \(p\) in
the harmonic sums gives
\[
\overline{\frac{p^sc_j}{Cj^m}}
=\epsilon_j
=\frac{L^mH_{\lfloor Lj/p\rfloor}^{(s)}
                    -H_{\lfloor j/p\rfloor}^{(s)}}C
\pmod p.                                               \tag{3}
\]
Here \(j/p<2\). The residue is zero for \(Lj<p\), and nonzero for
\(Lj>p\), by the fixed exclusions. If desired one may impose \(p>2L\)
in the finite statement; then all the displayed harmonic denominators
are automatically units. Without that extra condition, their denominators
are already accounted for by the fixed exceptional set.

When \(p\le K\), the unique zero-class node is \(j=p\). Its scaled
residue is instead
\[
\boxed{\epsilon_p=
\frac{L^mH_L^{(s)}-(1+L^{-1})/2}{C}\ne0\pmod p.}         \tag{4}
\]
The final term of (1) is essential in (4).
For \(L=2\), the three possible nonzero ordinary residues and the zero
residue are
\[
\frac{2^m}{C},\quad
\frac{2^m-1/2}{C},\quad
\frac{2^m-1/2+2^m3^{-s}}{C},\qquad
\epsilon_p=\frac{2^m-1/4}{C}.
\]
In particular the zero-class exclusion includes the fixed numerator
\(2^{s+1}-1\).

The number of nonzero diagonal residues is exactly
\[
b_p=K-\max(N,\lfloor p/L\rfloor).                        \tag{5}
\]
It includes the zero-class node whenever that node is present.

## 3. Valuations of the diagonal leading factors

Every \(P_{jk}\) is \(p\)-integral. In an ordinary squared residue class,
let \(\ell\) count all its original nodes in \(1,\dots,K\), and let
\(\delta\in\{0,1\}\) record its removed node in \(1,\dots,N\).
The canonical proportion ensures \(2N<p\), so there is at most one
removed node. Since \(p^2>2K\), each same-class ordinary node
difference has valuation exactly one. Consequently
\[
v_p(C_j)=\ell+4\delta-1.                                \tag{6}
\]
For \(p>K\), \(\ell\le2\), so this is at most five.
For \(K/2<p\le K\), \(\ell\le4\), so it is at most seven.
At the unique zero node \(j=p\), every other squared-node difference
and \(W(-p^2)\) is a unit; hence
\[
v_p(C_p)=m=s-1.
\]

Set \(A=p^s\mathcal C^{-1}B\). At any row with \(v_p(C_j)<s\),
its polynomial-part entries vanish modulo \(p\). Its only surviving
entry is its diagonal pole residue (3) or (4).

## 4. Clean coefficient bounds for \(s\ge9\), and for \(p>K\)

For every odd \(s\ge9\), (6) is strictly less than \(s\) everywhere.
For \(s=7\) the same is true throughout \(p>K\).
Thus \(A\) is integral and reduces to the full diagonal matrix
\(\operatorname{diag}(\epsilon_j)\). Coefficientwise,
\[
\overline{p^{sh}F_K(Z/p^s)}
=\prod_{j=N+1}^K(Z+\epsilon_j).
\]
The coefficient of \(Z^{h-b_p}\) is the product of all nonzero
\(\epsilon_j\), hence is a unit. Therefore
\[
\boxed{v_p([X^{h-b_p}]F_K)=-sb_p.}                       \tag{7}
\]
This proves an actual monic coefficient denominator. It is not a
claim that all other coefficients have no larger denominators.

## 5. The exceptional four-node classes at \(s=7\)

In \(K/2<p\le K\), equality \(v_p(C_j)=7\) occurs only in removed
classes containing four original nodes. Such a class has removed
index \(a\) with
\[
2p-K\le a\le N.
\]
Thus there are
\[
d_p=(K+N-2p+1)_+
\]
such classes. Each has the three surviving nodes
\(p-a,p+a,2p-a\). Let \(J\) be their union; \(|J|=3d_p\).
Every node of \(J\) exceeds \(p/2\), hence exceeds \(p/L\), so all
these exceptional rows are counted among the \(b_p\) bad rows.

The entire matrix \(A\) is still integral. Rows outside \(J\) reduce
diagonally; rows inside \(J\) may couple to every column.
Ordering the outside rows first nevertheless gives a block-triangular
reduction of \(ZI+A\). Its determinant is the outside diagonal product
times a monic residual characteristic polynomial of degree \(3d_p\).
No nonsingularity of that residual block is required.

The outside diagonal product has \(b_p-3d_p\) nonzero roots, counted
with multiplicity. Therefore the full reduced characteristic polynomial
has at least that many nonzero roots. Some coefficient with descending
index \(r\ge b_p-3d_p\) must be nonzero. By (2), its valuation in
\(F_K\) is exactly \(-7r\). If \(M_K\) denotes the coefficient-denominator
lcm of the monic \(F_K\), this proves
\[
\boxed{v_p(M_K)\ge7\max(0,b_p-3d_p),\qquad s=7.}         \tag{8}
\]

One can avoid root terminology altogether: if \(t_0\) is the order of
vanishing at \(Z=0\) of the monic residual characteristic polynomial,
then \(0\le t_0\le3d_p\). The full product has order \(h-b_p+t_0\),
so its coefficient of that smallest nonzero power is a unit. The exact
detected descending index is \(r=b_p-t_0\).

This argument includes all couplings to the complement. Counting only
an exceptional self-block without using the outside diagonal rows would
not have justified (8).

## 6. Prime-sum cost and the canonical obstruction

For \(\alpha=N/K=3/40\), define
\[
A_L(\alpha)=\int_{1/2}^{L}[1-\max(\alpha,u/L)]\,du.
\]
For \(L\ge2\), this is at least \(A_2(3/40)=9/16\):
on \([1/2,2]\), increasing \(L\) can only increase the integrand,
and the added interval has nonnegative integrand.

The prime number theorem applied to (7) therefore yields
\[
\liminf\frac{\log M_K}{K^2}\ge\frac{9s}{16},
\qquad s\ge9.
\]
For \(s=7\), the exceptional rows cost at most
\[
3\int_{1/2}^{(1+\alpha)/2}(1+\alpha-2u)\,du
=\frac{3\alpha^2}{4}.
\]
The \(+1\), floors, and fixed excluded primes affect only lower-order
terms. Thus
\[
\boxed{\liminf\frac{\log M_K}{K^2}
\ge7\left(\frac9{16}-\frac{3\alpha^2}{4}\right)
=\frac{25011}{6400}=3.90796875,\quad s=7.}               \tag{9}
\]

The positive density ratio \(\sigma_{L,s}/w_s\) tends to \(1-1/L\)
at zero and one at infinity, and is bounded away from zero on the
intervening compact sets. Gram comparison, and the factor \(C^h\)
in the leading coefficient, therefore transfer the previously audited
monic Jensen lower bound:
\[
\liminf\frac{\log|F_K(\zeta(s))|}{K^2}>-\frac{1789}{500}.
\]
Here the cited input is the canonical \(K=40n,N=3n\) calculation in
coprime-depletion-primitive-obstruction.md and its independent review.
The transfer costs only \(O_{s,L}(K)\).

Combining with (9) gives a uniform lower margin, for every fixed odd
\(s\ge7\) and fixed \(L\ge2\), of at least
\[
\frac{25011}{6400}-\frac{1789}{500}
=\frac{10559}{32000}>\frac14.
\]
Consequently, for all sufficiently large \(K=40n\),
\[
\boxed{|M_KF_K(\zeta(s))|>\exp(K^2/4).}                  \tag{10}
\]

The monic denominator \(M_K\) is a positive integer and \(M_KF_K\)
is primitive. Every nonzero rational scalar making the original
determinant polynomial integral therefore gives an integer multiple,
up to sign, of this primitive polynomial. Equation (10) proves actual
failure of scalar normalization for this specified family.

The result does not decide \(s=3,5\), changing proportions or numerator
powers, changing the pole grid, or allowing \(L\) to grow with \(K\).
It uses the prime number theorem and the separately certified Jensen
bound explicitly; the coefficient statements (7)–(8) are local.

## 7. Exact reduced-matrix checks

The standard-library checker check_single_depletion_monic.py generates
polynomial moments by cumulants of the uniform distribution on
\(0,\dots,L-1\), and independently checks the initial moments using
exact Bernoulli fractions. It constructs every critical row of the
actual reduced matrix, including its couplings to the complement,
and the full reduced characteristic polynomial from its block factors.

All 31 admissible cases pass, using \(L=2,3,5\), \(s=7,9,11\),
\(K=40,80\), and primes both below and above \(K\). Two additional
requested cases are correctly excluded by the fixed rational
step-value exceptions. Representative exceptional-block results are:

| \(K,N,s,L,p\) | exceptional rows | detected coefficient denominator exponent | proved lower bound |
|---|---:|---:|---:|
| \(80,6,7,2,41\) | 15 | 420 | 315 |
| \(80,6,7,2,43\) | 3 | 413 | 392 |
| \(80,6,7,3,41\) | 15 | 469 | 364 |

Each detected exponent is the exact valuation of a particular coefficient,
not a claim to have located the worst coefficient of the whole polynomial.
The all-prime proof uses the conservative last column and does not assume
the residual nonsingularity seen in these examples.
