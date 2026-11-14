# Apéry's recurrence, the Fano threefold V12, and the Gamma class

Status: primary-source identification, with a proved elementary dictionary
and exact algebraic checks. The geometric theorems are attributed to their
sources and have not been reproved here. Our three-variable Laurent
polynomial and its Beukers relative period are already published mathematics.

## 1. An exact match

The Fano threefold \(X=V_{12}\) is a sevenfold hyperplane section of one
component of the orthogonal Grassmannian \(\mathrm{OG}(5,10)\).
For \(H=c_1(X)\), its anticanonical degree is \(\int_XH^3=12\).
Golyshev–Kerr–Sasaki [GKS], §5.4, display its mirror polynomial

\[
 \phi(X_1,X_2,X_3)=
 \frac{(1-X_1)(1-X_2)(1-X_3)
 [(1-X_1)(1-X_2)-X_1X_2X_3]}{-X_1X_2X_3}.
\]

Changing the signs of all three variables gives **exactly**

\[
 \phi(-x,-y,-z)=\Lambda(x,y,z)
 =\frac{(1+x)(1+y)(1+z)[(1+x)(1+y)+xyz]}{xyz}.       \tag{1}
\]

There is no rescaling of the pencil parameter. Sign changes preserve
constant terms, so the elementary binomial derivation in
[our Laurent note](apery-laurent-reciprocal.md) gives

\[
 A_n=\operatorname{CT}\Lambda^n
 =\operatorname{CT}\phi^n
 =\sum_{k=0}^n\binom nk^2\binom{n+k}k^2.
\]

The period \(\Phi(t)=\sum A_nt^n\) satisfies

\[
 L\Phi=0,\quad
 L=\theta^3-t(2\theta+1)(17\theta^2+17\theta+5)
       +t^2(\theta+1)^3,\qquad\theta=t\,d/dt.       \tag{2}
\]

[GKS], §5.4, identifies the resolved fibers of \(1-t\phi=0\) as K3
surfaces of Picard rank 19, and displays this operator and period sequence.
Thus the recurrence is the third-order period equation for the varying
rank-three part of this K3 family's cohomology. The word “K3” alone does
not determine the coefficients \(17,5\).

## 2. An explicit connection to rational-curve counts

Put \(P(u)=(2u+1)(17u^2+17u+5)\), \(D=z\,d/dz\), and
\(\Psi(z)=\sum A_nz^n/n!\).
Dividing the Apéry recurrence by \(n!\) gives

\[
 (n+1)^4a_{n+1}-P(n)a_n+n^2a_{n-1}=0,\quad a_n=A_n/n!,
\]
\[
 \widetilde L\Psi=0,\qquad
 \widetilde L=D^4-zP(D)+z^2(D+1)^2.                \tag{3}
\]

Golyshev–Zagier [GZ], p.29 and §1.4, give the matrix representative

\[
 M=\begin{pmatrix}
 5&96&1692&12816\\
 1&12&216&1692\\
 0&1&12&96\\
 0&0&1&5
 \end{pmatrix}.                                   \tag{4}
\]

Its entries, subject to the scalar-shift convention below, arise from
genus-zero Gromov–Witten invariants: rational curves with incidence
conditions. The row-vector quantum differential system is

\[
 D(\xi_0,\xi_1,\xi_2,\xi_3)
 =(\xi_0,\xi_1,\xi_2,\xi_3)\,G(z),\qquad
 G_{ij}(z)=M_{ij}z^{j-i+1}.                        \tag{5}
\]

Eliminate the last three components using

\[
\begin{aligned}
 \xi_1&=(D-5z)\xi_0,\\
 \xi_2&=(D-12z)\xi_1-96z^2\xi_0,\\
 \xi_3&=(D-12z)\xi_2-216z^2\xi_1-1692z^3\xi_0,\\
 0&=(D-5z)\xi_3-96z^2\xi_2-1692z^3\xi_1-12816z^4\xi_0 .
\end{aligned}
\]

Substitution with \(Dz=z(D+1)\) gives exactly (3). This is a concrete
route from the enumerative matrix to Apéry's recurrence, once the
Gromov–Witten invariants in (4) are supplied.

**Scalar normalization:** [GZ], §1.2, treats \(M+cI\) as equivalent:
it multiplies \(\Psi\) by \(e^{cz}\) and translates the potential by \(c\).
Matrix (4) uses the Apéry/Eisenstein normalization. The unit-preserving
representative is \(M-5I\): its first column expresses
\(c_1\star1=c_1\). Its solution is \(e^{-5z}\Psi(z)\); the regularized
series is
\(\sum\operatorname{CT}(\Lambda-5)^nt^n=1+48t^2+600t^3+\cdots\).
These two conventions must not be silently mixed.

A particularly short spectral identity is

\[
 \det(uI-M)=u^2(u^2-34u+1).                        \tag{6}
\]

The largest eigenvalue of this nonnegative irreducible matrix is
\(\lambda=17+12\sqrt2\), exactly the positive minimum of \(\Lambda\)
and the exponential growth constant of \(A_n\). For \(M-5I\) the
eigenvalues shift by \(-5\), as they must. This provides a concrete
spectral interpretation of the reciprocal characteristic roots; it
does not compute the enumerative invariants from scratch.

## 3. The term −12t comes from a relative period

The positive integral is

\[
 I_n=\int_{(0,\infty)^3}\Lambda^{-n-1}
                 \frac{dx\,dy\,dz}{xyz}.
\]

Logistic substitutions \(a=x/(1+x)\), etc., give the Beukers cube
integral, with \(I_0=2\zeta(3)\) and \(I_1=10\zeta(3)-12\).
[GKS], §5.4, constructs a higher cycle on the mirror family and pairs
its regulator with the family of holomorphic two-forms. Its resulting
“truncated higher normal function” is explicitly this same integral:

\[
 V(t)=\sum_{n\ge0}I_nt^n,\qquad LV=-12t.             \tag{7}
\]

For the sign convention it is safest to compare their displayed cube
integral with ours, rather than omit the orientation in their
negative-octant integral.

A closed-cycle period satisfies a homogeneous equation; a regulator
or relative-cycle period can have a boundary contribution.
[GKS], Theorem 5.1, shows that for this pencil the inhomogeneous term
has degree at most one and vanishes at zero. Thus \(LV=ct\);
the first two integral values determine

\[
 c=I_1-5I_0=-12.                                  \tag{8}
\]

This is an advanced conceptual alternative to a long integration-by-parts
certificate. Our [divergence proof](geometric-divergence-proof.md)
supplies an elementary exact-form route to the same boundary contribution.

If \(B_0=0,B_1=1\) and \(B_n\) satisfies the recurrence for \(n\ge1\), then

\[
 L\sum B_nt^n=t,\qquad
 V=2\zeta(3)\Phi-12\sum B_nt^n,\qquad
 \widehat V=V/12,\quad \widehat V(0)=\zeta(3)/6.     \tag{9}
\]

The last number is the Apéry constant in these sources. It is neither
\(2\zeta(3)\) nor the cubic Gamma coefficient below.
The literal 12 depends on the normalizations of the cycle, form and
operator: scaling \(V\) by \(1/12\) makes the forcing \(-t\).
Its equality with the degree of \(V_{12}\) does not by itself derive
(8) from the Fano degree.

Why a rational multiple of \(\zeta(3)\)? In this explicit example
[GKS], Example 6.5, puts the limiting cycle in
\(W_{-2}H_M^3(\widetilde X_0,\mathbb Q(3))\). Its bottom-piece
description leads to
\(\mathrm{CH}^3(\operatorname{Spec}\mathbb Q,5)\), whose Borel
regulator is proportional to \(\zeta(3)\). The filtration index is
\(-2\); it must not be confused with the motivic degree or Tate twist 3.
The general existence discussion in §6 invokes a
conjecture, but §5.4 constructs this particular cycle explicitly.
We use the explicit example, not conjectural general existence.

### A localized source correction and a stronger exact dictionary

The downloaded arXiv v3 has a sign error in the displayed birational map
in §5.4. This was checked in the rendered PDF, not inferred from OCR.
Writing \(C=(1-x)(1-y)\), \(W=C-xyz\), the printed map has coordinates
\((z/(1-z),-C/W,x/(1-x))\). At \((-1,-1,-1)\) the original value
is \(\phi=40\), but the transformed value is \(2349/40\), not \(1/40\).
Consequently that printed map does not cover \(t\mapsto1/t\).

Changing all three coordinate signs repairs the algebra:

\[
 \iota(x,y,z)=\left(-\frac z{1-z},\frac C W,
                                      -\frac x{1-x}\right),\qquad
 \iota^2=1,\quad \phi\circ\iota=\phi^{-1}.           \tag{9a}
\]

These are rational identities, on the open set where their denominators
and \(\phi\) are nonzero. Moreover, for \(\Omega=dx\,dy\,dz/(xyz)\),

\[
 \iota^*\Omega=\phi^{-1}\Omega.                    \tag{9b}
\]

Proof: \(1-\iota_1=1/(1-z)\),
\(1-\iota_2=-xyz/W\), \(1-\iota_3=1/(1-x)\).
These identities give (9a) by substitution. The logarithmic Jacobian
has only one possible determinant term involving its first and third
rows, and equals
\(-xyz/[C(1-z)W]=1/\phi\), proving (9b).

Conjugating \(\iota\) by the sign change in (1) gives the corresponding
exact reciprocal map for our polynomial:

\[
 T(x,y,z)=\left(-\frac z{1+z},
       -\frac{(1+x)(1+y)}{(1+x)(1+y)+xyz},
       -\frac x{1+x}\right),\quad
 \Lambda\circ T=\Lambda^{-1},\quad T^*\Omega=\Lambda^{-1}\Omega.
                                                               \tag{9c}
\]

In particular,
\(T^*(\Lambda^n\Omega)=\Lambda^{-n-1}\Omega\).
This is a pointwise rational explanation for the exponent \(-n-1\)
in our inverse-power integral. It does not identify its noncompact
integration chain with the torus defining a constant term.
The two chains still have different boundaries and different equations.

We record the typo only for this retrieved source version; no claim is
made about other versions. The mirror polynomial, displayed cube
integral and its differential equation match our independently established
identities. The regulator interpretation is a cited geometric result,
not a proof newly certified by repairing this one rational map.

## 4. Chern numbers give the exact Gamma coefficient

For Chern roots \(\delta_i\), set

\[
 \widehat\Gamma_X=\prod_i\Gamma(1+\delta_i),\qquad
 \log\Gamma(1+u)=-\gamma u+
       \sum_{k\ge2}\frac{(-1)^k\zeta(k)}k u^k .
\]

Complex dimension three truncates this calculation at degree three.
[GZ], §1.5, supplies

\[
 \int H^3=12,\qquad \int Hc_2=24,\qquad
 \int c_3=\chi(X)=4-2h^{1,2}=-10,\quad h^{1,2}=7 .
\]

Even rational cohomology is generated by \(H\), so
\(c_2=2H^2,\ c_3=-5H^3/6\).
Newton's identities give
\(\sum\delta_i^2=-3H^2\) and \(\sum\delta_i^3=-15H^3/2\).
Exponentiating the displayed logarithm therefore proves

\[
 \widehat\Gamma_X
 =1-\gamma H+\frac{\gamma^2-3\zeta(2)}2H^2+
 \left(-\frac{\gamma^3}6+\frac32\gamma\zeta(2)
                  +\frac52\zeta(3)\right)H^3.      \tag{10}
\]

This calculation is elementary once the three Chern numbers are known.
The deep theorem, [GZ], Theorem 1/3, identifies these coefficients with
the limits \(\kappa_j=\lim_{z\to+\infty}\Psi_j(z)/\Psi_0(z)\) for the
normalized Frobenius basis of (3).

The modified class is

\[
 \widehat\Gamma'_X=\frac{\widehat\Gamma_X}{\Gamma(1+H)}
 =1-\zeta(2)c_2+\zeta(3)(Hc_2-c_3)
 =1-2\zeta(2)H^2+\frac{17}{6}\zeta(3)H^3.          \tag{11}
\]

Here the numerator 17 is topological: \(24-(-10)=34=2\cdot17\).
If \(S\) is an anticanonical K3 section, then
\(\chi(X\setminus S)=-10-24=-34\). Thus the cubic coefficient is
\(-\chi(X\setminus S)\zeta(3)/12\).
Its relation to the 17 in the recurrence uses the full dictionary;
an Euler characteristic alone does not determine a recurrence.

The Gamma factor also has a direct analytic explanation. Deform the
Apéry recurrence by replacing \(n\) with \(n+\varepsilon\), with
\(A_{-1}(\varepsilon)=0,A_0(\varepsilon)=1\).
The coefficients of the corresponding Borel-side deformation are
\(a_n(\varepsilon)=A_n(\varepsilon)/(1+\varepsilon)_n\).
For the Frobenius series including \(t^\varepsilon,z^\varepsilon\),
coefficientwise Laplace integration gives

\[
 \int_0^\infty e^{-u}\Psi(tu,\varepsilon)\,du
     =\Gamma(1+\varepsilon)\Phi(t,\varepsilon),     \tag{12}
\]

since the integral of \(e^{-u}u^{n+\varepsilon}\) is
\(\Gamma(n+1+\varepsilon)\).
This explains why the modified Gamma class divides by \(\Gamma(1+H)\).
[GZ], §2.1–2.2, proves coefficientwise through cubic order that

\[
 \lambda^{-\varepsilon}
 \lim_{n\to\infty}\frac{A_n(\varepsilon)}{A_n}
 =1-2\zeta(2)\varepsilon^2+
             \frac{17}{6}\zeta(3)\varepsilon^3+O(\varepsilon^4).
                                                               \tag{13}
\]

Multiplication by \(\Gamma(1+\varepsilon)\) gives (10).
This deformed limit must not be confused with \(B_n/A_n\to\zeta(3)/6\).

## 5. Scope and the higher-odd-zeta question

The precise explanations are distinct:

* The recurrence is both a K3 period equation and the regularization of a
  quantum differential equation encoding rational-curve counts.
* The forcing term is a relative-cycle/regulator boundary contribution;
  its shape follows from geometric singularity constraints.
* Odd zeta values appear in \(\log\Gamma\), with coefficients determined
  here by Chern numbers and a proved asymptotic theorem.

The relation to the Todd/index class comes from
\(\Gamma(1+u)\Gamma(1-u)=\pi u/\sin(\pi u)\); see [GZ], §1.5.
This product cancels the odd terms, so the Todd class alone does not
recover the odd-zeta information.

None of these geometric results supplies by itself the integer-denominator
bound or the nonzero small integer required by irrationality.
The Gamma-conjecture proof in [GZ] uses hypergeometric or modular methods;
citing it is not a shorter self-contained elementary irrationality proof.

For \(\zeta(5)\) the threefold Gamma class truncates too early. One needs
higher-dimensional geometry, and products such as \(\zeta(2)\zeta(3)\)
then occur as well. A fifth-weight regulator or Gamma-class appearance
does not isolate a single \(\zeta(5)\) linear form with controlled
denominators. That arithmetic isolation and decay remain a separate bridge.

## Verification and primary references

[Exact checker](../computations/check_fano_apery_dictionary.py) passes:
the signed mirror identity, the corrected reciprocal involution and
logarithmic Jacobian (and the printed-map counterexample),
constant terms through \(n=6\), elimination
of the 4-by-4 system, its characteristic polynomial, the Chern/Gamma
algebra, and 30 recurrence normalizations. These verify the elementary
dictionary, not the imported geometric theorems.

**[GZ]** V. Golyshev and D. Zagier, *Proof of the gamma conjecture for
Fano 3-folds of Picard rank 1*, Izvestiya Mathematics 80:1 (2016), 24–49,
DOI 10.1070/IM8343.
[Primary author-hosted PDF](https://people.mpim-bonn.mpg.de/zagier/files/tex/GammaConjecture/IZV_80_1_24.pdf).
Read pp.26–29, §§1.4–1.5 and §§2.1–2.2.

**[GKS]** V. Golyshev, M. Kerr and T. Sasaki, *Apéry extensions*,
arXiv:2009.14762v3 (10 July 2023), published DOI 10.1112/jlms.12825.
[Primary arXiv PDF](https://arxiv.org/pdf/2009.14762).
Read Theorem 5.1, §5.4/formulas (5.4)–(5.5), and §6/Example 6.5.

PDFs and extracted texts are saved under codex/sources.
The [source manifest](../sources/fano-apery-source-manifest.json)
records URLs, retrieval times, byte sizes and SHA-256 hashes.
The additional Galkin–Golyshev–Iritani paper there is background;
the claims in this note are sourced to [GZ] and [GKS].
