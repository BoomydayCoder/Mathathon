# R20 — "Infinitely many tensor-stable SL_2(p) Weil lattices?": what can be proved

Round 20 (single agent, ~90 min, 2026-11-25 18:25–19:55 UTC).  Labels: **PROVED** (complete argument written here),
**PROVED-MODULO-CITATION** (cited fact named), **CONDITIONAL**, **NUMERICAL**, **HEURISTIC**, **CONJECTURE**.
Nothing here was checked by a second agent or a human.

## 0. Verdict

* The Target (stability (BS) of `W_p` for infinitely many `p = 1 mod 4` with `(p+1)/2` prime) is **NOT proved and NOT
  disproved**.
* New structural result (Theorem 1, PROVED-MODULO-CITATION + NUMERICAL confirmation): for `p = 5 mod 8` the Weil
  lattice is an explicit **ideal lattice of the cyclotomic field `Q(zeta_{4p})`** ("two-prime Craig lattice"):

      W_p  ~  ( P1^{(m+1)/2} P2^{-(m-1)/2} ,  Tr_{Q(zeta_4p)/Q}(x xbar) / (2p) ),      m = (p-1)/4,

  `P1, P2` the two primes of `Q(zeta_4p)` over `p`.  `E_8` (p=5) and the Leech lattice (p=13, with the exact
  `s`-distribution 78624/45864/13104 of Round 5) come out of this formula in exact arithmetic, and so do all
  previously computed torus-fixed sublattices (`D_4^6` with `|s|=3` at p=37, `E_8^3` at p=61, the extremal
  40-dimensional lattice at p=61, `D_4^6` at p=157, Leech at p=277).  The torus-fixed lattices are the same
  construction in the subfields (Theorem 2), which re-proves Lemma B, the "sqrt(m) x unimodular" observation of
  Round 13 and the self-similarity of Round 16 arithmetically, with no curve and no group orbit.
* Consequences: (A) is **REDUCED-TO** a sharply stated inequality about one ideal of `Q(zeta_4p)` (Theorem 3), an
  unconditional bound `N(nu(w)) >= r^2 p^{-1+1/r}` (Theorem 4, off by a factor ~p), an exact and fast algorithm for
  every structured sublattice (Section 6; new exhaustive scans, no violator found for any `p >= 61`), and a proof
  that (B), (C) can **not** be deduced from (A) by Hermite/packing inequalities for `p > 73` (Section 5).
* **Partial result on (A) (Theorem 5, PROVED given Theorem 1):** for every `p = 5 mod 8` the vector condition holds
  on the explicit finite-index sublattice `Z[zeta_{4p}]` (= the Borel sublattice `Z[B]v_0`, index `p^{(p-3)/2}`), with
  equality exactly on the roots of unity (times units of `F`); every violator `x` must have `|N(x)| < 1` in the
  model `L_0`.  In Mordell-Weil language: every map in the Z-span of the Borel translates (and Frobenius twists) of
  the Galois map `C_p -> D -> E` has `N(q_F) >= r^2`, in particular degree `>= (p-1)/2`.
* Elementary restatement (Corollary 4.2, checked at p = 5, 13): `Hom(Jac C_p, E)` is the lattice of pairs of integer
  sequences on `F_p` that are mod `p` value tables of polynomials `f, f'` with `deg <= (3m-1)/2`,
  `deg(f - iota f') <= (m-1)/2`; `deg = V/(2p)` (`V` = total centred square sum), and (A) reads
  `V^2 - S^2/p >= p^2(p-1)^2` with `S` the Paley quadratic form.
* The ideal class of `W_p` in `Cl(Q(zeta_{4p}))` is the class of `P^{(p-1)/4}` (`P | p`), which has **order 1 or 2**
  (Stickelberger) — Remark 4.1.
* For `p = 1 mod 8` (73, 193, ...) no model and no computation was obtained (Section 7 says exactly why).
* Novelty was not checked.  G. Nebe, "Some cyclo-quaternionic lattices", J. Algebra 199 (1998), describes the faithful
  `SL_2(p)`-lattices of rank `2(p-1)` and very likely contains an equivalent description of `W_p`; neither this paper
  nor Gross (1990), Tiep (1997), Dummigan (1997) could be read in this session (two web searches returned only
  bibliographic data).  It is consistent with the known fact (quoted from memory, Bayer-Fluckiger) that the Leech
  lattice is an ideal lattice over `Q(zeta_52)`.  What is certainly specific to this project is the use of the model
  for the stability question (Theorems 3-5, Section 5, the scans).

Notation: `p = 1 mod 4` prime, `G = SL_2(p)`, `B = U x| T` the Borel subgroup, `u = [[1,1],[0,1]]`,
`d_t = diag(t, 1/t)`, so `d_t u d_t^{-1} = u^{t^2}`; `g` a primitive root mod p, `d = d_g`; `zeta = zeta_p`,
`K = Q(zeta)`, `F = Q(sqrt p) subset K`, `sigma = sigma_{g^2}` (generator of `Gal(K/F)`, order `(p-1)/2 = 2r`),
`r = (p-1)/4`; `W` the faithful irreducible `Q[G]`-module with character `2(psi+psi')`, `A = M_r(Delta)` the image
of `Q[G]` in `End_Q(W)`, `Delta = End_G(W)`.  For a vector `w`: `m_+-(w)` = squared norms of the `eta = +-sqrt p`
components, `nu(w)` the Hermitian norm, `N(nu(w)) = (p/4) m_+ m_-` (unimodular normalisation), threshold `r^2`.

---------------------------------------------------------------------------------------------------------------

## 1. The Borel subgroup already generates everything

**Lemma 1.1 (PROVED-MODULO-CITATION: character table of SL_2(p)).**  `W|_U ~ K^2` (U has no fixed vector), so `W`
is a 2-dimensional `K`-vector space via `zeta.w := u w`, and `d` acts `sigma`-semilinearly with
`d^{2r} = d_{-1} = -1`.  The image of `Q[B]` in `End_Q(W)` is all of `A`.  Hence `End_B(W) = Delta`, every
`B`-invariant bilinear form on `W` is `G`-invariant, and every `B`-linear map between `G`-lattices in `W` is `G`-linear.

*Proof.*  `psi(u^j)` are the two Gauss periods, so `psi|_U + psi'|_U` is the sum of all nontrivial characters and
`W|_U = K^2`.  The relations `d zeta d^{-1} = sigma(zeta)`, `d^{2r} = -1` make `W` a module over the cyclic algebra
`C = (K/F, sigma, -1) = sum_{j<2r} K d^j`, which is central simple over `F` of dimension `(2r)^2 = 4r^2`, and the
image of `Q[B]` in `End_Q(W)` (generated by `u` and `d`) is a homomorphic image of `C`.  A simple algebra acts
faithfully, so the image of `Q[B]` has `F`-dimension `4r^2 = dim_F M_r(Delta) = dim_F A`, and it is contained in `A`;
hence it is `A`.  The rest is formal (for forms: a `B`-invariant form is
`b_0(., c.)` with `b_0` a fixed positive `G`-invariant form and `c` commuting with `B`, i.e. `c in Delta`).  QED

**Lemma 1.2 (PROVED).**  Let `p = 5 mod 8` (r odd).  Put `L = K(i) = Q(zeta_{4p})`, let `sigma~` be the extension of
`sigma` to `L` with `sigma~(i) = i`, and let `B` act on `L` by `u.x = zeta x`, `d.x = i sigma~(x)`.  Then `W ~ L` as
`Q[B]`-modules (hence as `A`-modules), multiplication by `i` is an element `xi in Delta` with `xi^2 = -1`, and the
`G`-invariant forms are exactly `Tr_{L/Q}(alpha x ybar)`, `alpha in F`.

*Proof.*  `d zeta d^{-1} = sigma(zeta)` and `d^{2r} x = i^{2r} sigma~^{2r} x = (-1)^r x = -x` because `r` is odd.  So
`L` is a `C`-module of `F`-dimension `4r`, which is the dimension of the simple `C`-module (`C = M_r(Delta)`), and `W`
is the simple module.  Multiplication by `i` commutes with `zeta` and `d`.  The forms `Tr(alpha x ybar)`, `alpha in F`,
are `B`-invariant (`sigma~` commutes with complex conjugation and fixes `F`), they form a 2-dimensional space, and the
space of `G`-invariant forms is 2-dimensional (`= F`).  QED   (For `r` even `d^{2r} = +1` in this model: Section 7.)

## 2. Theorem 1: the explicit model (p = 5 mod 8)

Let `m = r = (p-1)/4` (odd), `iota in Z_p` with `iota^2 = -1`, `P1 = (1-zeta, i-iota)`, `P2 = (1-zeta, i+iota)` the
primes of `L` over `p` (`P1 P2 = (1-zeta)`, `P_j^{p-1} = (p)`, complex conjugation swaps them, `sigma~` fixes them).

**Theorem 1 (PROVED-MODULO-CITATION; NUMERICAL confirmation p = 5, 13 and sublattices for 29 primes).**
(a) Let `Lambda subset W = L` be a `G`-stable lattice which is stable under `xi` and satisfies
`Lambda_l = Z_l[zeta_{4p}]` for all primes `l != p` (in a suitable `B`-equivariant identification).  Then

      Lambda = P_a^{(1-3m)/2 + 2m j} P_b^{(1-m)/2 + 2m j'}     for some j, j' in Z,   {a,b} = {1,2},

and conversely every such ideal is a `G`-stable lattice.  (Which of `P1,P2` is `P_a` depends on `g` and on the sign of
`iota`; complex conjugation exchanges the two cases and is an isometry.)
(b) Each of them carries a unimodular `G`-invariant form `Tr(alpha x xbar)`, `alpha in F`, and every unimodular
`G`-lattice in `W` is even.  Up to the isometries `x -> sqrt(p) x` and `x -> (g(chi)/g(chibar)) x` (`chi` a quartic
character, `g` the Gauss sum; this element has absolute value 1 in every embedding and divisor `(P1/P2)^{+-2m}`) there
are two such lattices:

      L_0 = ( P1^{-(3m-1)/2} P2^{-(m-1)/2} , Tr(eps x xbar /(2 sqrt p)) ),   L_1 = ( P1^{(m+1)/2} P2^{-(m-1)/2} , Tr(x xbar)/(2p) ),

`L_1 = L_0 . delta` with `delta = sqrt p` on the `P1`-component and `1` on the `P2`-component (`delta in Delta_p^x`).
(c) The lattice `Z[G] v_0` used in Rounds 10-13 (with its unimodular member at `t = +-log eps`) is `L_0`; all
`G`-lattices of the form `Lambda I` (`I` an ideal of the maximal order) have the same values `R_U`, so **(BS) may be
tested on `L_1`**.
(d) In `L_1`:  `|x|^2 = (T_+ + T_-)/(2p)`, `m_+- = T_+-/(2p)`, `N(nu(x)) = T_+(x) T_-(x)/(16 p)`, where
`T_+-(x) = sum |tau(x)|^2` over the `p-1` embeddings `tau` of `L` with `tau(sqrt p) = +-sqrt p`.

*Proof.*  (i) *Primes l != p.*  The image of `Z_l[B]` in `A_l` is the crossed product order
`(Z_l[zeta]/Z_l[sqrt p], sigma, -1)`.  `Z_l[zeta]` is unramified (Galois as an extension of rings) over
`Z_l[sqrt p]`, so this order is Azumaya, hence maximal [citation: crossed products of Galois ring extensions are
Azumaya; Auslander-Goldman].  The image of `Z_l[G]` is an order containing it, so they are equal: at `l != p` *every*
`B`-lattice is a `G`-lattice.  In particular `Z_l[zeta_{4p}]` is `G`-stable.
(ii) *The prime p.*  `Z_p[i] = Z_p e_+ + Z_p e_-`, `e_+- = (1 +- i/iota)/2 in Z_p[xi]`.  So `Lambda_p = e_+ Lambda_p +
e_- Lambda_p`, both summands `G`-stable, each a fractional ideal `p^c` of the DVR `Z_p[zeta]` (`p = (1-zeta)`) inside
`e_+- L_p = Q_p(zeta)`, on which `d` acts by `x -> +-iota sigma(x)`.  As a `Q_p(sqrt p)[G]`-module `e_+- W_p` has
character `psi`.  The Brauer character of `psi` equals that of the simple module `Sym^{2m-1}(F_p^2)`: on the split
torus `sum_{j=0}^{2m-1} a^{2m-1-2j} = a^{-(2m-1)} (a^{4m}-1)/(a^2-1) = 0` for `a^2 != 1`, on the non-split torus the same
sum is `-b^{-(p+1)/2} = +-1` with the sign of the character table [citation: character table of SL_2(p); simple
modules `Sym^k`, `k < p`].  Hence the reduction mod `sqrt p` of any `G`-lattice `M` in `e_+- W_p` is simple, and the
`G`-lattices in `e_+- W_p` are exactly `sqrt(p)^j M` (every `G`-lattice is a `Z_p[sqrt p]`-module because
`eta in Z[U]`; if `M'` is another one, scale it so that `M' subset M`, `M' not subset sqrt(p) M`; the image of `M'` in
`M/sqrt p M` is a nonzero submodule, so `M' + sqrt p M = M` and `M' = M` by Nakayama).  Since `sqrt p Z_p[zeta] = p^{2m}`, the admissible exponents `c` form one class mod `2m`.
(iii) *Which class.*  `M/pM` (`p = (1-zeta)`) is the space of `U`-coinvariants of `M/sqrt p M = Sym^{2m-1}`.  With
`u e_1 = e_1`, `u e_2 = e_1+e_2` the image of `u-1` on `Sym^{2m-1}` is spanned by the monomials divisible by `e_1`, so the
coinvariants are spanned by `e_2^{2m-1}`, on which `d_t` acts by `t^{-(2m-1)}`.  On the other hand, for `M = p^c` in
the component where `d = iota' sigma` (`iota' = +-iota`): `d (1-zeta)^c = iota' (1-zeta^{g^2})^c = iota' g^{2c} (1-zeta)^c
mod p^{c+1}`.  So `iota' g^{2c} = g^{-(2m-1)} mod p`.  As `iota' = g^{m}` or `g^{-m} mod p` (the two components get the
two values), `c = (1-3m)/2` resp. `(1-m)/2 mod 2m`.  This gives (a); the converse holds because by (ii) `G`-stability of
`p^c` depends only on `c mod 2m` and at least one `G`- and `xi`-stable lattice exists (take any `G`-lattice and
multiply by a maximal order containing `xi`), and by (i) the components at `l != p` are irrelevant.
(iv) *Forms.*  `(A, Tr(alpha x xbar))` is unimodular iff `A Abar alpha D = (1)`, `D = 2 (P1P2)^{4m-1}` the different.
For `A = P_a^{c_a} P_b^{c_b}` as in (a) this reads `(alpha) = (1/2) (sqrt p)^{-(1+j+j')}`, solvable by a totally
positive `alpha in F` because `N(eps) = -1`.  For `j=j'=0`: `alpha = eps/(2 sqrt p)`; for `(j,j')=(1,0)`:
`alpha = 1/(2p)`.  Evenness: a unimodular lattice `M` has a characteristic class `w_0 + 2M` in `M/2M`, invariant under
all isometries; `u` has no nonzero fixed vector on `M/2M = (F_2[X]/Phi_p)^2` because `Phi_p(1) = p` is odd; so
`w_0 in 2M` and `M` is even.  Multiplication by `sqrt p` and by `rho = g(chi)/g(chibar)` commutes with `zeta`, `i` and `d`
(`sigma~ g(chi) = chibar(g^2) g(chi) = -g(chi)`), is `B`-linear hence `G`-linear (Lemma 1.1), and multiplies the form
by a constant; by Stickelberger `(g(chi)) = P_a^m P_b^{3m}`, so these maps shift `(j,j')` by `(1,1)` and `(1,-1)`: two
classes remain, represented by `L_0` and `L_1`.
(v) *(c).*  For `Lambda = Z[G]v_0` (`v_0` real and fixed by `C_m`): it is `xi`-stable (in the Weil model
`xi = rho(d_i)^{-1} o (complex conjugation of coordinates)`), `Z_l[G]v_0 = Z_l[B]v_0 = Z_l[zeta_{4p}] v_0` for `l != p`
by (i), and at `p` it is the smallest `G`-lattice containing `e_+- v_0`, i.e. `c_+- in (-2m, 0]`: exponents
`(-(3m-1)/2, -(m-1)/2)`.  [Check: the index of `Z[B]v_0` is then `p^{2m-1} = p^{(p-3)/2}`, which reproduces
`det_B(Lambda) = (r^2 p)^{p-1}`, the value computed in Round 10 for p = 5,...,61.]  Invariance of `R_U`: for a
`Delta`-subspace `U` and `delta in Delta_p^x`, `(Lambda delta) ∩ U = (Lambda ∩ U) delta` and `Lambda ∩ U` is locally free
over the maximal order, so numerator and denominator of `R_U` scale by the same power of `N(nrd delta)`.
(vi) *(d).*  In `L_1` the form is `Tr(x xbar)/(2p)`, the `eta = +-sqrt p` components are the embeddings over the two
places of `F`, so `m_+- = T_+-/(2p)` and `N(nu) = (p/4) m_+ m_-`.  QED

*Caveats.*  The proof uses the cited facts only through (i)-(iii).  Step (v) uses that the element called `xi` in the
Weil model preserves `Z[G]v_0`; this was not re-derived in detail here (it is not needed for (a), (b), (d) and for the
use of `L_1`).  I did not prove that `L_0` and `L_1` are isometric (for p = 5, 13 they are: `E_8`, Leech).  NUMERICAL (`t_L0.py`):
for (p,N) = (5,4), (13,12), (37,12), (61,12), (61,20) the `L_0`-type and `L_1`-type lattices are both even
unimodular and have identical tables of (norm, |s|) for all vectors of norm <= 4; `full61.py`: the full lattices
from the formula for p = 29, 37, 61 (ranks 56, 72, 120; both types) have integral even Gram matrix `G` with
`det G = 1` and `det E = p^{n/2}` exactly.

**Coordinates (PROVED from Theorem 1).**  Write `x = z/(1-zeta)^{(m-1)/2}`, `z = y_1 + i y_2`, `y_j in Z[zeta]`.  Then

      L_1 = { (y_1,y_2) in Z[zeta]^2 :  y_1 + iota y_2 in (1-zeta)^m },   |w|^2 = Tr_{K/Q}( (|y_1|^2+|y_2|^2) / (p |1-zeta|^{m-1}) ).

Scalar-form version (multiply by `sqrt p`, i.e. exponents `((5m+1)/2,(3m+1)/2)`): `sqrt(p) L_1` is the lattice of
pairs `(c,c')` of integer vectors indexed by `F_p`, each of sum 0, with Euclidean norm `|c|^2+|c'|^2`, such that mod `p`

      c = (f(j))_j,  c' = (f'(j))_j,   deg f, deg f' < (5m+1)/2,   deg(f + iota f') < (3m+1)/2      (Reed-Solomon type code),

and `eta` acts through the Paley matrix `Q = ((j-l)/p)` on each block (`m_+- = (|c|^2 +- c^T Q c/sqrt p)/2p` summed over
the two blocks).  `L_1` contains `A^{(a)} + i A^{(a)}` with index `p^m`, `A^{(a)}` the Craig lattice `(1-zeta)^a Z[zeta]`,
`a = (m+1)/2`, whose minimum is `>= 2a = (p+3)/4` (Craig; Conway-Sloane Ch. 8, cited).

## 3. Theorem 2: the tower of quotient lattices

**Theorem 2 (PROVED given Theorem 1; NUMERICAL confirmation).**  Let `p = 5 mod 8`, `k | m` (so `k` odd), `N = (p-1)/k
= 4m'`, `M_N` the subfield of `K` of degree `N`, `P_{1,N}, P_{2,N}` the primes of `M_N(i)` over `p`.  Then the fixed
lattice of the split-torus subgroup `C_k` in `L_1` is

      Fix(C_k) = L_1 ∩ M_N(i) = sqrt(k) x L_N(p),    L_N(p) := ( P_{1,N}^{(m'+1)/2} P_{2,N}^{-(m'-1)/2} , Tr_{M_N(i)/Q}(x xbar)/(2p) ),

an even unimodular lattice of rank `2N` of exactly the same shape with `m` replaced by `m' = N/4`.  In particular
`R_{Fix(C_k)} = 1` (Lemma B), and the condition `N(nu) >= r^2` restricted to `Fix(C_k)` is the condition
`N(nu') >= m'^2` on `L_N(p)` (self-similarity).  `L_4(p) = E_8` for every `p`; `L_12(p)` is the lattice called
`L_24(p)` in Round 13.

*Proof.*  `C_k = <d^N>` and `d^N x = i^N sigma~^N x = sigma~^N x` (`4 | N`), `sigma^N` has order `k`, so the fixed space
is `M_N(i)`.  For `x in M_N(i)`, `v_{P_j}(x) = k v_{P_{j,N}}(x)`, so `L_1 ∩ M_N(i) = P_{1,N}^{ceil((m+1)/2k)}
P_{2,N}^{ceil(-(m-1)/2k)}` and `ceil((m+1)/(2k)) = (m'+1)/2`, `ceil(-(m-1)/(2k)) = -(m'-1)/2` because `m = k m'`, `m'`
odd, `0 < 1/(2k) <= 1/2`.  The form restricts to `k Tr_{M_N(i)/Q}(x xbar)/(2p)`.  Unimodularity: `N(ideal)^2 |disc| /
(2p)^{2N} = p^2 4^N p^{2N-2}/(2p)^{2N} = 1`.  `det Fix(C_k) = k^{2N}` and `n/(2N) = k`, so `R_{Fix(C_k)} = k^{2N}/k^{2N} = 1`.  QED

Explicit coordinates (used by the code): `x = z/lam^{(m'-1)/2}`, `lam = N_{K/M_N}(1-zeta)`, `z = y_1 + i y_2`,
`y = sum_j a_j theta_j` in the Gaussian periods `theta_j = sum_{h in H} zeta^{g^j h}`; `z in P_{1,N}^{m'}` iff
`sum_j (a^{(1)}_j + iota a^{(2)}_j) gamma^{jl} = 0 mod p` for `0 <= l < m'`, `gamma = g^k`.

## 4. What Theorem 1 gives for (A)

**Theorem 3 (REDUCTION of (A); PROVED given Theorem 1).**  For `p = 5 mod 8` the vector condition (A) for the
member `L_1` of the family (all `w != 0` in `L_1` have `N(nu(w)) >= r^2`) is equivalent to

      (A')   T_+(x) T_-(x) >= p (p-1)^2    for every 0 != x in A_p := P1^{(m+1)/2} P2^{-(m-1)/2} subset Q(zeta_{4p}),

equivalently: for every totally positive `alpha in F (x) R` of norm 1 the ideal lattice `(A_p, Tr(alpha x xbar)/(2p))`
has minimum `>= (p-1)/sqrt p`.  Equivalently, with `x = x_1 + i x_2` and `x_+ := |x_1|^2 + |x_2|^2 in K^+`:
`N_{F/Q}( Tr_{K^+/F}(x_+) ) >= p r^2`.
(The vector condition for the member `L_0` — the lattice of Rounds 10-13 — is (A'') of Theorem 5.  (A') and (A'') are
both necessary for (BS), which is a property of the whole family by Theorem 1(c); I did not prove that (A') and (A'')
are equivalent to each other, because the principal lines of `L_0` and `L_1` need not be the same lines.)

**Theorem 4 (unconditional bound; PROVED given Theorem 1).**  For all `0 != w` in `L_1` (and in `L_0`, by Theorem 5(b)), `p = 5 mod 8`:

      N(nu(w)) >= r^2 p^{-1+1/r}        (needed: r^2).

*Proof.*  `x xbar in A_p conj(A_p) = P1P2 = (1-zeta)`, and `x_+ = (1/2)Tr_{L/K}(x xbar)` is the `K`-component of
`x xbar in (1-zeta)Z[zeta] + (1-zeta)Z[zeta] i`; it is real, so `x_+ in p_+ := (1-zeta)Z[zeta] ∩ K^+` and
`N_{K^+/Q}(x_+) >= p` (`x_+` is totally positive, nonzero).  `T_+- = 4 Tr_{K^+/F}(x_+)` at the two places of `F`, and by
AM-GM over the `r` conjugates, `T_+T_- >= 16 r^2 N_{K^+/Q}(x_+)^{1/r} >= (p-1)^2 p^{1/r}`; divide by `16p`.  QED

Remarks.  (1) Sanity: at p = 13 this gives `N >= 9 * 13^{-2/3} = 1.63` (truth 9, borderline); at p = 37:
`81 * 37^{-8/9} = 3.27` (truth <= 63 < 81).  The bound does not prove stability anywhere, as it must not at p = 37.
(2) It improves the bound `((p+3)/(4 b sqrt p))^2` of Round 11 (independent of the fundamental unit) but is still a
factor `p^{1-1/r}` short.  (3) AM-GM would give (A) from the *stronger* statement `N_{K^+/Q}(x_+) >= p^r`; that
statement is false at p = 13 (for Leech minimal vectors with `s = +-4`, `x_+ notin F`, so AM-GM is strict while
`N(nu) = r^2`) — HEURISTIC: for a generic rank-2 hermitian `Z[zeta_p]`-lattice it would hold for `p^{1/4} > pi`,
i.e. `p > 97`.  (4) The square-root-bound trick of quadratic residue codes (`x . conj(tau x) in (1-zeta)` for
`tau = sigma_n`, `n` a non-residue, and `|Tr(x conj(tau x))| <= 2 sqrt(T_+T_-)`) gives only `T_+T_- >= p^2/4` when the
trace is nonzero: no gain.  (5) Craig's argument bounds the Craig sublattice only, and the Hermite/unit-orbit
arguments use only the minimum, which can never suffice (Round 5, Lemma 3).  (6) Siegel-Weil over the genus of the
hermitian lattice shows that the mass-proportion of lattices with a violating vector is `<~ (pi e/sqrt p)^{p-1}`,
but the mass of the genus is super-exponentially large, so nothing follows for the single lattice `W_p`.

**Theorem 5 ((A) on an explicit sublattice of finite index; PROVED given Theorem 1).**  Work in the model
`L_0 = (A_0, Tr(eps x xbar/(2 sqrt p)))`, `A_0 = P1^{-(3m-1)/2} P2^{-(m-1)/2}`, which contains `O_L = Z[zeta_{4p}]` with
index `p^{2m-1} = p^{(p-3)/2}`.  There `N(nu(x)) = T_+(x)T_-(x)/16 = N_{F/Q}( (1/4) Tr_{L/F}(x xbar) )`, so
(a) the vector condition (A) for `L_0` is equivalent to

      (A'')   N_{F/Q}( Tr_{L/F}(x xbar) ) >= (p-1)^2 = N_{F/Q}( Tr_{L/F}(1) )     for all 0 != x in A_0

("1 minimises the F-norm of the relative trace form on the over-ideal `A_0` of `O_L`").
(b) For every `x != 0`: `N(nu(x)) >= r^2 p^{v(x)/r}`, where `v(x) := v_{p_+}(|x_1|^2+|x_2|^2) >= -(r-1)`.  In
particular `N(nu(x)) >= r^2` whenever `v_{P1}(x) + v_{P2}(x) >= -1` (then `x_+` is `p`-integral): for all `x in O_L`,
for all `x in P1^{-j}P2^{j} ∩ A_0` and `x in P1^{-j-1}P2^{j} ∩ A_0`, `j in Z`.  In the language of Corollary 4.2:
`N(q_F) >= r^2 p^{-floor((d_1+d_2)/2)/r}` with `d_1 = deg(f + iota f')`, `d_2 = deg(f - iota f')` (constants count as 0),
so (A) holds for all classes with `d_1 + d_2 <= 1`.
(c) On `O_L` equality `N(nu(x)) = r^2` holds exactly for `x in mu_{4p} . eps^Z` (the `U`- and `xi`-translates of the
torus-fixed borderline vector, `p` lines).
(d) Every violator of (A) satisfies `|N_{L/Q}(x)| < 1` and `v(x) < 0`.

*Proof.*  `m_+ = (eps/(2 sqrt p)) T_+`, `m_- = (1/(2 eps sqrt p)) T_-` (as `-eps' = 1/eps`), `N(nu) = (p/4) m_+m_-`.
(b) as in Theorem 4: `x_+` is totally positive in `p_+^{v}`, so `N_{K^+/Q}(x_+) >= p^v`, and AM-GM over the `r`
conjugates over `F`; `x xbar in A_0 conj(A_0) = (1-zeta)^{-(2m-1)}` gives `v >= -(m-1)`.  (c) Equality forces
`x_+ in F` and `N(x_+) = 1`, so `x_+ = eps^{2j}`; then `x' = eps^{-j} x` has `|x'_1|^2 + |x'_2|^2 = 1` in every
embedding, so `Tr_{K/Q}|x'_1|^2 + Tr_{K/Q}|x'_2|^2 = p-1`, while a nonzero `y in Z[zeta]` has `Tr|y|^2 >= p-1` with
equality iff `y` is a root of unity (AM-GM, Kronecker); hence one of `x'_1, x'_2` vanishes and the other is `+-zeta^a`.
(d) AM-GM over all `2(p-1)` embeddings: `T_+T_- >= (p-1)^2 |N(x)|^{2/(p-1)}`.  QED

Under the dictionary of Round 11 (M3) (unimodular point and Mordell-Weil point antipodal; PROVED there modulo
literature) the Mordell-Weil lattice is `(A_0, deg = Tr_{L/Q}(x xbar)/4)`, `q_F = Tr_{L/F}(x xbar)/4`; `x = 1` is the
Galois map `phi_0: C_p -> D -> E` of degree `(p-1)/2` (check: p = 5: degree 2; p = 13: degree 6).  So "(A) at the
Mordell-Weil point" says that the roots of unity stay minimal vectors of the trace form when `O_L` is enlarged to
`A_0`, and Theorem 5 proves `deg >= (p-1)/2` and `N(q_F) >= r^2` for all `phi` in `Z[B]phi_0 + Z[B] xi phi_0`.
(The identification of `O_L` with `Z[B]v_0` uses Theorem 1(v); inside the model the statement is unconditional.)
Sanity: the p = 37 violators (N = 63 < 81) must have `v <= -1`: no contradiction; p = 13: no violator exists, fine.

**Remark 4.1 (PROVED-MODULO-CITATION: Stickelberger).**  `A_p = P1^m (1-zeta)^{-(m-1)/2}`, so `[A_p] = [P1]^m` in
`Cl(L)`.  For a quartic character `chi`, `(g(chi)) = P1^m P2^{3m}` (or with `P1, P2` exchanged), hence
`P2^{2m} = (g(chi)/(1-zeta)^m)` is principal and `[A_p]^2 = 1`.  So `W_p` is the ideal lattice of a class of order
`<= 2`; if the class is trivial (e.g. p = 5, 13, where `h(L) = 1, 3` — the value 3 for `Q(zeta_52)` is quoted from
memory), `W_p ~ (Z[zeta_{4p}], Tr(beta y ybar))` for a
totally positive `beta` with `(beta) = (1-zeta)/(2p)`.  Whether `[P^{(p-1)/4}]` is trivial for all `p` is not known
to me.

**Corollary 4.2 (elementary form; PROVED given Theorem 1 and the Round-11 normalisation; NUMERICAL check p = 5, 13
in `mw_poly_model.py`).**  Let `iota in F_p`, `iota^2 = -1`.  Let `Lambda_poly` be the group of pairs
`(c,c') in Z^{F_p} x Z^{F_p}` modulo the two constant vectors, such that mod `p`

      c_k = f(k),  c'_k = f'(k),   f, f' in F_p[X],   deg f, deg f' <= (3m-1)/2,   deg(f - iota f') <= (m-1)/2.

Put `V = sum_k (c_k - mean c)^2 + sum_k (c'_k - mean c')^2` and `S = c^T Q c + c'^T Q c'` with the Paley matrix
`Q = ((j-l)/p)`.  Then `(Lambda_poly, V/(2p)) ~ (Hom(Jac C_p, E), deg)` (the model `L_0`, `x_l = (1/p) sum c^{(l)}_k zeta^k`),
`N(q_F) = (V^2 - S^2/p)/(16 p^2)`, and

      (A)  <=>  V^2 - S^2/p >= p^2 (p-1)^2  for all nonzero classes;
      "every nonconstant map C_p -> E has degree >= (p-1)/2"  <=>  V >= p(p-1)  for all nonzero classes,

with equality for `(p e_k, 0)`, `(0, p e_k)`.  (Proof: `y = sum c_k zeta^k in (1-zeta)^e`, `e <= p-1`, iff
`sum_k c_k k^l = 0 mod p` for `l < e`, iff `c mod p` is the value table of a polynomial of degree `<= p-1-e`;
`Tr(zeta^k) = p[k=0]-1`, `Tr(sqrt p zeta^k) = p (k/p)`.)  Check (exact Gram matrices): p = 13: rank 24, `2 deg`
integral and even with determinant `13^12`, minimal degree 6 attained by exactly 13104 maps, all with `N(q_F) = 9`;
26208 maps of degree 7, all with `N(q_F) = 9` (the types `(7, -+1)`): agrees with Rounds 5 and 11.

**Proposition 4.3 (PROVED given 4.2; degree bound on a larger sublattice).**  If `f` (or `f'`) is non-constant of
degree `<= D` then `V >= (p-D)^3/(12 D^2)`.  Consequently every nonzero `x in (1-zeta)^{-D} Z[zeta_{4p}] ∩ A_0` has
`deg(x) = Tr(x xbar)/4 >= (p-1)/2` as soon as `(p-D)^3 >= 12 D^2 p (p-1)` (roughly `D <= 0.28 sqrt p`).
*Proof.*  A non-constant polynomial of degree `<= D` takes each value at most `D` times, so for every real `t` at most
`D(2R+1)` indices have `|c_k - t| < R`; hence `sum (c_k-t)^2 >= int_0^{R_0} (p - D(2u+1)) 2u du = (p-D)^3/(12D^2)`,
`R_0 = (p/D-1)/2`.  If `f, f'` are both constant the class is `p` times an integral class and `V >= p(p-1)`.  QED
(For the full range `D = (3m-1)/2` this gives only `V >= 0.14 p`; the Castelnuovo-Severi bound is `V >= p(p+3)/2`.
The restriction of the degree question to `x_2 = 0` is the statement that the Craig lattice `A_{p-1}^{((3p+1)/8)}`
has minimal norm exactly `p-1`, attained by the Legendre-symbol vector `sum (a/p) zeta^a = sqrt p`; Craig's bound
gives `(3p+1)/4`.  I do not know whether this is in the literature.)

**Status of (A): REDUCED-TO (A') (equivalently (A''), or the elementary form 4.2); proved on the sublattice `O_L` only.**  Exactly where it is stuck: one needs a lower bound for
`T_+(x)T_-(x)` on the ideal `A_p` that beats the norm bound `N(x) >= p` by a factor `p^{1-1/r}`; all arguments above
are "generic in the ideal class" or use only the minimum, while (A') is false for neighbouring ideals (e.g. `P1` itself
contains `1-zeta` with `T_+T_- = 4p(p-1) < p(p-1)^2`) and for the very same ideal `A_p` at p = 37 and p = 29.
So a proof must use the specific exponent `(m+1)/2` (equivalently the `S`-symmetry) and the primality of `(p+1)/2`
in an essential way, and must still fail at p = 37.

## 5. (B) and (C) cannot be reduced to (A) by Hermite-type inequalities

**Proposition 5.1 (PROVED; uses only `det(wO) = N(nu(w))^4` from Round 5 Lemma 1).**
(B) is equivalent to:  `N(nu(w)) >= r^2 [W_p ∩ w Delta : w O]^{1/2}` for all `w != 0`.
For every line `U` there is `v in U` with `N(nu(v)) <= (p/4) N(n_U)` (Hermite, `gamma_8 = 2`, at any point of the
t-circle), so (A) with the constant `r^2` replaced by `(p/4) r^2` would imply (B); but
`min_w N(nu(w)) <= (p/16) gamma_{8r}^2 < (p/4) r^2` for every lattice of this type (`r >= 3`), so this route is empty.

**Proposition 5.2 (PROVED).**  Let `U` be a `Delta`-subspace of dimension `j`.  (i) Only `j <= r/2` matters:
`det(W∩U) = det(W∩U^perp)` and `(r/(r-j))^{r-j} <= (r/j)^j` for `j <= r/2`, so `R_{U^perp} <= R_U`.
(ii) If `N(nu(w)) >= c` for all `w`, then `det(W ∩ U)^{1/(4j)} >= 16 c/(p gamma_{8j}^2)`, i.e.
`R_U <= ( p gamma_{8j}^2 r^2 / (16 c j^2) )^{4j}`.  With `c = r^2` (i.e. exactly (A)) and `gamma_{8j} ~ 8j/(2 pi e)`
this is `R_U <= (p/72.9)^{4j(1+o(1))}`: it proves `R_U <= 1` only if `p <= 72` roughly, and for larger `p` the bound
is exponentially far from what is needed.  To get (C)_j from a vector bound one needs
`c >= r^2 p gamma_{8j}^2/(16 j^2)` (~ `0.055 p r^2`), which is the size of the *generic* minimum of `N(nu)`.
So **(A) alone can never suffice**; (B) and (C) are separate conditions on the heights
`H(U) = det(W_p ∩ U)^{1/4}` of all `Delta`-subspaces, and need `H(U) >= (r/j)^{2j}` directly.

*Proofs.*  5.1: `det(wO) = N(nu(w))^4` and `det(W∩U) = det(wO)/[W∩U:wO]^2`, `R_U = r^8/det(W∩U)`.  The Hermite step:
a rank-8 lattice of determinant `N(n_U)^4` (the same at every `t`) has a vector with `Q_t(v) <= 2 N(n_U)^{1/2}` and
`m_+m_- <= Q_t^2/4`.  5.2(i): `f(x) = x log(1/x)` satisfies `f(x) >= f(1-x)` on `[0,1/2]`.  (ii): Hermite in
dimension `8j` and `N(nu) = (p/4)m_+m_- <= (p/16) Q_t^2`.  QED

**Status of (B), (C): NO PROGRESS beyond this negative result.**  In the model of Theorem 1 a line is `x Delta`,
`Delta = F(i) + F(i) j`, `j(x) = g(chi) c'(x)` (`c'`: `i -> -i`, `j^2 = -p`), so (B) is the statement that the
`O`-module `A_p ∩ (F(i)x + F(i) g(chi) c'(x))` has determinant `>= r^8` for all `x`; not investigated further.

## 6. Exact computations with the model (NUMERICAL; exact integer Gram matrices)

Code: `ideal_model.py` builds the exact integer Gram matrices `G` (lattice) and `E` (`= G eta`) of `L_N(p)` from
Theorem 2 (Gaussian periods, construction A, exact rational arithmetic); `N(nu) = (p|w|^4 - 4 s^2)/16`, `s = w^T E w/2`,
is then an exact integer for every vector found.  `scan_exact.py`: (i) all vectors of norm <= 4 at the unimodular point
with their exact `(norm, |s|, N(nu))`; (ii) exhaustive search for vectors with `N(nu') < m'^2` over one full period of
the t-circle (`t in [-log eps, log eps]`, step 0.1, bound `2 sqrt(theta') cosh(0.05)`, `theta' = 4 m'^2/p`), with an
exact rational LLL at the start, exact integer base changes and 130-digit evaluation of
`Q_t = cosh(t) G + sinh(t) E/sqrt p` in the current basis (this removes the precision failures of Round 13 for large
`log eps`); the Fincke-Pohst enumeration itself is double precision on a well-conditioned reduced Gram matrix.
"fail = 0" below means that no enumeration was truncated or failed.  The fundamental unit is computed by continued
fractions.  What is exact: the lattice, `G`, `E`, the invariants of every vector found.  What is not certified: the
completeness of the floating-point enumeration (margin `1e-9`, as in Rounds 10-13).

Results (table in Section 8):
* **Model check.**  p = 5: `E_8`.  p = 13: even unimodular, no roots, 98280 pairs of norm 4 with
  `|s| = 0, 2, 4` for `39312, 45864, 13104` pairs — the Leech lattice with exactly the Round-5 distribution; t-circle:
  `min N(nu) = 9`, 19656 borderline pairs mod units (Round 10: 19656).  p = 37, N = 12: 144 roots, all `|s| = 3`,
  `N(nu') = 7 < 9` (the Round-10 violator).  p = 61: N = 12: 720 roots, `|s| in {1,3}`; N = 20: no roots, 19800 pairs of
  norm 4 (the extremal 40-dimensional lattice).  p = 157: 144 roots, `|s| = 1`.  p = 277: Leech.  All as in Rounds 10-13.
* **New: N = 12 (rank 24, threshold 9), exhaustive (fail = 0) for all 57 primes `p = 13 mod 24`, `p < 8900`, with
  `(p+1)/2` prime, and in addition for the 28 primes `p = 13 mod 24`, `p < 3000`, with `(p+1)/2` composite:**
  the only violator among these 85 primes is p = 37.  For every other prime the minimum of `N(nu')` is exactly 9 and
  is attained only by the 12 pairs (mod units) of the `E_8`-sublattice `Fix(C_{(p-1)/4})`.  This settles the cases left
  "partial / precision limited" in Round 13 (p = 421, 541, 661, 757; e.g. p = 661: 144 roots with `|s| = 21` exactly,
  `N(nu') = 55`).  Root systems for `p < 3750`, `(p+1)/2` prime: no roots (Leech) for p = 13, 277, 421, 541, 733, 757,
  1093, 1213, 1237, 1381, 1453, 1621, 2341, 2557, 2797, 3253, 3733; 144 roots for 37, 157, 661, 1933, 2917, 3517;
  720 roots for 61, 397, 613, 877, 997, 3061.  Near misses exist: vectors with `N(nu') = 12` (norm 4, e.g. `s = 46` at
  p = 541) occur for p = 541, 1237, 1381, 3733 — the condition is not satisfied with a growing margin on sublattices.
* **New: N = 20 (rank 40, threshold 25)** for p = 61, 421, 541, 661, 1381, 1621, 2341, 3061: no violator, minimum
  exactly 25 (forced).  **N = 28 (rank 56, threshold 49)** for p = 421, 757, 1093, 1933: no violator, minimum exactly
  49 (forced).  **N = 36 (rank 72, threshold 81)** for p = 2341, 2557, 2917 (t-circle only; feasible because the search
  radius `N/sqrt p` is far below the Gaussian-heuristic radius): no violator, minimum exactly 81 (forced).
* Not feasible with the plain enumerator: N = 36 for p <= 757, N = 44, 60, and the full lattices for p = 29, 37, 61
  (node limit hit already for roots in rank 72).  Note: the Round-2 enumerator `fp_enum.c` has a fixed buffer
  `int x[64]`; calls with rank > 64 smash the stack.  `fp_enum_big.c` is a copy with `x[1024]`, used here for all
  runs after 19:10 UTC.  The rank <= 56 results do not depend on this.
* Started but NOT finished when the round ended (timed out, no result, empty logs): N = 36 for p = 1621
  (about 600 t-steps in rank 72), N = 52 for p = 1093 (rank 104), and the targeted search `check29b.py` for the
  known p = 29 violator (norm 6, `s = +-11`, `N(nu) = 35 < 49`, reported in Round 10) in the full rank-56 lattice of the
  formula.  So the statement "(A') fails at p = 29" rests on Round 10, not on a computation in the new model; the
  reduced basis of the p = 29 lattice has diagonal `6^4 8^33 10^19`, consistent with minimum 6.
* Full lattices from the formula (`full61.py`): p = 29, 37, 61 (ranks 56, 72, 120): even, `det G = 1`, `det E = p^{n/2}`.

## 7. p = 1 mod 8: why nothing was obtained

For `r` even the model of Lemma 1.2 fails (`d^{2r} = +1`).  `xi` (`xi^2=-1`) still exists in `Delta`, so some lattice
of the family is an ideal of `Q(zeta_{4p})`, but `d = kappa sigma~` with `N_{L/F(i)}(kappa) = -1`, and `-1` is not the
relative norm of a root of unity or of a real unit; `kappa` is not a unit in general, the ideal has a non-trivial part
outside `p`, and I did not find it explicitly.  For the existing floating-point pipeline (`scan.py`) there is a second
obstacle: `W^{C}` has `Delta`-dimension `2^{a-2} >= 2` for the largest odd-order torus subgroup `C` (`2^a || p-1`), so no
rational vector generates a single copy of `W` (for p = 73 the natural `C_9`-fixed vector generates `W^2`); one needs
an explicit zero divisor of the cyclic algebra `(L'/F, sigma, -1)`, `[L':F] = 4`.  Predictions that can be tested once
a model exists: the smallest quotient curve is `Y^2 = X^9 - X` (genus 4, rank 16, threshold `N(nu') >= 4`) for
`p = 9 mod 16` (p = 73, 313, 457); roots can never violate there, the dangerous values are `N(nu') in {1,2,3}`
(`2` is a norm from `Q(sqrt 73)`: `(9+sqrt 73)/2`), i.e. maps `D -> E` with `q_F` of norm 1, 2, 3; `N = 1` would be a
degree-2 map and is excluded because `Y^2 = X^9 - X` has no bielliptic involution (its extra involutions have genus-2
quotients).  No borderline line is known for `p = 1 mod 8`.

*Two derivations that were NOT verified numerically (sketches, UNVERIFIED):*
(1) For `p = 9 mod 16` the same argument as Lemma 1.2 with `d = zeta_8 sigma~` on `Q(zeta_{8p})` gives a module of
`F`-dimension `8r`, i.e. `W^2`; a `G`-lattice of the form `W_p (x)_O L_2` (`L_2` of rank 2) should be the ideal
`prod_j P_j^{c_j}` of `Q(zeta_{8p})`, `P_j` the four primes over `p`, `c_j = (1 - 2r - e_j r/2)/2 mod 2r`,
`e_j in {1,3,5,7}` (the image of `zeta_8` in `Z_p` is `g^{e_j (p-1)/8}`).  To test (A) one must still isolate the
rank-one tensors.
(2) If `q = (p+1)/2` is prime the non-split torus acts freely, `xi` and `eta` commute with it, and `W` is
one-dimensional over the field `Q(zeta_{4q}, sqrt p)` of degree `2(p-1)`: `W_p` should also be an ideal lattice of
this second field.  The hypothesis "(p+1)/2 prime" enters only here: (A') as stated in Theorem 3 makes sense for all
`p = 5 mod 8`, and it is FALSE for p = 29 (Round 10 found vectors with `Gamma/theta = 5/7`), so any proof of (A') must
use the second structure (equivalently the element `S`), which is invisible in the Borel model.

## 8. Tables of numerical results

(generated by `make_table.py` from the `log_scan*.jsonl` files; "roots" = number of norm-2 vectors of the unimodular
lattice `L_N(p)` with the multiset of `|s|`; "min N" = smallest `N(nu')` found on the whole t-circle and the vector
realising it as (norm, s); "hits" = number of enumerated vectors below the search bound, both scan directions.)

<!--TABLE-->
Summary (° marks primes with (p+1)/2 composite, included only as additional tests of the quotient-curve lattices):

* N = 12 (rank 24, threshold 9): 85 primes searched exhaustively (fail = 0), of which 57 have (p+1)/2 prime; range p = 13..8893; violators: [37]; every other prime has minimum exactly 9 (forced borderline vectors only).  Primes: 13, 37, 61, 109°, 157, 181°, 229°, 277, 349°, 373°, 397, 421, 541, 613, 661, 709°, 733, 757, 829°, 853°, 877, 997, 1021°, 1069°, 1093, 1117°, 1213, 1237, 1381, 1429°, 1453, 1549°, 1597°, 1621, 1669°, 1693°, 1741°, 1789°, 1861°, 1933, 2029°, 2053°, 2221°, 2269°, 2293°, 2341, 2389°, 2437°, 2557, 2677°, 2749°, 2797, 2917, 3061, 3253, 3517, 3733, 4021, 4261, 4357, 4621, 4933, 5077, 5101, 5413, 5437, 5581, 5701, 6037, 6133, 6373, 6637, 6661, 6781, 6997, 7213, 7477, 7933, 8053, 8101, 8221, 8317, 8461, 8677, 8893.
* N = 20 (rank 40, threshold 25): 8 primes searched exhaustively (fail = 0), of which 8 have (p+1)/2 prime; range p = 61..3061; violators: none; every other prime has minimum exactly 25 (forced borderline vectors only).  Primes: 61, 421, 541, 661, 1381, 1621, 2341, 3061.
* N = 28 (rank 56, threshold 49): 4 primes searched exhaustively (fail = 0), of which 4 have (p+1)/2 prime; range p = 421..1933; violators: none; every other prime has minimum exactly 49 (forced borderline vectors only).  Primes: 421, 757, 1093, 1933.
* N = 36 (rank 72, threshold 81): 3 primes searched exhaustively (fail = 0), of which 3 have (p+1)/2 prime; range p = 2341..2917; violators: none; every other prime has minimum exactly 81 (forced borderline vectors only).  Primes: 2341, 2557, 2917.

| N (rank) | p | (p+1)/2 prime | k | roots (pairs by abs s) | norm-4 pairs | log eps | t-steps | fail | hits | min N(nu') on circle [norm, s] | threshold | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 (24) | 13 | yes | 1 | 0 {} | 98280 | 1.19 | 25 | 0 | 26208 | 9 [4, -4] | 9 | borderline only |
| 12 (24) | 37 | yes | 3 | 144 {'3': 72} | 96552 | 2.49 | 51 | 0 | 744 | 7 [2, -3] | 9 | **VIOLATION** |
| 12 (24) | 61 | yes | 5 | 720 {'3': 72, '1': 288} | 89640 | 3.66 | 75 | 0 | 24 | 9 [30, -117] | 9 | borderline only |
| 12 (24) | 109 | no | 9 | 0 {} | 98280 | 5.56 | 113 | 0 | 24 | 9 [150, -783] | 9 | borderline only |
| 12 (24) | 157 | yes | 13 | 144 {'1': 72} | 96552 | 5.36 | 109 | 0 | 24 | 9 [102, -639] | 9 | borderline only |
| 12 (24) | 181 | no | 15 | 0 {} | 98280 | 7.17 | 145 | 0 | 24 | 9 [582, -3915] | 9 | borderline only |
| 12 (24) | 229 | no | 19 | 0 {} | 98280 | 2.71 | 55 | 0 | 24 | 9 [6, -45] | 9 | borderline only |
| 12 (24) | 277 | yes | 23 | 0 {} | 98280 | 7.87 | 159 | 0 | 24 | 9 [942, -7839] | 9 | borderline only |
| 12 (24) | 349 | no | 29 | 720 {'3': 360} | 89640 | 9.82 | 197 | 0 | 24 | 9 [5916, -55260] | 9 | borderline only |
| 12 (24) | 373 | no | 31 | 720 {'3': 288, '15': 72} | 89640 | 9.23 | 185 | 0 | 24 | 9 [3180, -30708] | 9 | borderline only |
| 12 (24) | 397 | yes | 33 | 720 {'5': 72, '1': 288} | 89640 | 8.15 | 163 | 0 | 24 | 9 [1038, -10341] | 9 | borderline only |
| 12 (24) | 421 | yes | 35 | 0 {} | 98280 | 13.01 | 261 | 0 | 24 | 9 [130110, -1334817] | 9 | borderline only |
| 12 (24) | 541 | yes | 45 | 0 {} | 98280 | 14.15 | 283 | 0 | 24 | 9 [360222, -4189275] | 9 | borderline only |
| 12 (24) | 613 | yes | 51 | 720 {'3': 288, '9': 72} | 89640 | 11.50 | 231 | 0 | 24 | 9 [23934, -296289] | 9 | borderline only |
| 12 (24) | 661 | yes | 55 | 144 {'21': 72} | 96552 | 14.40 | 289 | 0 | 24 | 9 [417630, -5368617] | 9 | borderline only |
| 12 (24) | 709 | no | 59 | 0 {} | 98280 | 17.41 | 349 | 0 | 24 | 9 [8222604, -109471860] | 9 | borderline only |
| 12 (24) | 733 | yes | 61 | 0 {} | 98280 | 3.30 | 67 | 0 | 24 | 9 [6, -81] | 9 | borderline only |
| 12 (24) | 757 | yes | 63 | 0 {} | 98280 | 14.82 | 297 | 0 | 24 | 9 [597228, -8215956] | 9 | borderline only |
| 12 (24) | 829 | no | 69 | 0 {} | 98280 | 17.25 | 345 | 0 | 24 | 9 [6455580, -92935692] | 9 | borderline only |
| 12 (24) | 853 | no | 71 | 144 {'27': 72} | 96552 | 10.22 | 205 | 0 | 24 | 9 [5646, -82449] | 9 | borderline only |
| 12 (24) | 877 | yes | 73 | 720 {'1': 288, '21': 72} | 89640 | 13.09 | 263 | 0 | 24 | 9 [97788, -1447956] | 9 | borderline only |
| 12 (24) | 997 | yes | 83 | 720 {'1': 288, '7': 72} | 89640 | 12.04 | 241 | 0 | 24 | 9 [32268, -509436] | 9 | borderline only |
| 12 (24) | 1021 | no | 85 | 144 {'1': 72} | 96552 | 18.27 | 367 | 0 | 24 | 9 [16100958, -257237685] | 9 | borderline only |
| 12 (24) | 1069 | no | 89 | 144 {'25': 72} | 96552 | 18.49 | 371 | 0 | 24 | 9 [19603110, -320467383] | 9 | borderline only |
| 12 (24) | 1093 | yes | 91 | 0 {} | 98280 | 3.50 | 71 | 0 | 24 | 9 [6, -99] | 9 | borderline only |
| 12 (24) | 1117 | no | 93 | 0 {} | 98280 | 15.83 | 317 | 0 | 24 | 9 [1343670, -22453767] | 9 | borderline only |
| 12 (24) | 1213 | yes | 101 | 0 {} | 98280 | 18.77 | 377 | 0 | 24 | 9 [24325404, -423604404] | 9 | borderline only |
| 12 (24) | 1237 | yes | 103 | 0 {} | 98280 | 14.07 | 283 | 0 | 24 | 9 [220758, -3882141] | 9 | borderline only |
| 12 (24) | 1381 | yes | 115 | 0 {} | 98280 | 25.05 | 501 | 0 | 24 | 9 [12174108630, -226205945883] | 9 | borderline only |
| 12 (24) | 1429 | no | 119 | 144 {'27': 72} | 96552 | 5.24 | 105 | 0 | 24 | 9 [30, -567] | 9 | borderline only |
| 12 (24) | 1453 | yes | 121 | 0 {} | 98280 | 21.84 | 437 | 0 | 24 | 9 [481649766, -9179819991] | 9 | borderline only |
| 12 (24) | 1549 | no | 129 | 144 {'9': 72} | 96552 | 27.24 | 545 | 0 | 24 | 9 [103196513766, -2030770000665] | 9 | borderline only |
| 12 (24) | 1597 | no | 133 | 0 {} | 98280 | 18.43 | 369 | 0 | 24 | 9 [15111150, -301939533] | 9 | borderline only |
| 12 (24) | 1621 | yes | 135 | 0 {} | 98280 | 29.20 | 585 | 0 | 24 | 9 [718841301342, -14470866383625] | 9 | borderline only |
| 12 (24) | 1669 | no | 139 | 0 {} | 98280 | 20.98 | 421 | 0 | 24 | 9 [189949974, -3880050795] | 9 | borderline only |
| 12 (24) | 1693 | no | 141 | 144 {'15': 72} | 96552 | 14.13 | 283 | 0 | 24 | 9 [200190, -4118517] | 9 | borderline only |
| 12 (24) | 1741 | no | 145 | 720 {'5': 288, '21': 72} | 89640 | 25.18 | 505 | 0 | 24 | 9 [12350747478, -257669273025] | 9 | borderline only |
| 12 (24) | 1789 | no | 149 | 720 {'7': 288, '29': 72} | 89640 | 27.03 | 541 | 0 | 24 | 9 [77863159686, -1646672368545] | 9 | borderline only |
| 12 (24) | 1861 | no | 155 | 720 {'5': 288, '39': 72} | 89640 | 22.05 | 441 | 0 | 24 | 9 [523287948, -11287140300] | 9 | borderline only |
| 12 (24) | 1933 | yes | 161 | 144 {'25': 72} | 96552 | 20.52 | 411 | 0 | 24 | 9 [110875206, -2437363881] | 9 | borderline only |
| 12 (24) | 2029 | no | 169 | 0 {} | 98280 | 3.81 | 77 | 0 | 24 | 9 [6, -135] | 9 | borderline only |
| 12 (24) | 2053 | no | 171 | 720 {'9': 72, '7': 288} | 89640 | 24.06 | 483 | 0 | 24 | 9 [3732164166, -84552259491] | 9 | borderline only |
| 12 (24) | 2221 | no | 185 | 0 {} | 98280 | 30.05 | 601 | 0 | 24 | 9 [1429823174460, -33691986298692] | 9 | borderline only |
| 12 (24) | 2269 | no | 189 | 144 {'15': 72} | 96552 | 33.13 | 663 | 0 | 24 | 9 [30865104613020, -735114527677692] | 9 | borderline only |
| 12 (24) | 2293 | no | 191 | 144 {'27': 72} | 96552 | 24.52 | 491 | 0 | 24 | 9 [5602842174, -134146831959] | 9 | borderline only |
| 12 (24) | 2341 | yes | 195 | 0 {} | 98280 | 30.48 | 611 | 0 | 24 | 9 [2144998190220, -51891669248508] | 9 | borderline only |
| 12 (24) | 2389 | no | 199 | 144 {'21': 72} | 96552 | 33.31 | 667 | 0 | 24 | 9 [35775158691414, -874298328486435] | 9 | borderline only |
| 12 (24) | 2437 | no | 203 | 144 {'45': 72} | 96552 | 21.42 | 429 | 0 | 24 | 9 [243677190, -6014681667] | 9 | borderline only |
| 12 (24) | 2557 | yes | 213 | 0 {} | 98280 | 9.55 | 191 | 0 | 24 | 9 [1662, -42021] | 9 | borderline only |
| 12 (24) | 2677 | no | 223 | 0 {} | 98280 | 8.24 | 165 | 0 | 24 | 9 [438, -11331] | 9 | borderline only |
| 12 (24) | 2749 | no | 229 | 720 {'5': 288, '19': 72} | 89640 | 25.49 | 511 | 0 | 24 | 9 [13435107870, -352207445067] | 9 | borderline only |
| 12 (24) | 2797 | yes | 233 | 0 {} | 98280 | 28.43 | 569 | 0 | 24 | 9 [252670324668, -6681446202564] | 9 | borderline only |
| 12 (24) | 2917 | yes | 243 | 144 {'35': 72} | 96552 | 4.68 | 95 | 0 | 24 | 9 [12, -324] | 9 | borderline only |
| 12 (24) | 3061 | yes | 255 | 720 {'31': 72, '1': 288} | 89640 | 40.05 | 801 | 0 | 24 | 9 [26746649251671990, -739896634097224317] | 9 | borderline only |
| 12 (24) | 3253 | yes | 271 | 0 {} | 98280 | 4.04 | 81 | 0 | 24 | 9 [6, -171] | 9 | borderline only |
| 12 (24) | 3517 | yes | 293 | 144 {'15': 72} | 96552 | 17.71 | 355 | 0 | 24 | 9 [4948590, -146736333] | 9 | borderline only |
| 12 (24) | 3733 | yes | 311 | 0 {} | 98280 | 29.75 | 595 | 0 | 24 | 9 [815109590454, -24900897597939] | 9 | borderline only |
| 12 (24) | 4021 | yes | 335 | 0 {} | 98280 | 33.89 | 679 | 0 | 24 | 9 [49541232082092, -1570738338493020] | 9 | borderline only |
| 12 (24) | 4261 | yes | 355 | 720 {'57': 72, '1': 288} | 89640 | 41.44 | 829 | 0 | 24 | 9 [91538447978800110, -2987647229012483817] | 9 | borderline only |
| 12 (24) | 4357 | yes | 363 | 144 {'3': 72} | 96552 | 4.88 | 99 | 0 | 24 | 9 [12, -396] | 9 | borderline only |
| 12 (24) | 4621 | yes | 385 | 144 {'9': 72} | 96552 | 47.68 | 955 | 0 | 24 | 9 [45003360852726825558, -1529617827860629724385] | 9 | borderline only |
| 12 (24) | 4933 | yes | 411 | 0 {} | 98280 | 7.78 | 157 | 0 | 24 | 9 [204, -7164] | 9 | borderline only |
| 12 (24) | 5077 | yes | 423 | 144 {'21': 72} | 96552 | 28.14 | 563 | 0 | 24 | 9 [139820573100, -4981322548908] | 9 | borderline only |
| 12 (24) | 5101 | yes | 425 | 0 {} | 98280 | 37.41 | 749 | 0 | 24 | 9 [1477452312880062, -52760771609202045] | 9 | borderline only |
| 12 (24) | 5413 | yes | 451 | 0 {} | 98280 | 39.14 | 783 | 0 | 24 | 9 [8138635831993614, -299392305104461401] | 9 | borderline only |
| 12 (24) | 5437 | yes | 453 | 0 {} | 98280 | 29.39 | 589 | 0 | 24 | 9 [470464013742, -17345070933741] | 9 | borderline only |
| 12 (24) | 5581 | yes | 465 | 144 {'23': 72} | 96552 | 55.65 | 1113 | 0 | 24 | 9 [118141565713260054459510, -4412947266679509458432433] | 9 | borderline only |
| 12 (24) | 5701 | yes | 475 | 0 {} | 98280 | 42.78 | 857 | 0 | 24 | 9 [301153073806793310, -11369276408313753183] | 9 | borderline only |
| 12 (24) | 6037 | yes | 503 | 720 {'11': 288, '65': 72} | 89640 | 29.15 | 585 | 0 | 24 | 9 [354056424222, -13754761715529] | 9 | borderline only |
| 12 (24) | 6133 | yes | 511 | 720 {'17': 72, '13': 288} | 89640 | 12.08 | 243 | 0 | 24 | 9 [13494, -528381] | 9 | borderline only |
| 12 (24) | 6373 | yes | 531 | 720 {'13': 288, '9': 72} | 89640 | 39.70 | 795 | 0 | 24 | 9 [13076864670172620, -521970060188705292] | 9 | borderline only |
| 12 (24) | 6637 | yes | 553 | 144 {'9': 72} | 96552 | 14.18 | 285 | 0 | 24 | 9 [105510, -4297833] | 9 | borderline only |
| 12 (24) | 6661 | yes | 555 | 0 {} | 98280 | 51.46 | 1031 | 0 | 24 | 9 [1648510262195378440902, -67271540991115458213405] | 9 | borderline only |
| 12 (24) | 6781 | yes | 565 | 0 {} | 98280 | 49.61 | 993 | 0 | 24 | 9 [254698874795421666510, -10486822154141207896317] | 9 | borderline only |
| 12 (24) | 6997 | yes | 583 | 0 {} | 98280 | 15.11 | 303 | 0 | 24 | 9 [262230, -10967517] | 9 | borderline only |
| 12 (24) | 7213 | yes | 601 | 144 {'25': 72} | 96552 | 17.80 | 357 | 0 | 24 | 9 [3816444, -162064116] | 9 | borderline only |
| 12 (24) | 7477 | yes | 623 | 0 {} | 98280 | 35.67 | 715 | 0 | 24 | 9 [215877312110238, -9333417578923431] | 9 | borderline only |
| 12 (24) | 7933 | yes | 661 | 720 {'13': 288, '45': 72} | 89640 | 41.60 | 833 | 0 | 24 | 9 [78551687634844926, -3498196907373464421] | 9 | borderline only |
| 12 (24) | 8053 | yes | 671 | 0 {} | 98280 | 34.30 | 687 | 0 | 24 | 9 [52914900072030, -2374252118536167] | 9 | borderline only |
| 12 (24) | 8101 | yes | 675 | 144 {'35': 72} | 96552 | 5.19 | 105 | 0 | 24 | 9 [12, -540] | 9 | borderline only |
| 12 (24) | 8221 | yes | 685 | 144 {'63': 72} | 96552 | 45.71 | 915 | 0 | 24 | 9 [4690220751721517142, -212630525397758047095] | 9 | borderline only |
| 12 (24) | 8317 | yes | 693 | 720 {'67': 72, '1': 288} | 89640 | 51.75 | 1035 | 0 | 24 | 9 [1956822388553444180718, -89228740684379690817171] | 9 | borderline only |
| 12 (24) | 8461 | yes | 705 | 720 {'15': 288, '27': 72} | 89640 | 42.76 | 857 | 0 | 24 | 9 [242526392558861052, -11154236765152385340] | 9 | borderline only |
| 12 (24) | 8677 | yes | 723 | 0 {} | 98280 | 45.23 | 905 | 0 | 24 | 9 [2841836203943944782, -132359112293198141799] | 9 | borderline only |
| 12 (24) | 8893 | yes | 741 | 720 {'77': 72, '13': 288} | 89640 | 43.26 | 867 | 0 | 24 | 9 [391722652946412534, -18470252714442098679] | 9 | borderline only |
| 20 (40) | 61 | yes | 3 | 0 {} | 19800 | 3.66 | 75 | 0 | 24 | 25 [50, -195] | 25 | borderline only |
| 20 (40) | 421 | yes | 21 | 0 {} | 19800 | 13.01 | 261 | 0 | 24 | 25 [216850, -2224695] | 25 | borderline only |
| 20 (40) | 541 | yes | 27 | 240 {'1': 120} | 74520 | 14.15 | 283 | 0 | 24 | 25 [600370, -6982125] | 25 | borderline only |
| 20 (40) | 661 | yes | 33 | 0 {} | 19800 | 14.40 | 289 | 0 | 24 | 25 [696050, -8947695] | 25 | borderline only |
| 20 (40) | 1381 | yes | 69 | 0 {} | 19800 | 25.05 | 501 | 0 | 24 | 25 [20290181050, -377009909805] | 25 | borderline only |
| 20 (40) | 1621 | yes | 81 | 240 {'11': 120} | 74520 | 29.20 | 585 | 0 | 24 | 25 [1198068835570, -24118110639375] | 25 | borderline only |
| 20 (40) | 2341 | yes | 117 | 0 {} | 19800 | 30.48 | 611 | 0 | 24 | 25 [3574996983700, -86486115414180] | 25 | borderline only |
| 20 (40) | 3061 | yes | 153 | 0 {} | 19800 | 40.05 | 801 | 0 | 24 | 25 [44577748752786650, -1233161056828707195] | 25 | borderline only |
| 28 (56) | 421 | yes | 15 | 0 {} | 0 | 13.01 | 261 | 0 | 24 | 49 [303590, -3114573] | 49 | borderline only |
| 28 (56) | 757 | yes | 27 | 0 {} | 0 | 14.82 | 297 | 0 | 24 | 49 [1393532, -19170564] | 49 | borderline only |
| 28 (56) | 1093 | yes | 39 | 0 {} | 0 | 3.50 | 71 | 0 | 24 | 49 [14, -231] | 49 | borderline only |
| 28 (56) | 1933 | yes | 69 | n/c | n/c | 20.52 | 411 | 0 | 24 | 49 [258708814, -5687182389] | 49 | borderline only |
| 36 (72) | 2341 | yes | 65 | n/c | n/c | 30.48 | 611 | 0 | 24 | 81 [6434994570660, -155675007745524] | 81 | borderline only |
| 36 (72) | 2557 | yes | 71 | n/c | n/c | 9.55 | 191 | 0 | 24 | 81 [4986, -126063] | 81 | borderline only |
| 36 (72) | 2917 | yes | 81 | n/c | n/c | 4.68 | 95 | 0 | 24 | 81 [36, -972] | 81 | borderline only |
<!--/TABLE-->

## 9. Where exactly things are stuck, and next steps

1. (A): one needs `T_+(x)T_-(x) >= p(p-1)^2` on the ideal `A_p`; proved only for `x` with `x xbar` `p`-integral
   (Theorem 5).  Every available inequality (AM-GM + integrality, Craig, Hermite, unit orbits, the quadratic-residue
   trick, Siegel-Weil over the genus) is insensitive to the exponent `(m+1)/2` and to the primality of `(p+1)/2`, while
   the statement is false without them (p = 29, p = 37).  Most promising: make `S` (or the non-split torus of order
   `p+1`) explicit in the ideal model — `W_p` then carries two cyclotomic structures (over `Q(zeta_{4p})` and over
   `Q(zeta_{2(p+1)}, i, sqrt p)`), and (A) should be attacked with both, as in Tits' and Craig-type descriptions of the
   Leech lattice.  A second concrete sub-problem: is the Craig lattice `A_{p-1}^{((3p+1)/8)}` of minimal norm exactly
   `p-1` (Proposition 4.3)?  A third: is the class `[P^{(p-1)/4}]` trivial (Remark 4.1)?
2. Numerics: with BKZ + pruned enumeration the exact Gram matrices of `ideal_model.py` make rank 72-120 sublattices
   (`N = 36, 44, 60`) and possibly the full `W_61` accessible; the search radius is only `(8.5/sqrt p)` times the
   Gaussian-heuristic radius squared for the sublattice of the largest primes.
3. (B), (C): compute, in the model, the lattices `A_p ∩ x Delta` (`Delta = F(i) + F(i) g(chi) c'`) for the short vectors
   found, i.e. test non-principal lines in `L_12(p)`, `L_20(p)`: by Hermite a line with `R_U > 1` in `L_N(p)` contains a
   vector with `N(nu') < (p/4) m'^2`, so the enumeration radius must be raised by the factor `sqrt(p)/2`.
4. `p = 1 mod 8`: find the explicit zero divisor / the ideal of `Q(zeta_{8p})` (sketch (1) above) and run the same scan
   for `Y^2 = X^9 - X` (p = 73, 313, 457).

## 10. Files

* `proof.md` — this file.
* `ideal_model.py` — exact construction of `L_N(p)` (both types `L_1`, `L_0`), exact determinant, basic analysis.
* `scan_exact.py` — exact scan (roots, short vectors, exhaustive t-circle search).  `make_table.py` — table.
* `mw_poly_model.py`, `log_mw_poly_model.txt` — the elementary polynomial form (Corollary 4.2) checked at p = 5, 13.
* `t_L0.py` — comparison of the `L_0`- and `L_1`-type lattices; `t_lll.py` — timing of the exact LLL;
  `full61.py`, `log_full_lattices.txt` — full lattices for p = 29, 37, 61 from the formula.
* `fp_enum_big.c/.so` — copy of the Round-2 enumerator with a larger coordinate buffer.
* Logs: `log_verify_model_A.txt` (first verification run, old t-circle routine), `log_scanA_N12.jsonl`,
  `log_scanA2_N12.jsonl`, `log_scanA3_N12.jsonl` (N = 12), `log_scanB*_N20.jsonl` (N = 20), `log_scanC*_*.jsonl`
  (N = 28, 36), `log_scanD*` (rank >= 88, failed: see Section 6), `log_scanE*_N12_allp.jsonl` (N = 12 for the primes
  `p = 13 mod 24` with `(p+1)/2` composite), `job*.txt` (argument lists), `ideal_model_results_*.json`.
