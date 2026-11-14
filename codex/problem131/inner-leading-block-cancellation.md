# An actual leading-block cancellation for the ζ(7) functional

Status: exact local functional identities and checked rational matrix ranks.
Section 5 gives a proposed global valuation improvement with its detailed
lifting argument for independent audit. It is much too small to close the
current ζ(7) budget. This is not an irrationality proof or a historical
novelty claim.

## 1. Derive the block from the actual residue substitution

Fix an ordinary class 0<a<p/2 and write x=K/p, z=a/p. Under the substitution
of the source variable v=a+pu, the near integer poles are precisely

\[
 I_K=[-\lfloor x+z\rfloor,\,\lfloor x-z\rfloor]\cap\mathbb Z,
 \quad I_N=[-\lfloor\alpha x+z\rfloor,\,\lfloor\alpha x-z\rfloor]
                 \cap\mathbb Z.
\]

Set U_K(u)=∏_{j∈I_K}(u−j), U_N(u)=∏_{j∈I_N}(u−j). A near factor
j²−v² contributes p(−2a)(u−j_local), up to a unit congruent to one modulo p
after extracting its constant. A CRT row with order i at this class
contributes p^i(−2a)^i u^i times a common unit. Thus after the row valuation
factors are removed, the self-source block is, up to common and diagonal units,

\[
 B_{ij}=\tau_{s,0}\left(u^{i+j}\frac{U_N(u)^r}{U_K(u)}\right). \tag{1}
\]

Here τ_{s,0}(1/(u-j))=H_{d(j)}^{(s)} and
τ_s(P)=L(P^{(s−2)})/(s−1)!. This uses Y=p^sX+C_{p,s}≡0 mod p from the
distribution formula. It is a computation of the leading reduction, not a
replacement of the whole p-adic functional by τ_{s,0}.

There is an extra reduction issue beyond the source's integrality lemma.
That lemma allows the analytic functional to lose one power of p, so one
cannot simply discard pU₁ without checking its degree. In the fixed patterns
here, the first expansion coefficient has degree at most deg U₀+1=O(M),
because one first-order far-factor or near-unit correction adds at most one
power of u. This is below p+1 for sufficiently large p. Polynomial division
and the unit differences between near poles therefore make τ(U₁/U_K)
integral. All terms p^jU_j with j≥2 remain divisible by p even after the
possible one-power loss. This justifies using (1) for the leading reduction.

The paired source v=p−a+pu equals −(a+p(−1−u)). The pulled-back function
v^sR(−v²) is odd, while τ_s(g(−1−u))=−τ_s(g(u)). Therefore the two sources
**add**, giving a factor 2; they do not cancel by reflection.

Other classes can contribute to border entries of the fully scaled matrix,
especially when extra rows were assigned between tied class sizes. For this
reason a singular self block alone does not prove a determinant gain. We also
test the rectangular block with one extra moment column.

## 2. The annihilator is explicit

The functional has two zero pole values:
τ_{s,0}(1/u)=τ_{s,0}(1/(u+1))=0. It also annihilates polynomials of degree
below s−2. Therefore

\[
 \boxed{\tau_{s,0}\left(\frac{u^j}{u(u+1)}\right)=0,
                  \qquad 0\le j\le s-1.} \tag{2}
\]

Indeed polynomial division has quotient degree at most j−2<s−2, and its
remaining poles are 0 and −1, whose functional values vanish.

For s=7, r=6, x=3, and z>9/40, I_N is empty and
I_K={−3,−2,−1,0,1,2}. The base allocation has L=6 rows. Put

\[
 q(u)=\frac{U_K(u)}{u(u+1)}=(u+3)(u+2)(u-1)(u-2).
\]

The coefficient vector of q, padded to length six, annihilates the moment
columns u^j for 0≤j≤6. Thus it lies in the left kernel of the actual leading
6×7 block. Its degree is four, strictly below the top row degree five.

More generally, when I_N is empty and I_K contains 0,−1, the polynomials
u^k U_K/[u(u+1)] give explicit kernel vectors whenever both their row degree
fits and k plus the tested column degree is at most s−1. This yields a direct
rank bound, rather than merely a numerical determinant coincidence.

## 3. Exact finite comparisons

`python codex/computations/check_inner_leading_blocks.py` uses exact rational
Bernoulli numbers, polynomial division, harmonic sums and matrix ranks. It
corrects SymPy's B₁=+1/2 convention to the required B₁=−1/2. Six representative
floor patterns are tested for each of (s,r)=(5,6),(7,6),(7,8), with both the
base row count and one extra row. Selected results are:

| s,r | x,z | I_N | L | self rank | rank with extra column |
|---|---|---|---:|---:|---:|
| 5,6 | 3,2/5 | empty | 6 | 6 | 6 |
| 7,6 | 3,2/5 | empty | 6 | 4 | 5 |
| 7,6 | 3,2/5 | empty | 7 | 6 | 7 |
| 7,8 | 3,2/5 | empty | 7 | 6 | 7 |
| 7,6 | 3,1/10 | {0} | 3 | 3 | 3 |
| 7,6 | 31/10,2/5 | empty | 7 | 6 | 7 |
| 7,6 | 4,2/5 | empty | 9 | 8 | 9 |

The larger odd-dimensional self blocks can be singular by reflection parity,
but their border columns restore full row rank. The robust 6×7 defect is the
new feature relevant to a possible global improvement. Nonzero rational minors
remain nonzero modulo all sufficiently large primes within a fixed floor
pattern; no conclusion for all patterns is inferred from this finite list.

## 4. Couplings to other classes: why one extra column suffices here

Retain α=3/40, λ=37/40, r=6 and H=λ+3α=23/20. In a class with no N-pole,
ℓ_K=6 and no extra row, the weights for s=7 are w_i=i−6, 0≤i≤5.
All are integers. At its own source, every row belonging to another CRT block
has vanishing order L=6. Thus any potentially leading external column is a
unit multiple of the moment column u⁶. Identity (2) kills it as well.

At any other ordinary source c, a row from our selected class has order L_c.
Compared to its assigned weight, the excess half-weight is

\[
 Z_c-Z_a+(L_a-i)+(\ell_a-\ell_c)/2.
\]

Here Z_a=T, Z_c≥T, ℓ_c is 6 or 7, and q uses only i≤4 while L_a=6.
The excess is therefore at least 3/2. Such sources cannot create a nonzero
leading entry against an integer-parity column. The reserved zero source has
the strict surplus supplied by its large L₀ in the source proof. A zero-block
row at our ordinary source again contributes a column of degree six.

This is more than self-block singularity: the explicit vector annihilates all
leading entries in the integer-weight sector of the global scaled matrix.

## 5. A proposed prime gain on an interval, not just at one point

For

\[
 3<x<\frac{70}{23},
\]

the base allocation has T=6. Classes with z<x−3 have ℓ_K=7 and receive extras
first. All other classes have ℓ_K=6; ties may be assigned arbitrarily without
altering the old row-weight sum. Choose these ties so that every **nonextra**
class has I_N empty. There are enough such classes: their available measure
is 1/2−αx, while the number of nonextra classes divided by p is

\[
 g(x)=\frac72-Hx=\frac72-\frac{23x}{20}.
\]

Thus p g(x)+O(1) disjoint classes have the robust null vector above. The
O(1) includes the reserved zero block and endpoint rounding; M stays fixed.

The half-integer scaling needs care. Factor the original row valuations in
the ramified extension with √p. Each selected kernel lies entirely among rows
whose weights are integers. After unimodular row/column operations within
that parity sector, its entries against other integer-parity columns are
divisible by p; entries against half-integer-parity columns have valuation at
least 1/2 automatically. Entries between two selected kernel rows are
divisible by p. Consequently every selected row can receive an additional
half-weight 1/2, giving **one full p-power per kernel row** in the determinant
bound. Counting only reduction rank over the ramified residue field would
give a weaker half-power and miss this parity argument.

Subject to checking the stated local comparison against the source's entire
finite-K basis, the proposed improved inner envelope is

\[
 \gamma^{in}_{p,7,\rm new}\ge\gamma^{in}_{p,7,\rm old}
                              +p\,g(K/p)+O_M(1)
\]

on this interval. The corresponding asymptotic arithmetic improvement is

\[
 \int_3^{70/23}\frac{g(x)}{x^3}\,dx
 =\frac{(6H-7)^2}{252}=\boxed{\frac1{25200}}.
\]

This is a small but structurally meaningful target: the row-weight bound is
not sharp for the actual higher-odd functional, although it was sharp for an
arbitrary matrix with those weights. Independent review of the global lifting
step remains appropriate before promoting this to a fully checked local lemma.

## 6. Scope and next test

The proposed saving is about 0.00003968, compared with a remaining budget gap
of order 0.9 in the bounded ζ(7) experiments. It does not come close to a proof.
Its significance is the exact mechanism: vanishing Bernoulli polynomial
moments and the two zero harmonic pole values create extra leading rank loss.
The next useful extension would classify these annihilators on other floor
patterns and quantify their coupling to border rows before doing more searches.
Adding r=s+1 did not preserve this particular robust defect in the tested
pattern, because it changed the allocated row count from six to seven.
