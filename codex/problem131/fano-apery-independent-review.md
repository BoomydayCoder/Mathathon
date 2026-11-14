# Independent review: Fano V12 / Apéry / Gamma dictionary

Reviewer: simplify agent, independently of the dictionary author.
Status: accepted for the elementary algebra and the stated, attributed
primary-source dictionary, with one terminology refinement communicated to
the author. Deep geometric theorems are checked as source attributions;
this audit does not reprove mirror symmetry, regulator theory or the Gamma
conjecture.

## Material examined

- `fano-v12-apery-gamma-dictionary.md` in full.
- `check_fano_apery_dictionary.py`, independently rerun: PASS.
- Archived Golyshev–Zagier 2016 text, notably the matrix on p.29,
  §§1.4–1.5 and §§2.1–2.2.
- Archived Golyshev–Kerr–Sasaki arXiv:2009.14762v3 text, Theorem 5.1,
  §5.4 and Example 6.5. The rendered PDF page 28 was also viewed directly,
  including all three signs in its displayed birational map.
- New independent full-Jacobian and real-domain checker
  `check_reciprocal_region.py`: PASS.

## Proof checks

**Reciprocal map.** For C=(1+x)(1+y), Q=C+xyz,
\(T=(-z/(1+z),-C/Q,-x/(1+x))\), direct substitution gives
\(T^2=1\) and \(\Lambda(T)\Lambda=1\). Independently computing the full
ordinary Jacobian and then multiplying by \(xyz/(T_1T_2T_3)\) gives exactly
\(1/\Lambda\). This checks the sign and exponent shift rather than relying
only on the author's one-permutation logarithmic determinant.

The printed GKS v3 map really has the three signs stated in the dictionary.
Its value at (−1,−1,−1) transforms φ=40 into 2349/40, not 1/40.
Changing all three coordinate signs produces the stated involution.
This is a valid localized correction to the retrieved version. It is not
an audit of all later editions or of every step of the source's cycle
construction. The dictionary appropriately retains the geometric
interpretation as an attributed result rather than claiming to certify it
merely from this correction.

**Domains and orientations.** Rational identities alone do not justify
identifying real integration chains. The separate note
[reciprocal-region-beukers.md](reciprocal-region-beukers.md) supplies complete
inverse maps, positive orientation, nonnegative exhaustion and convergence.
It proves the exact region and cube interpretation requested by the parent.
The real chain is not the constant-term torus, and no claim equates those
cycles. The source's negative-octant sign is best checked against its
explicit positive cube integral, as the dictionary does.

**Quantum matrix.** The row-system convention is essential and is used
consistently. Eliminating ξ1,ξ2,ξ3 gives
\(D^4-z(2D+1)(17D^2+17D+5)+z^2(D+1)^2\).
This matches dividing the Apéry recurrence by n!, including the n²
coefficient of the preceding Borel coefficient. The scalar shift by cI
adds czI to G(z), hence multiplies the row solution by exp(cz).
Regularization translates Λ by c. Thus the unit-preserving M−5I and the
Apéry/Eisenstein M have genuinely different spectral normalizations; the
note correctly distinguishes them. The characteristic polynomial and
nonnegative irreducibility justify the Perron value 17+12√2 for M.

**Relative period and forcing.** GKS §5.4 displays the identical positive
cube integral, its first two values, and LV=−12t. Theorem 5.1 applies with
d=2,h=1 and the stated nonsingularity at zero, so the source's conclusion
that LV is a multiple of t has the cited hypotheses. Its coefficient is
I1−5I0=−12. The dictionary's B1=1 normalization is distinct from the
common B1=6 convention; the equations V=2ζ(3)A−12B and V(0)/12=ζ(3)/6
are correct in the convention actually declared.

**Chern/Gamma coefficients.** From H³=12, Hc2=24 and c3=−10 one obtains
c2=2H²,c3=−5H³/6. Newton sums are −3H² and −15H³/2.
Consequently the cubic log-Gamma coefficient is +5ζ(3)H³/2,
and exponentiation gives exactly (10), including the mixed γζ(2) term.
Dividing by Γ(1+H) removes γ and gives
\(1-2ζ(2)H²+(17/6)ζ(3)H³\). These are the source's coefficients.
The Laplace identity is correctly stated coefficientwise for the
Frobenius series including t^ε; its factor is Γ(1+ε), not its inverse.
The quoted deformed-sequence limit is only needed through cubic order,
as the dictionary specifies. It is not the rational companion limit B/A.

## Terminology refinement

Example 6.5 literally places the limiting cycle in
\(W_{-2}H_M^3(\widetilde X_0,\mathbb Q(3))\), whose bottom-piece description
leads to \(\mathrm{CH}^3(\operatorname{Spec}\mathbb Q,5)\) and its Borel
regulator. The phrase “weight-three piece” in the initial dictionary draft
could confuse motivic cohomological degree 3 with the displayed filtration
index −2. I recommended the exact notation, or “the motivic piece whose
regulator is Borel's ζ(3) regulator.” This does not change the ζ(3)Q
conclusion. The recommendation was sent to the author rather than editing
the actively maintained source file concurrently.

## What the new domain calculation adds

For R={(a,b,c)∈(0,1)³:(1−a)(1−b)<abc}, one has
\[
 I_n=\int_R\Lambda(-a,-b,-c)^n\frac{da\,db\,dc}{abc}.
\]
Writing D=1−w+uvw and (a,b,c)=(w,(1−w)/D,u) proves directly
that this is the Beukers cube integral. The positive minimum of Λ becomes
the positive maximum of Φ=Λ(−a,−b,−c), and
\(0<I_n\le2ζ(3)(17-12\sqrt2)^n\).
This is a proved rational-domain explanation of the reciprocal exponent,
not a claim of a new period, a new irrationality result, or historical
novelty. Controlled integer denominators remain a separate ingredient.
