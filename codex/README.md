# Mathathon independent research

Work in progress on `codex/independent-research`. Read [CHECKPOINT.md](CHECKPOINT.md)
for the latest state and [HANDOFF.md](HANDOFF.md) for remaining work. Results derived
in this session are not claims of historical novelty.

## Three reading routes

**Shortest complete proof, then recurrence motivation.** Start with the
[short elementary hybrid proof](problem131/short-elementary-hybrid-proof.md).
Rational poles supply the denominators; two integrations by parts supply a
positive integral; AM–GM and an elementary lcm estimate close the contradiction.
It requires no recurrence, PNT, or complex contours. The
[discovery guide](problem131/discovery-guide.md) then derives the square from
its constraints and explains the recurrence
through a quadratic Newton grid and finite balance, and gives a complete elementary
irrationality proof. Its main transfer now uses
[two reviewed mixed-product telescopes](problem131/apery-mixed-numerator-independent-review.md):
one gives the two-index recurrence, the other its adjacent balance, and their
combination explains \(17=6\cdot3-1\), \(5=6\cdot1-1\). The source at index zero
is explicit; the older diagonal certificate remains an appendix.
A single positive Laurent polynomial connects its
integer coefficients to the reciprocal-power Beukers integral; one AM–GM
minimum explains their reciprocal scales. The resistor chain gives a further
exact interpretation of the error. For the original paper's sequence of steps, use
[Zudilin explained](problem131/zudilin-explained.md). The
[40-hour research map](RESEARCH_MAP.md) selects a manageable 80/20 project
and offers a choice of one extension.

The fully real
[positive resolvent–Legendre proof](problem131/resolvent-legendre-pade-proof.md)
is a complete alternative, now accepted in its
[full-proof audit](problem131/resolvent-legendre-pade-independent-review.md).
A [separately audited Chebyshev-product argument](problem131/half-shift-elementary-independent-review.md)
now supplies its integral identity using real calculus alone.
Its opening design narrative starts from a positive Stieltjes error, finds
rational Bernoulli moments, transports Legendre polynomials through an exact
Hilbert Gram isometry, and obtains \(A_n\) by evaluation on the negative
half-square grid. Direct arithmetic and Delannoy growth finish the proof.
The underlying Wilson–Padé irrationality argument is explicitly attributed
to Prévost (1996); the reconstruction supplies elementary versions of its
norm, denominator, and growth steps. Use the mixed telescopes above for the
official recurrence-motivation question, which this alternative proof does
not need to solve.

**Follow the zeta-five simplification without the prime number theorem.**
Start with the [PNT-free guide](problem131/pnt-free-zeta5-guide.md).
It connects residue-class covariance to Bernoulli polynomials and superfactorials,
constructs an explicit normalizer, and completes the rationality contradiction.
A 25-term rational certificate controls every prime individually; both the
small-prime and large-prime repairs are included. The arithmetic rate is below
1.379, while the reviewed equilibrium estimate supplies decay faster than 1.38
in the quadratic exponent. A reviewed primitive-part argument keeps the final
degree at \(37n\) when \(K=40n\); the certificate denominator \(40000\) is only
auxiliary bookkeeping. Explicit floor rounding gives the same degree without
any added upper arithmetic cost. This computer-assisted alternative retains Fauzan's
determinant and local arithmetic construction, including its refined rank bound.
For fuller analytic motivation, read the [equilibrium guide](problem131/fauzan-equilibrium-guide.md).
The earlier [plain-rank route](problem131/fauzan-equilibrium-guide.md#12-spend-a-little-analytic-margin-to-simplify-the-local-arithmetic-too)
remains available when simpler local arithmetic is preferred; it uses the PNT.
The [optional outer-kernel improvement](problem131/outer-gain-budget-update.md)
adds a local lemma for a larger margin.
An [independently audited eleven-factor variant](problem131/eleven-factor-normalizer.md)
uses that extra lemma pointwise, keeping the PNT-free proof. It has fewer factors
but less pleasant coefficients and needs the stronger analytic certificate;
the 25-term route remains the primary exposition.

**Choose an extension: positive approximation or cyclotomic arithmetic.**
The [rank-one construction](problem131/rank-one-zeta-subspace-approximants.md)
gives strictly monotone rational bounds for every fixed odd zeta value.
Its optimized errors share their lower and upper exponential rates across
fixed odd indices, and an
[independently reviewed elementary trial](problem131/rank-one-elementary-trial-independent-review.md)
proves \(\limsup E_{h,s}^{1/h}<25/12544<1/500\), without assuming an exact
optimized rate. Denominator control remains separate.
A second [independently reviewed construction](problem131/higher-order-positive-stieltjes-independent-review.md)
uses a positive kernel of order m for zeta(2m+1), with explicit densities and
rational moments. It produces convergent rational lower bounds to each pure odd
zeta value, with no lower-zeta mixture. Its
[theorem and proof](problem131/higher-order-positive-stieltjes-zeta.md) explain why
ordinary kernels lose positivity and how increasing the kernel order restores it.
The [first twenty denominator experiments](problem131/pure-higher-kernel-pade-experiment.md)
show serious arithmetic costs; analytic convergence is not an irrationality proof.
The stronger [nested interval theorem](problem131/nested-higher-zeta-gauss-radau.md)
has also passed its [independent review](problem131/nested-higher-zeta-radau-independent-review.md):
explicit rational lower bounds rise and upper bounds fall to every odd zeta value.
The [real-order extension](problem131/real-order-positive-stieltjes-zeta.md)
connects these measures through beta products and a classical Brownian-exit law.
For an arithmetic alternative, use the
[cyclotomic depletion connection](problem131/cyclotomic-depletion-connection.md):
positive density combinations acquire integral polynomial moments, at the
cost of slower decay. The
[research results](problem131/research-results.md)
collect precise theorem statements, proofs, obstructions, and links. The common
analytic budget works for every fixed odd index, but the arithmetic normalizations
needed for all-odd irrationality remain unproved. The negative results identify
specific constructions that fail; they are not global impossibility claims.

A concrete new connection for a short presentation is the
[log2 endpoint proof](problem131/log2-endpoint-catalan-arithmetic.md), now
[independently reviewed](problem131/log2-catalan-independent-review.md).
The same Wilson–Legendre family joins logarithms to zeta3; its logarithmic
numerator factors into Catalan numbers and beta integrals, so a single lcm
clears it. This supplies an explicit arithmetic link to the Catalan project.

## Secondary project: exact prime-power Fuss–Catalan valuations

The [problem #137 discovery guide](problem137/discovery-guide.md) gives a standalone
combinatorial proof of the exact valuation and its leading unit modulo the prime.
Expand each vertex into a perfect prime-ary block, classify the smallest group
orbits, and count their depth-constrained skeletons by an explicit cycle-splicing
bijection with inverse. The proof uses no factorial-valuation formula or Lagrange
inversion. It addresses the exactness gap named in the [official statement](sources/problem137-official.md);
[the primary-source comparison](problem137/konvalinka-comparison.md) credits the
binary Deutsch–Sagan argument and distinguishes Konvalinka’s existing arithmetic
formula and tree-action divisibility proof. Historical novelty and weighted
extensions remain unresolved.

## Connections worth following

| Question | Read |
|---|---|
| Why precisely this square and this recurrence? | [Design classification](problem131/zudilin-design-classification.md); [two mixed telescopes and their adjacent balance](problem131/apery-mixed-numerator-balance.md) |
| Why does a growing integer solution create a small positive error? | [The same Laurent polynomial in positive and inverse powers](problem131/apery-laurent-reciprocal.md), with [independent review](problem131/apery-laurent-independent-review.md); [resistor chain](problem131/discrete-green-motivation.md) |
| Can advanced machinery make the structure clearer? | [Elliptic periods and modular branch cancellation](problem131/apery-modular-elliptic.md); [completed independent review](problem131/apery-modular-review.md) |
| How can the same moments control real size and prime denominators? | [Bernoulli local/global dictionary](problem131/bernoulli-local-global-dictionary.md); [local arithmetic reconstruction](problem131/fauzan-local-arithmetic-review.md) |
| Why do prime-allocation floors suggest a factorial product? | [Centered indicators, Brownian-bridge covariance, and superfactorials](problem131/floor-covariance-barnes-design.md); [PNT-free assembly](problem131/pnt-free-zeta5-guide.md) |
| What do positive matrices reveal—and where do they fail? | [All-odd monotone rank-one bounds](problem131/rank-one-zeta-independent-review.md); [moving-pole Wilson–Apéry identity](problem131/rank-one-apery-bridge-independent-review.md) |
| What connects this to prime-adic Fuss–Catalan questions? | [Carry cancellation and digit factorization](problem131/apery-digit-arithmetic.md); [elementary character rule](problem131/apery-cm-character.md), now [independently reviewed](problem131/apery-cm-independent-review.md) |
| What transfers to other irrationality proofs? | [Shared exponential, pi-squared, and zeta coefficients](problem131/shared-coefficients-exp-pi-zeta.md); [zeta-two/zeta-three balance](problem131/zeta2-zeta3-balance-unification.md); [powers of pi](problem131/pi-powers.md) |

The new [monomial-kernel design](problem131/monomial-kernel-design.md) makes the
common arithmetic pairing explicit: diagonal monomials select the zeta constant,
while off-diagonal terms are rational. Its [dedicated review](problem131/monomial-kernel-review.md)
is complete; choosing partners with sufficient nonzero decay remains a separate
problem. Follow-up periodic-weight and mixed-partner constructions have separate
review status.

For the geometric interpretation, [the Legendre/relative-period note](problem131/laurent-legendre-relative-periods.md)
records the necessary boundary term: the real error generating function satisfies
an inhomogeneous equation. The Laurent critical values alone are not a derivation
of the full recurrence. The [reviewed geometric divergence proof](problem131/geometric-divergence-proof.md)
now supplies the full exact-form certificate and computes that boundary directly;
its large rational primitives remain a presentation challenge.

The [reciprocal-region bridge](problem131/reciprocal-region-beukers.md)
offers a shorter geometric picture: an explicit birational involution sends
the positive octant to a bounded Beukers region, and a rational cube
parametrization recovers the integral and its exact Jacobian. The Laurent
minimum becomes the Beukers-kernel maximum. The core Laurent polynomial,
Beukers relative period, and inhomogeneous period equation are **published**
in Golyshev–Kerr–Sasaki, *Apéry extensions*, §5.4; the
[Fano dictionary](problem131/fano-v12-apery-gamma-dictionary.md) records the
exact match and correction of a printed birational-map coordinate. Their
reconstruction here is not a new mirror-model construction.

For higher odd indices, the
[reviewed Euler-operator hierarchy](problem131/odd-hierarchy-independent-review.md)
expresses all polynomial moment functionals through \(\mu_3\), with fixed-order
operators that preserve integer coefficients after the stated scalar clearing.
This does not preserve pointwise positivity or bound the denominators of a
solved moment system. The new
[balanced hyperbolic normality proof](problem131/hyperbolic-balanced-normality.md)
has now passed its [full independent audit](problem131/hyperbolic-balanced-normality-independent-review.md):
every rectangular balanced index has a unique monic polynomial with simple
positive roots. This gives no denominator-growth or error-rate theorem.
Its Bernoulli construction's factorial constant-term formula remains finite
evidence, not a theorem.

One nearby distinction is already settled: the moving-pole Wilson bounds equal
the adjacent-Apéry quotients for every index and inherit Apéry's arithmetic.
Equality with the original \(D_{2h}\) rank-one determinant roots is still an
unproved norm identity, despite finite exact agreement.

Detailed higher-derivative, truncation, and prime-band investigations are indexed
in [research results](problem131/research-results.md); their distinctions between
actual failure and insufficient estimates matter.

The [primary-source comparison](problem131/rodrigues-primary-source-comparison.md)
identifies the higher-Rodrigues pure forms, finite Barnes decomposition, and
positive central period as exact specializations of Nesterenko's 2003 work.
The underlying polynomial is also a classical confluent Jacobi–Piñeiro
specialization. Our explicit Laurent geometry and certified asymptotic analyses
are developments within that established framework; their historical priority
has not been established. The
[completed finite-saddle review](problem131/finite-rodrigues-saddles-independent-review.md)
now covers the oscillatory analytic forms for \(\zeta(5),\zeta(7),\zeta(9),\zeta(11)\).
Their available arithmetic estimates still do not yield irrationality proofs.

## What has been verified

The elementary zeta-three guide is a complete written proof. The advanced modular
route now includes a reviewed direct global covering argument, with standard
modular-form and covering-theory inputs stated explicitly.
The resolvent–Legendre alternative also has a completed full-proof review;
its positive error, cubic-lcm numerator arithmetic, and growth estimates are
all proved, with the short Fourier contour calculation stated as an input
derived in the text.

For the zeta-five alternative, the [completed assembly audit](problem131/fauzan-simplified-proof-audit.md)
checks the shared positive moments, local arithmetic, real determinant transfer,
and final contradiction. The new [normalizer review](problem131/pnt-free-normalizer-independent-review.md)
checks the superfactorial replacement for prime averaging, including every higher
prime power and eventual coefficientwise integrality. Its exact sparse certificate
proves the arithmetic rate is below 6893/5000; the [rational energy review](problem131/rational-energy-review.md)
certifies the analytic magnitude exceeds 2771/2000. The
[primitive-part and rounding review](problem131/prime-exponent-rounding-independent-review.md)
removes the unnecessary determinant power. The
[new guide](problem131/pnt-free-zeta5-guide.md) assembles these dependencies.
The older plain-rank route and its exp(−13n²) bound remain valid with their own
normalizer and PNT input; the outer-kernel option gives exp(−26n²).
These are computer-assisted derivations retaining Fauzan's construction, not
claims of elementary proofs, historical novelty, formal verification, or external
acceptance.
No all-odd irrationality theorem is claimed. See [source provenance](sources/odd-zeta-provenance.md).

Separately, the pinned Romik Lean formalization **was built locally**, and its
ordinary-series statement and axiom checks passed. The completed kernel replay
checked all 71 project modules against their imported environments. The compiled
theorem retains the explicit prime number theorem axiom, along with Lean's
standard logical axioms. This was not a fresh replay of all Mathlib constants
and does not formalize our new equilibrium, PNT-free normalizer, or Python
certificates. The [verification report](problem131/lean-local-verification.md)
links the exact commit, successful build log, axiom output, and durable replay
result.

Symbolic scripts check identities; finite experiments alone do not prove all-index
statements. Reconstructed motivation explains how one might design a proof, not
what its historical authors thought. All work is under `codex/`; no Claude-owned
files have been changed. The exact [#131](sources/problem131-official.md) and
[#137](sources/problem137-official.md) statements were retrieved from the public
API at 08:46 UTC on 14 November 2026; raw JSON and SHA256 provenance are saved.
Earlier HTML snapshots are application shells, not statements. The prime-power
tree theorem matches #137's arity scope; the primary-source comparison with
Konvalinka is now recorded in [the comparison note](problem137/konvalinka-comparison.md).

## Reproduce

Most scripts use Python 3 with SymPy and mpmath; their headers state dependencies.
The [Newton closure check](computations/check_newton_diagonal_closure.py) and
[modular review check](computations/check_apery_modular_review.py) use exact identities.
The [mixed-product check](computations/check_apery_mixed_numerator_balance.py)
verifies the two telescopes symbolically; its
[independent product check](computations/check_apery_mixed_independent.py)
also tests every small boundary case. The
[reciprocal-region check](computations/check_reciprocal_region.py) verifies both
inverse maps and their Jacobians. The
[elementary rank-one trial certificate](computations/check_rank_one_elementary_trial.py)
uses exact rational bounds, with no floating-point acceptance test.
The [independent resolvent–Padé checker](computations/check_resolvent_pade_independent.py)
uses standard-library rational arithmetic to check the polynomial isometry,
associated numerator, and even-integer denominator lemma.
The [independent Laurent check](computations/check_apery_laurent_independent_review.py)
verifies the reciprocal measure and covariance separately from the main derivation.
The sparse superfactorial certificate, its independent normalizer checks, and
the endpoint and rational-energy certificates use only the Python standard
library and local interval primitives. From the repository root:

```bash
python codex/computations/certify_equilibrium_endpoints.py
python codex/computations/check_rational_energy_primitives.py
python codex/computations/certify_equilibrium_energy_rational.py
python codex/computations/check_floor_covariance_barnes.py
python codex/computations/check_sparse_barnes_certificate.py
python codex/computations/check_pnt_free_normalizer_review.py
python codex/computations/check_prime_exponent_rounding_review.py
```

The sparse check verifies 561 exact breakpoints and propagates the inequality
to the infinite tail by an exact period identity. The independent normalizer
check tests the factorial identities and repair bounds. For the older PNT route,
run [check_table_free_arithmetic.py](computations/check_table_free_arithmetic.py).
The [rational energy review](problem131/rational-energy-review.md) explains every
remainder and outward-rounding bound; the saved JSON records exact dyadic endpoints.
Older [Carlson](computations/certify_equilibrium_energy.py) and
[virial](computations/certify_equilibrium_energy_elementary.py) checks use
`python-flint`/Arb and provide additional cross-checks. The rational certificate
shares the endpoint box and proved potential identities; it is not an independent
reconstruction of those inputs.

For the secondary project, [the cycle-splicing check](problem137/check_cycle_splicing.py)
exhaustively tests the bijection and its inverse; [the independent orbit check](problem137/check_prime_power_orbit_review.py)
counts skeleton depth profiles separately from the closed formula.
