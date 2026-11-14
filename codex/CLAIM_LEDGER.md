# Claim ledger

All 'proved' labels refer to arguments written and checked in this workspace,
not historical novelty. Standard external inputs are identified in each note.
Finite computation is supplementary evidence only.

## Current status — 14 November 2026

This section supersedes pending-review or verification-status wording in the
earlier stage snapshots below. Those snapshots preserve the development history;
they are not a list of currently running tasks. No historical novelty or
all-odd irrationality theorem is claimed.

| Current result | Proof / completed review | Exact scope |
|---|---|---|
| PNT-free zeta-five normalization keeps degree \(h=37K/40\) | [Guide](problem131/pnt-free-zeta5-guide.md); [normalizer review](problem131/pnt-free-normalizer-independent-review.md); [floor-rounding and primitive-part review](problem131/prime-exponent-rounding-independent-review.md) | Rational majorant denominator \(40000\) is auxiliary only. Primitive parts remove the final determinant power; explicit floor rounding adds no upper cost. Both prime repairs and source local lemmas remain necessary. |
| Optional eleven-factor PNT-free normalizer, cost \(<6923/5000\) | [Independent audit and complete table](problem131/eleven-factor-normalizer.md); [independent exact checker](computations/check_eleven_barnes_independent.py) | Adds the reviewed outer Hermite-kernel lemma; requires analytic magnitude \(2771/2000\), giving margin \(9/10000\). \(C=113,M=10^7\); denominator \(3545916\) is auxiliary. Main 25-term route retained. |
| Pinned existing Lean project builds and all 71 project modules replay successfully | [Verification report](problem131/lean-local-verification.md); [compiled axiom check](verification/romik-lean/independent-check.log); [durable replay result](verification/romik-lean/kernel-replay-result.json) | Explicit PNT axiom, plus standard Lean axioms. Replay is against imported environments, not a fresh replay of all Mathlib declarations. This verifies the pinned source formalization, not our new PNT-free proof or Python certificates. |
| Adaptive adjacent-Legendre boundary repair has sharp positive asymptotic \(\pi n^{-3}\log n\,16^{-n}\) | [Construction](problem131/adaptive-boundary-cancellation.md); [completed independent review](problem131/adaptive-boundary-independent-review.md) | Eventually defined rational tuning cancels the total boundary sum. Sufficient denominators are proved; finite primitive-arithmetic samples do not prove an asymptotic obstruction. |
| Signed single-\(\zeta(7)\) model, finite-window determinant nonvanishing, and total-variation obstruction | [Construction](problem131/signed-single-zeta7-model.md); [independent review](problem131/signed-zeta7-independent-review.md) | Applies the primary general nonvanishing proposition after checking its hypotheses, not its large-\(s\) theorem at \(s=7\). The positive lower bound concerns the specified total-variation majorant and multiplier, not the signed determinant or all primitive normalizations. |
| Euler-depleted positive functional has sharp polynomial denominator \(2^{1+v_2(s-1)}\) | [Construction](problem131/euler-depleted-positive-moments.md); [completed review](problem131/euler-depletion-independent-review.md) | Odd \(s\ge3\); pole identities, local divided-difference losses and determinant scaling reviewed. The polynomial distribution exponent is still \(p^{1-s}\); no improved determinant normalizer follows. Numerical determinant tables remain diagnostics. |
| Integer depletion has least uniform polynomial denominator \(\prod_{p\mid\ell}p^{1+v_p(s-1)}\) | [Cumulant construction](problem131/integer-depletion-cumulants.md); [completed independent review](problem131/integer-depletion-independent-review.md) | Odd \(s\ge3\), integer \(\ell\ge2\). Discrete-uniform cumulants explain prime support; Lucas witnesses prove sharpness. Restoring the real decay by the stated rescaling creates actual growing moment denominators along explicit subsequences, not a general determinant obstruction. |
| Higher-Rodrigues \(\zeta(5)\) contour has sharp rate \(\lambda^{-1/2}\) and infinitely many terms of each sign | [Laurent construction](problem131/higher-rodrigues-laurent-geometry.md); [completed \(m=2\) review](problem131/higher-rodrigues-asymptotic-independent-review.md) | Saddle selection and oscillating asymptotic proved; quadratures only corroborate. The available sufficient denominator estimate does not close irrationality and is not an actual primitive-denominator lower bound. |
| All-\(m\) Eulerian contour decomposition and central-component theorem | [Construction](problem131/eulerian-rodrigues-period-decomposition.md); [completed independent core review](problem131/eulerian-periods-independent-review.md) | Each fixed \(m\): exact decomposition, conditional component saddle theorem, unconditional central asymptotic and positive central Mellin convolution. No uniform estimate in growing \(m\) or full-form all-odd decay theorem. |
| Certified full-form oscillatory asymptotics for \(m=2,3,4,5\) | [Completed finite-saddle review](problem131/finite-rodrigues-saddles-independent-review.md) plus core component theorem | Pure forms for \(\zeta(5),\zeta(7),\zeta(9),\zeta(11)\); unique dominant conjugate pair, exponential limsup below one, infinitely many terms of each sign. No irrationality conclusion, all-index nonvanishing, or \(m=6\) extension. |

The [primary-source comparison](problem131/rodrigues-primary-source-comparison.md)
now identifies the rational Rodrigues family, pure derivative forms, finite
Barnes decomposition, and positive central period as exact specializations of
Nesterenko's 2003 results. The polynomial is also a classical confluent
Jacobi–Piñeiro specialization. These constructions are not historically new.
For the higher-Rodrigues family, the specific Laurent geometry and certified
asymptotic constants remain historically unattributed, not established as
novel. The base Apéry three-variable Laurent polynomial and its reciprocal
Beukers period are an exact published V12 mirror model: see the
[primary-source dictionary](problem131/fano-v12-apery-gamma-dictionary.md)
to Golyshev–Kerr–Sasaki, *Apéry extensions*, §5.4.

## Earlier established results and stage snapshots

| Claim | Proof | Independent review / checks | Limits |
|---|---|---|---|
| Zudilin's seven lemmas and arithmetic | zudilin-explained.md | independent-review.md; exact generic certificates | Reconstruction of published mathematics |
| Short recurrence certificate forced by pole matching | simplification.md | generic symbolic rank and identity | Uniqueness only in stated ansatz |
| Complete recurrence-only irrationality route | apery-asymptotics.md | parent checked arithmetic, Casoratian and growth; check_apery_asymptotics.py | Classical ingredients reorganized; no novelty claim |
| Beukers integral equals Zudilin form | beukers-explained.md | independent Mellin derivations; exact kernel checks | Primary original Russian translation and later author exposition retrieved |
| Four-factorial elementary lcm rate | denominator-design.md | exact period-six table and integer rate comparisons | Sufficient asymptotic bound, not proof D_n<3^n |
| Apéry generating function is square of second-order ODE solution | symmetric-square.md | generic rational identities and formal uniqueness | No modular parametrization or new irrationality follows by itself |
| All higher odd derivatives have nonzero boundary forcing | forced-recurrence-obstruction.md | independent log coefficient and sign audits | Only the unmodified derivative family |
| Full derivative forms decay polynomially, with exact leading term | log-kernel-extension.md | odd-zeta-independent-review.md | Fixed derivative order as n grows |
| No fixed polynomial in H_n over Q(n) repairs exponential decay | forcing-correction.md | pole-propagation proof and exact residual checks | Growing-length corrections excluded from theorem |
| Truncated third derivative has rate16^-n and constant pi/2 | truncated-derivative.md | independently reproduced by two agents; truncation-review.md | Eventual positivity, not yet all-n positivity |
| Truncation introduces primes above n with valuation−3 | truncated-derivative.md | two independent arithmetic reviews and finite checks | Specified n=(p+1)/2, p!=17; not a full denominator-growth theorem |
| All odd zeta values admit positive single-constant moment model | odd-moment-generalization.md | independent-review.md | Hermite formula input; arithmetic remains unresolved |
| Compact-support moment criterion for irrationality | determinant-capacity-recipe.md | odd-zeta-independent-review.md | Does not apply unchanged to n-dependent unbounded weights |
| Bernoulli large-prime cancellation rule | odd-arithmetic-barrier.md | exact local tests and matrix normalization checks | von Staudt–Clausen input |
| Positive weighted rational family can cancel zeta3 | orthogonality-design.md | corrected residue formula independently checked | Tuning positivity checked finitely, asymptotic arithmetic unresolved |
| Explicit classical proofs for pi powers and exp(r) | pi-powers.md | pi-powers-review.md | Exposes transcendence machinery, not low-prerequisite new proof |

## Source-verification status

Fauzan's 2026 preprint states irrationality of zeta5. The pinned Romik
formalization has now passed a local 8995-job build, two statement checks,
the compiled axiom report, and replay of all 71 project modules. Its explicit
mathematical axiom is the PNT. Exact versions, logs and replay limits are in
[the verification report](problem131/lean-local-verification.md) and
[source provenance](sources/odd-zeta-provenance.md). This separate formal check
does not formalize the new constructions in this research branch.

## Fifth-stage independently checked components

| Claim | Proof | Review | Limit |
|---|---|---|---|
| Positive one-interval equilibrium measure and global potential maximum | fauzan-analytic-simplification.md; general-equilibrium-family.md | equilibrium-independent-review.md | Standard logarithmic-potential framework; no novelty claim |
| Endpoint existence and energy < -69/50 | equilibrium-endpoint-certificate.md; fauzan-equilibrium-certificate.md | endpoint-independent-review.md; equilibrium-energy-review.md | Two energy reductions share Arb; not formal verification |
| Same leading analytic bound for every fixed odd index | common-odd-zeta-analytic-budget.md | common-odd-zeta-budget-review.md | Arithmetic normalizer unproved for other indices |
| Functional odd-index shift and simple-pole obstruction | odd-functional-shift.md | exact symbolic and prime-band script | Uniform pole normalization cost is not determinant cost |
| Inhomogeneous recurrence for W_n | W-recurrence.md | exact recurrence checks and independent prime audit | Does not imply density of nonvanishing modulo varying primes |

## Sixth-stage developments

| Claim | Proof | Review | Limit |
|---|---|---|---|
| Zudilin saturated interpolation and optimal difference lattice | zudilin-design-classification.md | zudilin-design-independent-review.md | Specified pole/zero/degree constraints |
| Dilation breaks exponential decay, exact asymptotic | zudilin-design-classification.md | same review, including second term | Fixed integer dilation d>=2 |
| Sharp discrete inequality and all Green entries | discrete-green-motivation.md | discrete-green-independent-review.md | Identification with zeta3 uses separate bounded form |
| Table-free and plain-rank arithmetic budget | fauzan-arithmetic-simplification.md | fauzan-arithmetic-simplification-review.md | New relaxed normalizer; full assembly audit pending |
| Local determinant arithmetic including small primes | fauzan-local-arithmetic-review.md | independent reconstruction and3335 outer cases | Full proof assembly separately reviewed |
| Bernoulli boundary-rule and general odd distribution dictionary | bernoulli-local-global-dictionary.md | bernoulli-dictionary-review.md | Does not imply determinant valuation bound |
| All-real roots and exact inertia of zeta pencils | real-rooted-zeta-pencils.md | editorial proof review; deeper localizing review active | No denominator control or convergence rate |
| Fixed odd-index arithmetic envelope obstruction | odd-inner-prime-shift.md | exact checks; proof review pending | Obstruction to sufficient estimate, not actual normalizers |
| Fuss-Catalan forced carries and finite digit transducer | problem137/valuation.md |36120 exact cases; proof review pending | Official problem statement unverified |

## Seventh-stage proof status

| Claim | Proof | Review | Limit |
|---|---|---|---|
| Full simplified computer-assisted zeta5 derivation | fauzan-simplified-proof-audit.md | end-to-end review and all certificates rerun | Retains Fauzan construction; not novelty or external consensus |
| Exact rational-only energy upper bound -2771/2000 | rational-energy-certificate.md | rational-energy-review.md; independent primitive/rerun checks | Shares proved equilibrium formulas/endpoint box |
| Direct derivation of17and5 from kernel closure | newton-bispectral-motivation.md | newton-symmetric-square-review.md | Identification/denominators remain separate proof steps |
| Fully motivated modular/elliptic zeta3 route | apery-modular-elliptic.md | apery-modular-review.md | Standard modular-curve facts explicit |
| Fixed-pole optimized rational bounds and actual denominator failure | one-pole-schur-structure.md | one-pole-independent-review.md | Does not rule out all fixed-denominator constructions |
| General Lucas factorization and square/Delta-square digit polynomial | apery-digit-arithmetic.md | apery-digit-square-review.md | Character choice final review ongoing |
| Prime-power arity group valuation and leading unit | problem137/prime-power-sylow-audit.md | problem137/prime-power-orbit-review.md | Uses multivariable Lagrange inversion; arbitraryarityoutside theorem |


## Eighth-checkpoint snapshot and corrections

This historical addendum superseded earlier access and pending-review labels;
the current summary at the top of this ledger now takes precedence.
Earlier entries are retained as a record of what had been checked at that time,
not as the current verdict. All proof paths below are relative to `codex/`.

- Exact official statements are now saved in
  [problem131-official.md](sources/problem131-official.md) and
  [problem137-official.md](sources/problem137-official.md), with raw JSON and
  provenance. Application-shell downloads are not statements.
- The local arithmetic and full simplified zeta-five assembly are no longer
  pending. [The end-to-end audit](problem131/fauzan-simplified-proof-audit.md)
  found no unresolved bridge in the alternative derivation. It retains
  Fauzan's determinant/local construction and named standard inputs; this is
  not external acceptance, historical novelty, or formal verification of our
  alternative proof. The separate source formalization's local verification
  is recorded in the current summary above.
- The plain-rank route uses a **new normalizer** \(\widetilde m\), with outer
  lower exponent \(v_p(S_K)+2\sum w_i-r_p\). The earlier source normalizer with
  the refined zero-weight correction is not interchangeable with it.
- The character rule is no longer an experimental Legendre-symbol proposal.
  Its elementary differential-involution proof has completed independent review.
  The older limited CM review is not the current boundary of verification.
- The #137 counting proof no longer needs Lagrange–Good inversion. An explicit
  cycle-splicing bijection and inverse replace that step. Algebraic counting
  audits remain independent checks. The exact official problem matches the
  prime-power theorem; arbitrary arity and weighted extensions remain outside it.

| Current claim | Proof and review | Exact scope / remaining limit |
|---|---|---|
| Rational-only equilibrium certificate \(J<-2771/2000\), with exact outward enclosure \([-1.385510140104,-1.385509421992]\) | [Certificate](problem131/rational-energy-certificate.md); [independent primitive and assembly review](problem131/rational-energy-review.md) | Standard-library integer/dyadic arithmetic; shares the proved endpoint enclosure and potential formulas. Removes Arb from this scalar certificate, not advanced inputs from the whole proof. |
| Same simplified integer polynomials satisfy \(0<\widetilde Q_{40n}(\zeta(5))<e^{-13n^2}\) eventually | [Complete assembly](problem131/fauzan-simplified-proof-audit.md) plus [stronger certificate](problem131/rational-energy-certificate.md) | Arithmetic budget \(2203747/1600000\); exact margin \(13053/1600000\), giving rate \(13.053n^2\) before lower-order errors. No explicit starting index, new irrationality measure, or all-odd theorem. |
| The same finite coefficients give endpoint-polynomial proofs for \(e^r\), \(\pi^2\), and a tangent construction | [Shared coefficients](problem131/shared-coefficients-exp-pi-zeta.md); [independent review](problem131/shared-coefficients-review.md); [tangent consequence](problem131/shared-coefficients-tangent.md) | \(r\ne0\) rational for exponential/tangent claims. Tangent nonvanishing uses the Casoratian, not a false positivity assertion. The elementary \(\pi^2\) argument alone does not prove all nonzero rational powers of \(\pi\) irrational. |
| Mixed Newton/binomial kernel gives the zeta-two recurrence with coefficients 11 and 3; common finite balance gives the zeta-three closure | [Balance construction](problem131/zeta2-zeta3-balance-unification.md); [independent review](problem131/zeta-balance-independent-review.md); [Newton closure review](problem131/newton-symmetric-square-review.md) | The degree-at-most-two classification is exactly for the stated single balance polynomial and yields pairs \((1,1),(2,1),(2,2)\); it does not classify every possible recurrence or auxiliary family. |
| For every fixed positive \(r,s\), binomial-weight modes are \(Nx_*+O(1)\) and logarithmic growth is \(N\log[(1+x_*)^s/(1-x_*)^r]+O(\log N)\), where \((1-x_*)^r(1+x_*)^s=x_*^{r+s}\) | [General theorem and proof](problem131/zeta-balance-independent-review.md), Section 3 | No restriction \(s\le r\) is needed. This concerns denominator-sequence growth, not identification or decay of a companion linear form. |
| For \(s\ge2\), the monomial logarithmic pairing of integer polynomials isolates \((s-1)(\sum p_kq_k)\zeta(s)\) with rational remainder cleared by \(D_n^s\) | [Monomial-kernel theorem](problem131/monomial-kernel-design.md) | [Independent review complete](problem131/monomial-kernel-review.md), including finite-partial-sum cancellation for s=2. Arbitrary-partner positivity and sufficiently fast decay are not claimed. |
| For \(p\ge5\), \(P_p=\Delta^{(1-(\frac{-6}{p}))/2}Q_p^2\), where \(\Delta=1-34t+t^2\), \(Q_p(0)=1\), \(Q_p\) is squarefree and coprime to \(\Delta\) | [Character construction](problem131/apery-cm-character.md); [completed independent review](problem131/apery-cm-independent-review.md); [root-multiplicity proof](problem131/apery-digit-square-review.md) | Main proof uses finite sums, the exact differential pullback, and a deck involution; CM is optional interpretation with standard geometric inputs. No historical novelty assertion. |
| Exact prime-power Fuss–Catalan valuation and leading unit by finite group action and a bijection | [Standalone guide](problem137/discovery-guide.md); [orbit review](problem137/prime-power-orbit-review.md); [cycle splicing](problem137/cyclic-type-cycle-splicing.md) | \(q=p^e\), \(L=(q-1)n+1=\sum d_jp^j\), \(D=(\sum d_j-1)/(p-1)\); \(v_p(T_q(n))=D\) and \(p^{-D}T_q(n)\equiv(-1)^D/\prod d_j!\pmod p\). No factorial-valuation theorem or inversion formula in the final proof. |
| The typed skeleton count is \(\prod_i(pk_i)!/[\prod_i k_i!\prod_jd_j!]\), including zero internal-type cases | [Explicit forward/inverse bijection](problem137/cyclic-type-cycle-splicing.md); [independent algebraic check](problem137/typed-count-matrix-tree-audit.md) | Fixed missing-parent name in arbitrary assignments corresponds to a variable root name in trees. Exact checks cover 24,242 assignments and surjectivity; a separate depth-profile computation covers 234 cases, 28 with zero internal types. These checks supplement the universal proof. |
| The #137 proof addresses the exactness step distinguished in Konvalinka's primary paper | [Primary-source comparison](problem137/konvalinka-comparison.md); [official statement](sources/problem137-official.md) | Proposition 2 gives arithmetic exactness; Proposition 7 gives tree-action divisibility. Deutsch–Sagan owns the established binary orbit argument. This comparison does not establish historical novelty of the prime-power refinement. |

### Questions recorded at the eighth checkpoint

- **Actual local nullvectors:** odd-zeta agent is deriving leading local-block
  nullvectors to test determinant-specific valuation gains. A failure of a
  sufficient entrywise budget is not a lower bound on the actual normalizer.
- **Half-reflected \(W_n\) sums:** connections agent is seeking useful positive
  representations. A signed identity or finite nonvanishing tests do not prove
  positivity or nonvanishing modulo a varying prime.
- **Monomial extensions:** the original pairing theorem's independent review is
  complete. Periodic-weight and mixed-partner constructions are being developed
  separately; their review status does not follow from the original theorem.

No proof of irrationality for all odd zeta values, weighted Fuss–Catalan
extension, historical novelty claim, or guaranteed usage/runtime budget follows
from this ledger.

## Ninth-stage independently checked discoveries

| Claim | Proof | Review | Scope |
|---|---|---|---|
| One Laurent polynomial gives Apéry constant terms and reciprocal Beukers periods | apery-laurent-reciprocal.md; fano-v12-apery-gamma-dictionary.md | apery-laurent-independent-review.md; exact primary-source match | Core model and relative period already published in Golyshev–Kerr–Sasaki §5.4; independently derived here, not new |
| Bernoulli numbers are logarithmic derivatives of a positive period tower | bernoulli-period-tower.md; laurent-deformation-arithmetic.md | bernoulli-period-independent-review.md | Higher odd constants appear as multiple integration constants |
| Rational Euler projection isolates the top constant and gives a positive Bose tail | bernoulli-period-projectors.md | exact operator checks and parent derivation | Preserves Taylor denominators, not rational-q evaluated arithmetic |
| Multiplicity lift restores exponentially small single-ζ5 forms | multiplicity-lifted-zeta-family.md | multiplicity-lifted-independent-review.md | Infinitely often nonzero, but actual new-prime denominators obstruct the old normalization |
| Actual ζ7 inner-kernel cancellation improves a sufficient budget by1/25200 | inner-leading-block-cancellation.md | inner-leading-block-independent-review.md; inner-leading-block-adversarial-review.md | Small gain, not an irrationality proof or exact normalizer asymptotic |
| Laurent one-variable periods are Legendre P/Q with a boundary correction | laurent-legendre-relative-periods.md | exact symbolic checks; parent inspection | Full reduction subsequently completed and reviewed in geometric-divergence-proof.md; real series satisfies L F=-12t |

The general Laurent-deformation geometry and higher-pole constructions retain
their individual review labels. The local Lean build was still active at this
stage; it and the project replay have since completed successfully. Neither
formalizes our simplified proof.

| Further checked claim | Proof / evidence | Limit |
|---|---|---|
| Actual ζ5 outer correction kernel gains K+3N-3p+1 | outer-polynomial-correction-kernel.md; outer-correction-adversarial-review.md | Specified outer band; original source dependencies remain |
| New normalizer budget2190247/1600000, eventual exp(-26n²) | outer-gain-budget-update.md and completed assembly dependencies | Stronger margin, optional extra lemma, no novelty claim |
| Higher positive Laurent amplitudes mix lower zeta values | positive-period-amplitude-mixing.md; positive-period-amplitude-review.md | Precisely the displayed family; not a general impossibility theorem |
| Existing pinned Lean theorem compiles with ordinary-series target | verification/romik-lean/build.log; independent-check.log; repository-check.log | Explicit PNT axiom; project replay subsequently completed for all 71 modules; not our new proof formalized |

## Tenth-stage proof status

| Claim | Proof / verification | Limit |
|---|---|---|
| PNT-free explicit all-prime normalizer | barnes-normalizer-without-pnt.md; pnt-free-normalizer-independent-review.md | Retains source local lemmas; eventual integrality, no practical first K |
| Sparse25-term profile and patched arithmetic exponent<1.3786 | check_sparse_barnes_certificate.py; sparse_barnes_verification.json | Exact finite computation, no LP solver in verification |
| Complete PNT-free ζ5 alternative route | pnt-free-zeta5-guide.md plus completed component reviews | Not formally verified; computer-assisted; not historical novelty |
| Explicit Hermite ideal is the entire selected outer kernel | outer-correction-Hermite-ideal.md; independent review | Optimal only within selected C row space, not exact determinant valuation |
| Odd-prime next normalized Fuss–Catalan digit | problem137/next-digit-orbits.md; next-digit-independent-review.md | p odd and arityp^e; no modulo-p³ or novelty claim |
| Adaptive adjacent-Legendre pureζ5 repair | adaptive-boundary-cancellation.md; adaptive-boundary-independent-review.md | Independent analytic review completed; arithmetic samples not a failure theorem |

The first optional replay had incomplete result capture and is not used as a
certificate. The second replay completed at 10:28:40 UTC on 14 November 2026,
exit status 0, covering exactly 71 project modules against imported environments.
Its durable record supports the current verification statement above; it is
not a fresh replay of all Mathlib constants.

## Fifteenth-stage checked constructions

All paths below are under problem131/. “Checked” records proof review, not
historical priority. Earlier stage-specific open labels are superseded here.

| Claim | Proof and independent review | Scope |
|---|---|---|
| Full Apéry recurrence from two short mixed telescopes | apery-mixed-numerator-balance.md; apery-mixed-numerator-independent-review.md | Includes n=0 forcing and associated numerator; no large diagonal certificate needed |
| Exact V12 mirror dictionary and corrected reciprocal involution | fano-v12-apery-gamma-dictionary.md; fano-apery-independent-review.md | Core mirror/period construction published; localized source-map correction verified |
| Bounded reciprocal region gives Beukers cube | reciprocal-region-beukers.md; check_reciprocal_region.py | Correct real domain and orientation; not a torus-cycle identification |
| Shifted-Jacobi mixed recurrence and positive Hurwitz3 forms | shifted-jacobi-mixed-apery-family.md; shifted-jacobi-independent-review.md | Real shifts a>=0; rational normalizers; no all-Hurwitz irrationality result |
| Rational monotone one-sided bounds for all odd zeta values | rank-one-zeta-subspace-approximants.md; rank-one-zeta-independent-review.md | Arithmetic control remains open |
| Universal optimized exponential rates and elementary bound below1/500 | rank-one-universal-error-rate.md; rank-one-elementary-trial-rate.md; their independent reviews | No adjacent-Apéry conjecture used; no exact optimized rate asserted |
| Moving-pole Wilson bounds equal adjacent-Apéry quotient | rank-one-apery-bridge-draft.md; rank-one-apery-bridge-independent-review.md | Original D_(2h) Schur equality still conjectural |
| Euler hierarchy with sharp normalizer and exact higher pole order | odd-measure-differential-hierarchy.md; odd-hierarchy-independent-review.md | Coefficient transfer only; target moments require higher zeta jets |
| Positive hyperbolic decomposition and explicit cotangent-root polynomial | odd-zeta-positive-hyperbolic-decomposition.md; hyperbolic-multiple-orthogonality.md; independent reviews | Individual component moments need not be rational |
| All balanced multiple-orthogonal polynomials are normal and rational with positive simple roots | hyperbolic-balanced-normality.md; hyperbolic-balanced-normality-independent-review.md | No unbalanced theorem, general coefficient formula, denominator rate or irrationality claim |

The balanced-polynomial constant-term formula is a conjecture. Its Bernoulli
antidifference reduction has a written proof awaiting separate full review.
The initial pure-linear-form experiments do not establish an asymptotic failure
theorem. No historical novelty assertion follows from an unsuccessful search.
