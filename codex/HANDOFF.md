# Handoff

Active autonomous research; user asks to keep going without permission pauses.
Target originally Nov14 2026 08:00 Los Angeles (16:00 UTC). User authorizes
remaining credits, no automatic resets. No balance/meter API is available.

## Read first

- problem131/discovery-guide.md: complete motivated recurrence-only zeta3 proof
  and integral supplement.
- problem131/research-results.md: theorem-level discoveries and limitations.
- problem131/fauzan-equilibrium-guide.md: new analytic simplification with
  an explicit positive equilibrium measure replacing sixteen trial measures.
- CHECKPOINT.md and CLAIM_LEDGER.md: status and proof distinctions.

## Strongest recent developments

The equilibrium measure has a one-sign-change positivity proof, globally
maximizes the required potential on its support, and has uniquely determined
endpoints. certify_equilibrium_endpoints.py uses exact dyadic interval arithmetic.
certify_equilibrium_energy.py and certify_equilibrium_energy_elementary.py use
two distinct reductions with Arb and certify energy objective < -69/50. They
share Arb and the endpoint enclosure; they are not independent software stacks.
Run energy scripts with PYTHONPATH=/workspace/scratch/odd-zeta-python.

common-odd-zeta-analytic-budget.md proves the same leading analytic exponent
for each fixed odd s. It does NOT prove a matching arithmetic estimate.
The source's existing arithmetic bound together with our analytic replacement
would strengthen its final squeeze, but full local arithmetic is being audited.
odd-functional-shift.md proves the differential shift and explains its pole and
uniform-normalization cost; the determinant may have further cancellations.

## Ongoing tasks

1. Odd-zeta agent: simpler arithmetic integral bound <1.38.
2. Simplify agent: local-prime/rank-allocation arithmetic audit.
3. Zudilin exposition agent: further design classification/minimality.
4. Parent: discrete Green-function/electrical interpretation of Apéry recurrence.

## Conventions and source limits

Zudilin F_n=2(A_n*zeta(3)-B_n), A0=1,A1=5,B0=0,B1=6.
Earlier auxiliary-rational-function.md uses half this normalization.
Original recurrence's final term is u_(n-1), as determined by its proof.
Exact SAIR problem statements were not retrieved; HTML shells are not statements.
Fauzan September2026 preprint and formalization reports were retrieved; no local
Lean build or complete independent preprint validation. Earlier language that
zeta5 is open was outdated relative to current sources.

All claims from other agents require review. No historical novelty claimed.
Work only in codex/, branch codex/independent-research. No Claude files changed;
a read-only tree search found no published HANDOFF on research-claude at e364df8.

## Latest continuation

Read sixth checkpoint. Connections agent is performing full simplified-zeta5
assembly audit. IMPORTANT: relaxed plain-rank proof uses new tilde-m normalizer
with outer exponent2sum(weights)-rank; do not use the source's stronger
min(rank,zero-count) exponent while claiming that refinement was eliminated.
The new budget is2203747/1600000 and final eventual exponent4n².
Integral agent developed primary-source modular/elliptic route; exposition
agent now independently auditing it. Simplify agent developing localizing pencils
and exact enclosing intervals; odd-zeta agent exploring changed proportions.
No new all-odd irrationality theorem. No historical novelty claim.

## Seventh-stage strongest results and continuation

- discovery-guide now derives17and5 fromNewton interpolation and detailedbalance.
- fauzan-simplified-proof-audit.md is the completed assembly; no unresolved
  bridge found at mathematical-review level, not formal verification/consensus.
- rational-energy-certificate.md and rational-energy-review.md replace the
  numerical package dependency with a standard-library certificate and give
  eventualexp(-13n²) for the simplified normalizer.
- one-pole-schur-structure.md + independentreview: exactWilson endpoint formulas
  and a genuine reduced-denominator obstruction.
- apery-modular-elliptic.md + review: complete advanced alternative proof using
  standardmodular-curve facts, and a precise restrictedhigher-weight obstruction.
- apery-cm-character.md now uses an elementary ODE deck-symmetry proof as main
  route; CM is optional interpretation. Independentfinalreview inprogress.
- problem137/prime-power-sylow-audit.md and prime-power-orbit-review.md contain
  the fully reviewedgroup proof, typedcounting andleading-unit formula.

Runtime is still active. User authorizes remaining credits, no automatic resets;
no usage meter API exists. Keep bounded discovery/review tasks substantive.
