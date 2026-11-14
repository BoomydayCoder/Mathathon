# Problem #137 checkpoint

The exact official statement is retrieved in ../sources/problem137-official.md (raw JSON/provenance nearby). It asks precisely for a combinatorial exact-valuation proof at prime-power arity, beyond Konvalinka's orbit divisibility bound. Initial access failures are resolved.

## Current mathematical result

For q=p^e, L=(q-1)n+1=sum_j d_j p^j and D=(sum_j d_j-1)/(p-1), the developed proof gives

v_p(T_q(n))=D,
p^{-D}T_q(n) congruent (-1)^D/product_j d_j! modulo p.

The construction expands q-vertices into perfect p-blocks, acts by the full Sylow p-group, and counts minimum orbits using depth-constrained digit-labeled skeletons. cyclic-type-cycle-splicing.md supplies an explicit type-preserving functional-digraph bijection for their enumeration, eliminating the initial Lagrange–Good machinery. prime-power-orbit-review.md and prime-power-sylow-audit.md give independent proof audits, including zero-type and single-leaf cases. typed-count-matrix-tree-audit.md independently verifies the same enumeration by Matrix–Tree.

A useful simpler entry is tree-pgroup-proof.md: fixed points prove the unit-index classification, and a cyclic-four-child counterexample at n=28 explains why the full Sylow action is needed. prime-power-arity.md gives an independent arithmetic check. valuation.md contains the earlier forced-carry formula and finite-state transducer, now background rather than the endpoint of the project. discovery-guide.md is the concise project-facing route.

## Literature status

konvalinka-comparison.md reads the primary paper and pinpoints Proposition2's arithmetic exactness versus Proposition7's orbit bound. The binary wreath-action argument is established Deutsch–Sagan literature, explicitly attributed. bounded-literature-comparison.md covers five further promising sources/chains: four full texts read in relevant parts and one MIT thesis record whose PDF returned429. Gao–Gu2021 gives substantial weighted prime-power congruences. No historical novelty claim is justified for the current proof assembly; its status remains unresolved.

Primary PDFs/text and hashes are recorded in ../sources/tree-primary-provenance.md. No Claude-owned files were read or modified by this agent.

## Proven weighted special case

depth-weighted-extension.md proves that ordinary-depth weights b_h=u modulo p, with u a p-unit, preserve valuation D and multiply the leading unit by u^n. This differs from the cited left-/non-right-edge path-weight models; those weights generally vary within symmetry orbits and require averaging/difference hypotheses.

## Reproducible checks already passed

* check_valuations.py:36,120 carry/transducer checks plus1,812 prime-power unit-index checks.
* check_tree_orbits.py:recursive cyclic/Sylow examples and the q4,n28 counterexample.
* check_prime_power_sylow.py:3,216 exact valuation/leading-unit checks and exact minimum-orbit enumeration for small q2,3,4,8,9 trees.
* check_prime_power_orbit_review.py:independent depth-profile skeleton counts, including zero internal types (connections agent).
* check_cycle_splicing.py:24,242 slot configurations plus inverse/surjectivity checks (connections agent).
* check_depth_weights.py:254 independent weighted-recursion checks.

Apéry connections are in ../problem131/apery-digit-arithmetic.md and its computation script: carries kill cross-digit binomial terms, yielding Lucas factorization, while tree quotients/orbits explain a different digit-sum mechanism here.

Next useful work is exposition, a bounded specialist novelty check, or a carefully defined weighted extension. Do not reopen a claim of incompleteness merely because the early notes originally stopped at integrality; the later cycle-splicing proof completes the required enumeration. This agent has not committed changes; the parent manages the shared branch and root checkpoint.
