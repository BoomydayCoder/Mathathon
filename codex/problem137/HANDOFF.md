# Secondary project checkpoint

Read valuation.md first. It contains complete proofs of forced-carry cancellation and the full-ary-tree cycle quotient, a p-unit corollary, explicit limitations, and concrete construction directions. No historical novelty claim. Official problem #137 exact wording remains blocked/unverified; these results must not be described as answering its precise requirements yet.

Run `python codex/problem137/check_valuations.py` from repository root: 36,120 exact checks cover valuation, forced carry, and finite-transducer formulas. Source lead: Raney (1960), DOI verified via Crossref, full text not read. No Claude files read or modified. This agent made no commit.

Next priority: obtain exact statement and primary literature, then judge whether a tree action proving valuation is required. Current proof is elementary arithmetic plus a separate combinatorial integrality proof; the missing bridge is explicitly identified.

Follow-up: prime-power-arity.md gives a complete classification: for q=p^e, p-unit full q-ary tree counts occur exactly at base-q repunit indices (including zero). check_valuations.py now verifies this as well. The broader Apéry bridge, polynomial Lucas factorization and reflection symmetry, plus formal transducer audit are in ../problem131/apery-digit-arithmetic.md and ../computations/check_apery_digits.py.

Tree-action follow-up: tree-pgroup-proof.md proves the p-unit classification by finite wreath-product fixed points and gives an exact prime-arity orbit proof plus a cyclic-q4 counterexample. prime-power-sylow-audit.md audits the parent-derived full Sylow extension: exact valuations for all q=p^e via depth-constrained skeletons and the explicitly specialized Lagrange–Good coefficient theorem. The latter theorem is stated as machinery rather than silently assumed proved. check_tree_orbits.py supplies recursive small-size orbit enumerations and the q4,n28 counterexample.
