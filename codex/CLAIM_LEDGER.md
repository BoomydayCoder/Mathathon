# Claim ledger

All 'proved' labels refer to arguments written and checked in this workspace,
not historical novelty. Standard external inputs are identified in each note.
Finite computation is supplementary evidence only.

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

## Current-source claims, not independently established here

Fauzan's September2026 preprint states irrationality of zeta5; inspected formal
repositories report completed proofs. We have not built Lean or audited all
31 pages. Exact versions and inspected assumptions are in
sources/odd-zeta-provenance.md. Current-source reporting must not be confused
with our independent mathematical results above.
