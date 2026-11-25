# R10 Weil lattices -- Builder A notes (running)

Start 08:24:47 UTC, hard limit 09:06:47 UTC.

## Code
- `weil.py`  : generic pipeline (Weil rep odd part in all embeddings, orbit BFS, exact Z-basis of Lambda_0 = Z[G]v0,
               integer Gram matrices B, E, exact det, order test, fixed lattices, t-grid Fincke-Pohst, Gamma).
- `lat.c`    : floating LLL (modified Gram-Schmidt on actual vectors, exact int64 transformation) + BKZ with enumeration.
- `fp_enum.c`: Fincke-Pohst (copied from R2).
- `run1.py p`: full pipeline for one prime; `search_full.py`: BKZ search on the full lattice over the t-grid;
  `verify37.py p`: independent exact verification of a violator incl. explicit rank-8 sublattice O.w.

## How the Z-basis is obtained (exactness)
orbit vectors (floats, 1e-14) -> choose 4h independent orbit vectors M (pivoted QR) -> integer Gram G_M (rounded, residual 1e-12)
-> exact adjugate (fraction-free Gauss-Jordan, verified G_M*adj = det*I exactly) -> exact M-coordinates of ~4n orbit vectors
-> exact HNF modulo det(G_M) (python ints) -> chain of prime-index steps with exact modular solves + float LLL (unimodular, int64)
-> all orbit vectors have integer coordinates (float check, <1e-10), S and T map the basis to integer unimodular combinations,
and EXACT identity det(B) * index^2 = det(G_M).  Hence the saved basis is a Z-basis of Z[G]v0.
