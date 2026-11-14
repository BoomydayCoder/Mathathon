# Odd-zeta source provenance

Retrieved 2026-11-13 America/Los_Angeles (2026-11-14 UTC). Hashes identify the bytes actually saved. No original-paper proof or Lean build is certified merely by successful retrieval.

| Local file in codex/sources | Source URL / derivation | Version | SHA-256 |
|---|---|---|---|
| fauzan-zeta5.pdf | https://zenodo.org/records/22826419/files/ZETA5_IS_IRRATIONAL.pdf | v1; paper dated 2026-09-17 | a42c05ccc46b9a7c36428913453c4617202c59109713ba20f464ab9ec640d0da |
| rivoal-0008051.pdf | https://arxiv.org/pdf/math/0008051 | v1; PDF identifies 2000-08-07 | f9921534f4a2fc32b4579ed01e9791db3b8b2b8545f669385efd4957eaf927c6 |
| fischler-lai-sprang-zudilin-2026.pdf | https://arxiv.org/pdf/2611.16156 | v1; PDF identifies 2026-11-12 | 7cc92106416e6347820130c02140ea695b0b0f6a0b3b8ec12c0cf5256729bfd9 |
| fauzan-zeta5.txt | Local pdftotext -layout extraction of fauzan-zeta5.pdf | v1; paper dated 2026-09-17 | 7f1188c58fa57c7287de49bd7d11bef5fa2bafae22482b2cf9deaaa895566304 |
| rivoal-0008051.txt | Local pdftotext -layout extraction of rivoal-0008051.pdf | v1; PDF identifies 2000-08-07 | 666395596080eb366561e7e1a7297d820a8be8ee69d63904c4e1e7b9c8c0eb88 |
| fischler-lai-sprang-zudilin-2026.txt | Local pdftotext -layout extraction of fischler-lai-sprang-zudilin-2026.pdf | v1; PDF identifies 2026-11-12 | 979e080cf210487ef02128087a81ebd4e2c15816f1bde48b774535ea9e54c66c |
| romik-README.md | https://raw.githubusercontent.com/danromik/zeta5-irrationality/fc2fefa77076e22e7e80b4c101ded5a98d53e14d/README.md | fc2fefa77076e22e7e80b4c101ded5a98d53e14d | ccb3028962a9bb48d8a04d05c9847aa367e5cd4831fe94ad54b02cde853c93f2 |
| romik-Zeta5-Axioms.lean | https://raw.githubusercontent.com/danromik/zeta5-irrationality/fc2fefa77076e22e7e80b4c101ded5a98d53e14d/Zeta5/Axioms.lean | fc2fefa77076e22e7e80b4c101ded5a98d53e14d | 654d5be5931bf03de1eb3170b4f2a04f88d5e3e2bd48656e521753398639b324 |
| romik-Zeta5-Interface.lean | https://raw.githubusercontent.com/danromik/zeta5-irrationality/fc2fefa77076e22e7e80b4c101ded5a98d53e14d/Zeta5/Interface.lean | fc2fefa77076e22e7e80b4c101ded5a98d53e14d | 1f0e948f43a06fe866265c26c6358a619af0c5dc4aa5ce22ed213c2646306113 |
| romik-cert-Check.lean | https://raw.githubusercontent.com/danromik/zeta5-irrationality/fc2fefa77076e22e7e80b4c101ded5a98d53e14d/cert/Check.lean | fc2fefa77076e22e7e80b4c101ded5a98d53e14d | 2ed7594c96bc1920139ea7172f266f6d3abc196cb135df4799e218cd452021d2 |

## Verification status

* Fauzan: full PDF retrieved; introduction, theorem statements and moment construction read. Its central arithmetic/energy estimates have not been completely independently checked.
* Romik: **local build and axiom inspection reproduced on 2026-11-14 UTC.** The full pinned commit compiled with Lean 4.34.0 and its pinned Mathlib (8995 jobs, exit 0). Both the repository check and an independently written, fully qualified real-series check compiled. The actual axiom report is exactly `propext`, `Classical.choice`, `Quot.sound`, and `Zeta5.Axioms.chebyshev_theta_asymptotic`. The latter is the standard PNT statement θ(x)/x→1. Thus the ordinary ζ(5) irrationality statement has been locally checked conditional on that external theorem. A separate project kernel replay is running; no completed replay claim is made yet. See `../problem131/lean-local-verification.md` and logs/metadata in `../verification/romik-lean/`. This verifies the pinned formalization, not the new Python equilibrium certificate developed in this research session.
* Rivoal: Theorem 1 checked against retrieved PDF; full proof not rederived.
* Fischler–Lai–Sprang–Zudilin: current theorem statement and explicit use of Fauzan inspected; full proof not rederived.

## Correction to the earlier task framing

The earlier statement that ζ(5) irrationality remains open is outdated as a description of current reported research. Fauzan has a 2026 primary proof claim, Romik reports a Lean formalization conditional only on the standard PNT, and a November 2026 paper by Fischler, Lai, Sprang and Zudilin explicitly uses the method. Our own verification is limited as above; do not substitute either unconditional endorsement or blanket dismissal for those facts. The project goal should be an independently verified simplification, not an unsupported first-proof claim.
