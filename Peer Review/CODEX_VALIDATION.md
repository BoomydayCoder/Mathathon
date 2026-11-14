# Codex packaging validation

The 26 check commands listed in the six Codex READMEs all exited successfully on 14 November 2026. Their complete output and per-script exit status are saved under each package's `verification/`. Finite checks do not replace the mathematical proofs or human review.

Run `python3 "Peer Review/CODEX_REPRODUCE.py"` from the repository root to rerun them. This updates the verification files. Python 3, SymPy and mpmath are required; the exact environment used is in `CODEX_ENVIRONMENT.json`.

All relative Markdown file links in the six packages were checked for existing targets. Additional dependencies link to the immutable Codex research commit. SHA-256 hashes in each source manifest distinguish original source files from copies with relocated links.

The main explanation is a nine-page PDF, built from its included LaTeX source. The final build had no LaTeX warnings, overfull boxes or missing-character messages. Text extraction succeeded; rendered title and mathematical pages were inspected. This is document validation, not an independent human mathematical audit.

Before publication, the change was checked to contain only new files under `Peer Review/`. Claude's pre-existing files were preserved byte for byte. The original Codex research checkout remains on its research branch.
