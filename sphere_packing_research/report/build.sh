#!/bin/bash
# Build the final report PDF (XeLaTeX + pandoc for the markdown appendices)
set -e
cd "$(dirname "$0")"
pandoc -f markdown+tex_math_dollars+lists_without_preceding_blankline-smart -t latex --wrap=preserve ../rounds/R6_quaternionic/proof_revised.md -o app_proof.tex
pandoc -f markdown+tex_math_dollars+lists_without_preceding_blankline-smart -t latex --wrap=preserve ../research_log.md -o app_log.tex
# demote headings inside appendices so they do not clash with the main sectioning
sed -i 's/\\section{/\\subsection*{/; s/\\subsection{/\\subsubsection*{/; s/\\subsubsection{/\\paragraph{/' app_proof.tex
sed -i 's/\\section{/\\subsection*{/; s/\\subsection{/\\subsubsection*{/; s/\\subsubsection{/\\paragraph{/' app_log.tex
# rounds 14/15: proofs written as markdown+LaTeX by sub-agents (included as appendices when present)
for pair in "../rounds/R14_leech_tensor/proof.md:app_proof14.tex" "../rounds/R15_rogers/proof.md:app_proof15.tex"; do
  src="${pair%%:*}"; dst="${pair##*:}"
  if [ -f "$src" ]; then
    pandoc -f markdown+tex_math_dollars+lists_without_preceding_blankline-smart -t latex --wrap=preserve "$src" -o "$dst"
    sed -i 's/\\section{/\\subsection*{/; s/\\subsection{/\\subsubsection*{/; s/\\subsubsection{/\\paragraph{/' "$dst"
  else
    : > "$dst"
  fi
done
./build_notes.sh
xelatex -interaction=nonstopmode -halt-on-error main.tex > build1.log 2>&1 || { tail -40 build1.log; exit 1; }
xelatex -interaction=nonstopmode -halt-on-error main.tex > build2.log 2>&1 || { tail -40 build2.log; exit 1; }
grep -c "Missing character" main.log || true
grep -n "undefined" main.log | head -20 || true
