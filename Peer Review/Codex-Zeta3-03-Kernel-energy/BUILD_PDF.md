# Build the main explanation

`main-explanation.pdf` is the main explanatory document. Its editable source is `main-explanation.tex`.

From this directory, with a standard TeX Live installation including the packages declared in the source:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main-explanation.tex
pdflatex -interaction=nonstopmode -halt-on-error main-explanation.tex
```

If a changed table of contents shifts pagination, run once more until references stabilize. The submitted PDF has nine pages. The source compiles without LaTeX warnings, overfull boxes or missing characters. This typeset synthesis is new; the source proof notes and their internal AI reviews remain available alongside it for comparison.
