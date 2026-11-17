# Likelihood estimation of Gaussian mixtures with one mode per component

First manuscript draft, 17 November 2026.

## Files

- manuscript.tex: complete editable manuscript, including proofs.
- references.tex: author–year reference list used by the manuscript.
- First_Draft.pdf: compiled manuscript.
- Author_Note.md and Author_Note.pdf: separate notes on the scientific scope and empirical work still needed.
- figures/: the density-contour figure and its plotted coordinates.
- code/: figure generation and the existing numerical evidence.

## Building

Compile with XeLaTeX twice, or until the cross-reference rerun messages disappear:

    xelatex -interaction=nonstopmode -halt-on-error manuscript.tex
    xelatex -interaction=nonstopmode -halt-on-error manuscript.tex

The source uses standard LaTeX packages and includes the bibliography directly, so a BibTeX run is unnecessary. Compilation produces manuscript.pdf; the supplied First_Draft.pdf is the same final document under a descriptive filename.

The author's name is retained in this working copy. The draft is for discussion and development, not an anonymous submitted manuscript. The substantive application and larger comparative experiments are not yet present. See the separate author note for the precise next work.
