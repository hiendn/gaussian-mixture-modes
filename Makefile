.PHONY: setup reproduce verify manuscript
setup:
	Rscript scripts/bootstrap.R
reproduce:
	Rscript scripts/run_all.R
verify:
	Rscript scripts/verify.R
manuscript:
	latexmk -cd -xelatex -interaction=nonstopmode -halt-on-error manuscript/manuscript.tex
