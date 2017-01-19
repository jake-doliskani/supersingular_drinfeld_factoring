all:
	pdflatex supersingular_drinfeld_factoring.tex
	bibtex   supersingular_drinfeld_factoring
	pdflatex supersingular_drinfeld_factoring.tex
	pdflatex supersingular_drinfeld_factoring.tex


