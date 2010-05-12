default: bitree.pdf

notes: notes.tex prob.bib
	pdflatex notes
	bibtex notes
	pdflatex notes
	pdflatex notes

bitree.pdf: bitree.asy	
	asy -f pdf -noV bitree.asy

clean:
	rm -rf *.aux
	rm -rf *.bbl