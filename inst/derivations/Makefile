RM = rm -f

all:	sweave pdf compress clean

sweave: 
	"$(R_HOME)/bin/R" CMD Sweave EMDerivation.Rnw
	"$(R_HOME)/bin/R" CMD Sweave Residuals.Rnw

pdf:	
	pdflatex EMDerivation.tex
	pdflatex Residuals.tex
	bibtex EMDerivation
	bibtex Residuals
	pdflatex EMDerivation.tex
	pdflatex Residuals.tex
	pdflatex EMDerivation.tex
	pdflatex Residuals.tex
	pdflatex EMDerivation.tex
	pdflatex Residuals.tex

	
compress:
	Rscript -e "tools::compactPDF('EMDerivation.pdf', gs_quality = 'ebook')"
	Rscript -e "tools::compactPDF('Residuals.pdf', gs_quality = 'ebook')"

clean:
	$(RM) *.log *.aux *.ilg *.ind *.idx *.blg *.bbl *.out *.Rout *.toc *.lof *.lot
	$(RM) *.tex


https://iugo-cafe.org/