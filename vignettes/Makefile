RM = rm -f

all:	sweave stangle makeRfiles pdf

sweave: 
	"$(R_HOME)/bin/R" CMD Sweave EMDerivation.Rnw
	"$(R_HOME)/bin/R" CMD Sweave Quick_Start.Rnw
	"$(R_HOME)/bin/R" CMD Sweave UserGuide.Rnw
	"$(R_HOME)/bin/R" CMD Sweave Quick_Examples.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_1.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_2.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_3.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_4.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_5.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_6.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_7.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_8.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Case_Study_dlm1.xRnw
	"$(R_HOME)/bin/R" CMD Sweave ./tex/MARp.xRnw
	"$(R_HOME)/bin/R" CMD Sweave ./tex/Covariates.xRnw
	"$(R_HOME)/bin/R" CMD Sweave Manual.xRnw
	 
stangle: 
	"$(R_HOME)/bin/R" CMD Stangle EMDerivation.Rnw
	"$(R_HOME)/bin/R" CMD Stangle Quick_Start.Rnw
	"$(R_HOME)/bin/R" CMD Stangle UserGuide.Rnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_1.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_2.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_3.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_4.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_5.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_6.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_7.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_8.xRnw
	"$(R_HOME)/bin/R" CMD Stangle Case_Study_dlm1.xRnw
	"$(R_HOME)/bin/R" CMD Stangle ./tex/Covariates.xRnw
	"$(R_HOME)/bin/R" CMD Stangle ./tex/MARp.xRnw

makeRfiles:
	cat ./figures/CS1--Cs1_*.R > ../inst/doc/Chapter_PVA.R
	cat ./figures/CS2--Cs2_*.R > ../inst/doc/Chapter_SealTrend.R
	cat ./figures/CS3--Cs*_*.R > ../inst/doc/Chapter_SealPopStructure.R
	cat ./figures/CS5--Cs5_*.R > ../inst/doc/Chapter_AnimalTracking.R
	cat ./figures/CS4--Cs*.R > ../inst/doc/Chapter_DFA.R
	cat ./figures/CS6--Cs*.R > ../inst/doc/Chapter_StructuralBreaks.R
	cat ./figures/CS7--Cs*.R > ../inst/doc/Chapter_SpeciesInteractions.R
	cat ./figures/CS8--Cs*.R > ../inst/doc/Chapter_CombiningTrendData.R
	cat ./figures/CSDLM--Cs_*.R > ../inst/doc/Chapter_UnivariateDLM.R
	cat ./figures/Covar--Covar_*.R > ../inst/doc/Chapter_Covariates.R
	cat ./figures/ARMAp--Cs_*.R > ../inst/doc/Chapter_MARp.R

pdf:	
	texi2dvi --pdf --clean --verbose EMDerivation.tex
	texi2dvi --pdf --clean --verbose Quick_Start.tex
	texi2dvi --pdf --clean --verbose UserGuide.tex
	cat Manual.xRnw.tex > Manual.tex
	texi2dvi --pdf --clean --verbose Manual.tex
	cat Manual.pdf > UserGuide.pdf
	rm Manual.pdf
	rm *.xRnw.tex
	rm Manual.tex
	make clean

clean:
	$(RM) *.log *.aux *.ilg *.ind *.idx *.blg *.bbl *.out *.Rout *.toc *.lof *.lot Rplots.ps Rplots.pdf
	$(RM) *.sty *.bst *.cls
	$(RM) *.tex
	$(RM) *.xRnw *.xRnw.R
	$(RM) ./tables/*.*
	$(RM) ./tex/*.tex
	$(RM) ./tex/*.bib
	$(RM) ./tex/*.xRnw
	$(RM) ./figures/*.*
	$(RM) *.Rdata *RData
	rm -rf ./tex
	rm -rf ./figures
	rm -rf ./tables
	