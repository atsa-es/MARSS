RM = rm -f

all:	sweave stangle makeRfiles pdf

sweave: 
	"$(R_HOME)/bin/R" CMD Sweave EMDerivation.Rnw
	"$(R_HOME)/bin/R" CMD Sweave Quick_Start.Rnw
	"$(R_HOME)/bin/R" CMD Sweave UserGuide.Rnw
	"$(R_HOME)/bin/R" CMD Sweave Residuals.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Quick_Examples.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_1.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_2.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_3.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_4.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_5.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_6.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_7.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_8.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_dlm1.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Case_Study_mlr.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Setting_Inits.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/MARp.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Covariates.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Manual.Rnw
	 
stangle: 
	Rscript -e "Stangle('./manual_files/Case_Study_1.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_2.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_3.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_4.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_5.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_6.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_7.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_8.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_dlm1.Rnw')"
	Rscript -e "Stangle('./manual_files/Case_Study_mlr.Rnw')"
	Rscript -e "Stangle('./manual_files/Setting_Inits.Rnw')"
	Rscript -e "Stangle('./manual_files/Covariates.Rnw')"
	Rscript -e "Stangle('./manual_files/MARp.Rnw')"

makeRfiles:
	cat ./figures/CS1--Cs1_*.R > ../inst/doc/Chapter_PVA.R
	cat ./figures/CS2--Cs2_*.R > ../inst/doc/Chapter_SealTrend.R
	cat ./figures/CS3--Cs*_*.R > ../inst/doc/Chapter_SealPopStructure.R
	cat ./figures/CS5--Cs*.R > ../inst/doc/Chapter_AnimalTracking.R
	cat ./figures/CS4--Cs*.R > ../inst/doc/Chapter_DFA.R
	cat ./figures/CS6--Cs*.R > ../inst/doc/Chapter_StructuralBreaks.R
	cat ./figures/CS7--Cs*.R > ../inst/doc/Chapter_SpeciesInteractions.R
	cat ./figures/CS8--Cs*.R > ../inst/doc/Chapter_CombiningTrendData.R
	cat ./figures/CSDLM--Cs_*.R > ../inst/doc/Chapter_UnivariateDLM.R
	cat ./figures/Covar--Covar_*.R > ../inst/doc/Chapter_Covariates.R
	cat ./figures/MLR--Cs_*.R > ../inst/doc/Chapter_MLR.R
	cat ./figures/MCI--Cs_mci_*.R > ../inst/doc/Chapter_inits.R
	cat ./figures/ARMAp--Cs_*.R > ../inst/doc/Chapter_MARp.R

pdf:	
	pdflatex EMDerivation.tex
	pdflatex Quick_Start.tex
	pdflatex Manual.tex
	pdflatex Residuals.tex
	bibtex EMDerivation
	bibtex Manual
	bibtex Residuals
	pdflatex EMDerivation.tex
	pdflatex Quick_Start.tex
	pdflatex Residuals.tex
	pdflatex Manual.tex
	pdflatex EMDerivation.tex
	pdflatex Quick_Start.tex
	pdflatex Residuals.tex
	pdflatex Manual.tex
	makeindex Manual
	pdflatex Manual.tex
	cat Manual.pdf > UserGuide.pdf

clean:
	$(RM) *.log *.aux *.ilg *.ind *.idx *.blg *.bbl *.out *.Rout *.toc *.lof *.lot Rplots.ps Rplots.pdf
	$(RM) *.sty *.bst *.cls
	$(RM) *.tex
	$(RM) *.R
	$(RM) ../README.md ../_config.yml ../favicon.ico ../logo.jpg ../logot.png
	$(RM) ../_layouts/*.*
	$(RM) Manual.pdf
	$(RM) ./tables/*.*
	$(RM) ./tex/*.*
	$(RM) ./figures/*.*
	$(RM) ./manual_files/*.*
	$(RM) debug_manual.R debugger_template.Rnw
	rm -rf ./manual_files
	rm -rf ./tex
	rm -rf ./figures
	rm -rf ./tables
	rm -rf ../_layouts
	$(RM) ../R/toLatex_marssMODEL.R ../R/predict_marssMLE.R ../man/predict_marssMLE.Rd
	