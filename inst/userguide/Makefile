RM = rm -f

all:	cleandoc sweave stangle makeRfiles pdf compact clean

cleandoc: 
		rm -r -f ./doc/
		rm -d -f ./doc/
		mkdir ./doc/

sweave: 
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
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Structural_TS.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/KFAS.Rnw
	"$(R_HOME)/bin/R" CMD Sweave ./manual_files/Manual.Rnw
	 
stangle: 
	Rscript -e "Stangle('./manual_files/Quick_Examples.Rnw')"
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
	Rscript -e "Stangle('./manual_files/Structural_TS.Rnw')"
	Rscript -e "Stangle('./manual_files/KFAS.Rnw')"

makeRfiles:
	cat ./figures/QE--Cs*.R > ./doc/Quick_Examples.R
	cat ./figures/CS1--Cs1_*.R > ./doc/Chapter_PVA.R
	cat ./figures/CS2--Cs2_*.R > ./doc/Chapter_SealTrend.R
	cat ./figures/CS3--Cs*_*.R > ./doc/Chapter_SealPopStructure.R
	cat ./figures/CS4--Cs*.R > ./doc/Chapter_DFA.R
	cat ./figures/CS5--Cs*.R > ./doc/Chapter_AnimalTracking.R
	cat ./figures/CS6--Cs*.R > ./doc/Chapter_StructuralBreaks.R
	cat ./figures/CS7--Cs*.R > ./doc/Chapter_SpeciesInteractions.R
	cat ./figures/CS8--Cs*.R > ./doc/Chapter_CombiningTrendData.R
	cat ./figures/CSDLM--Cs_*.R > ./doc/Chapter_UnivariateDLM.R
	cat ./figures/Covar--Covar_*.R > ./doc/Chapter_Covariates.R
	cat ./figures/MLR--Cs_*.R > ./doc/Chapter_MLR.R
	cat ./figures/MCI--Cs_mci_*.R > ./doc/Chapter_inits.R
	cat ./figures/ARMAp--Cs_*.R > ./doc/Chapter_MARp.R
	cat ./figures/STS--Cs*.R > ./doc/Chapter_Structural_TS.R
	cat ./figures/KFAS--Cs*.R > ./doc/Chapter_KFAS.R

pdf:	
	pdflatex Manual.tex
	bibtex Manual
	pdflatex Manual.tex
	pdflatex Manual.tex
	makeindex Manual
	pdflatex Manual.tex
	
compact:
	Rscript -e "tools::compactPDF('Manual.pdf', gs_quality = 'ebook')"
	cat Manual.pdf > ./doc/UserGuide.pdf

clean:
	$(RM) *.log *.aux *.ilg *.ind *.idx *.blg *.bbl *.out *.Rout *.toc *.lof *.lot Rplots.ps Rplots.pdf
	$(RM) *.tex
	$(RM) *.R
	$(RM) Manual.pdf

