# This is a utility function to make the manual from a folder
# used for debugging
make_manual = function(){
# Having trouble with Makefile on Windows.  Switch to function
library(here)
  
  rm(list=ls())
  
  fils = c(
    "EMDerivation.Rnw",
    "Quick_Start.Rnw",
    "UserGuide.Rnw",
    "./manual_files/Quick_Examples.Rnw",
    "./manual_files/Case_Study_1.Rnw",
    "./manual_files/Case_Study_2.Rnw",
    "./manual_files/Case_Study_3.Rnw",
    "./manual_files/Case_Study_4.Rnw",
    "./manual_files/Case_Study_5.Rnw",
    "./manual_files/Case_Study_6.Rnw",
    "./manual_files/Case_Study_7.Rnw",
    "./manual_files/Case_Study_8.Rnw",
    "./manual_files/Case_Study_dlm1.Rnw",
    "./manual_files/Case_Study_mlr.Rnw",
    "./manual_files/Setting_Inits.Rnw",
    "./manual_files/MARp.Rnw",
    "./manual_files/Covariates.Rnw",
    "./manual_files/Manual.Rnw")
    
  #Sweave and Stangle
  for(fil in fils[11]){
    Sweave(here("vignettes", fil))
    Stangle(here("vignettes", fil))
    created.objs = ls(); created.objs=created.objs[created.objs!="fils"]
    rm(list=created.objs)
  }
  
  #make R files
  shell("cat ./figures/CS1--Cs1_*.R > ../inst/doc/Chapter_PVA.R")
  shell("cat ./figures/CS2--Cs2_*.R > ../inst/doc/Chapter_SealTrend.R")
  shell("cat ./figures/CS3--Cs*_*.R > ../inst/doc/Chapter_SealPopStructure.R")
  shell("cat ./figures/CS5--Cs5_*.R > ../inst/doc/Chapter_AnimalTracking.R")
  shell("cat ./figures/CS4--Cs*.R > ../inst/doc/Chapter_DFA.R")
  shell("cat ./figures/CS6--Cs*.R > ../inst/doc/Chapter_StructuralBreaks.R")
  shell("cat ./figures/CS7--Cs*.R > ../inst/doc/Chapter_SpeciesInteractions.R")
  shell("cat ./figures/CS8--Cs*.R > ../inst/doc/Chapter_CombiningTrendData.R")
  shell("cat ./figures/CSDLM--Cs_*.R > ../inst/doc/Chapter_UnivariateDLM.R")
  shell("cat ./figures/Covar--Covar_*.R > ../inst/doc/Chapter_Covariates.R")
  shell("cat ./figures/MLR--Cs_*.R > ../inst/doc/Chapter_MLR.R")
  shell("cat ./figures/MCI--Cs_mci_*.R > ../inst/doc/Chapter_inits.R")
  shell("cat ./figures/ARMAp--Cs_*.R > ../inst/doc/Chapter_MARp.R")

#make pdfs
  texfils = c("EMDerivation","Quick_Start","UserGuide","Manual")
  for(fil in texfils){
    system(paste0("pdflatex ", fil, ".tex"))
    system(paste0("bibtex ", fil))
    system(paste0("pdflatex ", fil, ".tex"))
  }
  shell("cat Manual.pdf > UserGuide.pdf")

#clean:
  system("rm -f *.log *.aux *.ilg *.ind *.idx *.blg *.bbl *.out *.Rout *.toc *.lof *.lot Rplots.ps Rplots.pdf")
  system("rm -f *.sty *.bst *.cls *.tex *.Rnw Manual.pdf")
  system("rm -f ./tables/*.* ./tex/*.* ./figures/*.* ./manual_files/*.*")
  system("rm -rf ./manual_files ./tex ./figures ./tables")
}