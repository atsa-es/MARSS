RM = rm -f
version = 3.11.9
GitHubdir = ~/Documents/GitHub
tmpdir = ~/Dropbox

all:	movetovignettes build movetodoc rebuild

manual:
	cd ${GitHubdir}/MARSS/inst/userguide && make
	cd ${GitHubdir}/MARSS/inst/derivations && make
  
movetovignettes:
	rm -r -f ${GitHubdir}/MARSS/inst/doc/
	rm -d -f ${GitHubdir}/MARSS/inst/doc/
	mkdir ${GitHubdir}/MARSS/inst/doc/
	cp -p ${GitHubdir}/MARSS/inst/derivations/*.Rnw ${GitHubdir}/MARSS/vignettes/
	cp -p ${GitHubdir}/MARSS/inst/derivations/*.bib ${GitHubdir}/MARSS/vignettes/
	cp -p ${GitHubdir}/MARSS/inst/userguide/UserGuide.Rnw ${GitHubdir}/MARSS/vignettes/
	cp -p ${GitHubdir}/MARSS/inst/userguide/doc/UserGuide.pdf ${GitHubdir}/MARSS/inst/doc/

build:
	rm -f ${GitHubdir}/MARSS_${version}.tar.gz
	cd $(GitHubdir) && "$(R_HOME)/bin/R" CMD build ${GitHubdir}/MARSS
	rm -r -f ${tmpdir}/MARSS
	rm -f ${tmpdir}/MARSS_${version}.tar.gz
	cp -p $(GitHubdir)/MARSS_${version}.tar.gz ${tmpdir}/MARSS_${version}.tar.gz
	cd $(tmpdir) && tar -xvzf ${tmpdir}/MARSS_${version}.tar.gz

movetodoc:
	rm -f ${tmpdir}/MARSS/inst/doc/*.pdf
	cp -a -p ${GitHubdir}/MARSS/inst/derivations/doc/*.pdf ${tmpdir}/MARSS/inst/doc
	cp -a -p ${GitHubdir}/MARSS/inst/userguide/doc/. $(tmpdir)/MARSS/inst/doc/
	cp -a -p ${GitHubdir}/MARSS/inst/derivations/doc/. ${GitHubdir}/MARSS/inst/doc/
	cp -a -p ${GitHubdir}/MARSS/inst/userguide/doc/. $(GitHubdir)/MARSS/inst/doc/

rebuild:
	cd $(tmpdir) && R CMD build --no-build-vignettes MARSS
