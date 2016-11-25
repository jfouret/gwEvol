SHELL=sh
GITREPO=$(shell pwd)
GITREPOSED=$(shell pwd | sed 's/\//\\\//g')
GITVERSION=$(shell git describe --tags | sed 's/^v//g')
INSTALLPATH=/usr/local/bin/
all:
	mkdir -p bin
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/positiveSelection.py 		| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/gwEvol-PS
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/positiveSelectionReport.py 	| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/gwEvol-PS-report
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/allInOne.py 					| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/allInOne.py
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/getResults.py 				| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/getResults.py
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/positiveSelectionTest.R 		| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/positiveSelectionTest.R
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/geneGroupFilter.py 			| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/geneGroupFilter.py
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/checkAlign.py 				| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/checkAlign.py
	chmod 755 -R bin
install :
	ln -s ${GITREPO}/bin/gwEvol-PS ${INSTALLPATH}
	ln -s ${GITREPO}/bin/gwEvol-PS-report ${INSTALLPATH}
force_install :
	ln -sf ${GITREPO}/bin/gwEvol-PS ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/gwEvol-PS-report ${INSTALLPATH}
uninstall :
	rm ${INSTALLPATH}/gwEvol-PS
	rm ${INSTALLPATH}/gwEvol-PS-report
doc :
	Rscript -e "rmarkdown::render('README.md')"
clean :
	-rm -R bin
	-rm README.html
