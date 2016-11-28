### CAN BE CHANGED
# install path
INSTALLPATH=/export/bin
# Default queue name for PBS
QUEUE=batch

### DO NOT CHANGE
SHELL=bash
GITREPO=$(shell pwd)
GITREPOSED=$(shell pwd | sed 's/\//\\\//g')
GITVERSION=$(shell git describe --tags | sed 's/^v//g')

progs = gwEvol-paml gwEvol-report allInOne.py getResults.py positiveSelectionTest.R geneGroupFilter.py checkAlign.py
bins=$(addprefixbin/,$(progs))

### makefile core
all : $(bins)

$(bins) : bin 
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" $(subst bin/,scripts/,$@) | \
	sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" | \
	sed -e "s/SEDMATCHQUEUE/${QUEUE}/g"  > $@
	chmod 755 $@

bin :
	mkdir bin
	chmod 755 bin

.PHONY : install
install : $(bins)
	ln -s ${GITREPO}/bin/gwEvol-PS ${INSTALLPATH}
	ln -s ${GITREPO}/bin/gwEvol-PS-report ${INSTALLPATH}

.PHONY : force_install
force_install : $(bins)
	ln -sf ${GITREPO}/bin/gwEvol-PS ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/gwEvol-PS-report ${INSTALLPATH}

.PHONY : uninstall
uninstall : 
	rm ${INSTALLPATH}/gwEvol-PS
	rm ${INSTALLPATH}/gwEvol-PS-report

.PHONY : doc
doc : 
	Rscript -e "rmarkdown::render('README.md')"

.PHONY : clean
clean :
	-rm -R bin
	-rm README.html
