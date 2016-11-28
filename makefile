### CAN BE CHANGED
# install path
INSTALLPATH=/export/bin
# Default queue name for PBS
QUEUE=batch
#default path for ete3 software
ETE3=/export/source/archive/anaconda_ete/bin/

### DO NOT CHANGE
SHELL=bash

GITREPO=$(shell pwd)
GITREPOSED=$(shell pwd | sed 's/\//\\\//g')
GITETE3=$(shell echo ${ETE3} | sed 's/\//\\\//g')

GITVERSION=$(shell git describe --tags | sed 's/^v//g')

progs = gwEvol-paml gwEvol-report allInOne.py getResults.py positiveSelectionTest.R geneGroupFilter.py checkAlign.py
bins=$(addprefix bin/,$(progs))

prog_install = gwEvol-paml gwEvol-report
installs=$(addprefix ${INSTALLPATH}/,$(prog_install))

### makefile core
all : $(bins)

$(bins) : bin 
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" $(subst bin/,scripts/,$@) | \
	sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" | \
	sed -e "s/SEDMATCHETE3/${GITETE3}/g" | \
	sed -e "s/SEDMATCHQUEUE/${QUEUE}/g"  > $@
	chmod 755 $@

bin :
	mkdir bin
	chmod 755 bin

.PHONY : install
install : $(installs)
$(installs) : $(bins)
	ln -sf $(subst ${INSTALLPATH}/, ${GITREPO}/bin/, $@) ${INSTALLPATH}

.PHONY : uninstall
uninstall : 
	rm $(installs)

.PHONY : doc
doc : 
	Rscript -e "rmarkdown::render('README.md')"

.PHONY : clean
clean :
	-rm -R bin
	-rm README.html
