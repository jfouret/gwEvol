# INSTALLATION

## Fast installation

make all  
sudo make install  

## More

make clean 	==> cleanup the repository and reset the original state  
sudo make uninstall ==> uninstall executables from the $PATH  
sudo make force_install ==> overwriting old link in $PATH  
make doc ==> transform this document in html

# List of softwares

Executable			|Descrition
--------------------|--------------------
gwEvol-PS			| Run Positive selection model via PAML
gwEvol-PS-report	| Run all scripts to get the report and the analyze of a positive selection analysis

## gwEval-PS


# Usage

For help run executable with -h or --help  

# Dependencies

*Python 2.7 (or higher but never tested)
	*jupype
*R
	*markdown
	*ggplot2
	*knitr
*ete3
*PAML
*PBS serveur with qsub
