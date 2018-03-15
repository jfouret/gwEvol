# INSTALLATION

## Fast installation

```
make all  
sudo make install  
```
## More

* `make clean` 	==> cleanup the repository and reset the original state  
* `sudo make uninstall` ==> uninstall executables from the $PATH  
* `sudo make force_install` ==> overwriting old link in $PATH  

# List of softwares

Executable			|Descrition
--------------------|--------------------
gwEvol-paml			| Run Positive selection model via PAML
gwEvol-result	| Run all scripts to get the report and the analyze of a positive selection analysis

# Usage

For help run executable with -h or --help  

# Controling jobs execution
## Error during job completion related to the scheduler work

```

# number of model tested
model=5 
expected=$(echo "3+2*$model" |bc)

for pbs in $(ls schd/gwEvol_paml_*) ; do echo "$pbs" ; status=""; for dir in $(more $pbs | grep 'ete3 evol' | cut -f 1 -d ';' | cut -f 2 -d ' ' | uniq ) ; do status="${status}$(ls $dir | wc -l):" ; done ;echo "$status" ;done | grep -P ':(?!'"$expected"')\d+' -B 1 --no-group-separator 

```

## Error during the job completion related to the job script

```
for errFile in $(find -name '*.err') ; do [[ -s $errFile ]] && echo $errFile >> errFile.list; done

```
# Dependencies

* Python 2.7 (or higher but never tested)
	* upype
* R
	* markdown
	* ggplot2
	* knitr
* ete3 ( easy to install with conda )
* PAML
* TORQUE or SLURM scheduler PBS compatible


