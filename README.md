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


