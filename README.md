# vimco

**vimco** package proves a Bayesian variable selection method for GWAS data with multiple traits. Unlike in BVSR where each trait is analyzed seperately, **vimco** performs a joint analysis for the multiple traits, while accounting for correlation among the multiple traits.
 
## Installation
To install the development version of **vimco**, it's easiest to use the 'devtools' package. Note that REMI depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```{r, fig.show='hold', eval=FALSE}
library(devtools)
install_github("XingjieShi/VIMCO")
```
## Usage

```{r, fig.show='hold', eval=FALSE}
library(vimco)
?vimco
```
## Replicate simulation results in Shi et al. (2018)
All the simulation results can be reproduced by using the code at [simulation](https://github.com/XingjieShi/VIMCO/tree/master/simulation). 
Before running simulation to reproduce the results, please familarize yourself with **vimco** using 'VIMCO' vignette. Simulation results can be reproduced by following steps in the [**vimco** vignette](https://github.com/XingjieShi/VIMCO/blob/master/vignettes/vimco.Rmd).


## Reference
[VIMCO: Variational Inference for Multiple Correlated Outcomes in Genome-wide Association Studies](https://arxiv.org/abs/1807.10467)
