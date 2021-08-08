# covariate-powered-cross-weighted-multiple-testing

This repository serves as companion to the following paper:


> N. Ignatiadis and W. Huber (2021), **Covariate powered cross-weighted multiple testing**, Journal of the Royal Statistical Society: Series B (JRSS-B)

The paper is also available on [[arXiv]](https://arxiv.org/abs/1701.05179).


All numerical results and figures in the aforementioned paper are made third-party reproducible herein.

We note that the Bioconductor package [IHW](http://bioconductor.org/packages/IHW) provides a user-friendly implementation of the IHW-BH/Storey procedures with conditional distributions estimated with the Grenander estimator.

## Pre-rendered vignettes

First, below we provide links to pre-rendered vignettes that reproduce figures.

* hQTL data analysis example (Figs. 1 and 7) [[Vignette]](http://bioconductor.org/packages/devel/data/experiment/vignettes/IHWpaper/inst/doc/hqtl_IHW_BY.html)
* Simulations on grouped multiple testing (Figs. 2 and 3) [[Vignette]](http://htmlpreview.github.io/?https://github.com/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing/blob/master/vignettes/grouped_multiple_testing.html)
* Simulations on multiple testing with continuous covariates (Figs. 4 and 5) [[Vignette]](http://htmlpreview.github.io/?https://github.com/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing/blob/master/vignettes/betamix_simulations.html)
* Simulations on simultaneous two-sample testing (Fig. 6) [[Vignette]](http://htmlpreview.github.io/?https://github.com/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing/blob/master/vignettes/two_sample_testing.html)


## Subdirectories of this repository

### **IHWStatsPaper/**
[![R build status](https://github.com/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing/workflows/R-CMD-check/badge.svg)](https://github.com/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing/actions) 
[![codecov](https://codecov.io/gh/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing/branch/master/graph/badge.svg)](https://codecov.io/gh/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing)

This directory contains `IHWStatsPaper`, a R package wrapping/implementing the different methods compared, the simulation functions, as well as the benchmarking code. It can be installed as follows.
```r
devtools::install_github("Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing",
                         subdir="IHWStatsPaper")
```

The simulations in the paper were implemented based on the following [commit](https://github.com/Huber-group-EMBL/covariate-powered-cross-weighted-multiple-testing/tree/81784986c1975c476df3bb82317daa112047dca3) of the package.

### **scripts/**
R scripts that run the simulations  (see in folder for details).
 
### **precomputed_results/**
Simulation results from scripts in the *scripts* directory that have been pre-computed/cached.

### **vignettes/**
R markdown documents that produce the pre-rendered vignettes above upon being compiled (using files saved in the *precomputed_results* folder). The only exception is the vignette for the hQTL data analysis example, which we include as part of the Bioconductor package [IHWpaper](http://bioconductor.org/packages/devel/IHWpaper).
