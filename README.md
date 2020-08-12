# IHWStatsPaper

This repository serves as companion to the following preprint:


> *Covariate-powered cross-weighted multiple testing*, N. Ignatiadis and W.Huber [[arXiv]](https://arxiv.org/abs/1701.05179)


All numerical results and figures in the aforementioned paper are made fully third-party reproducible herein.

We note that the Bioconductor package [IHW](http://bioconductor.org/packages/IHW) provides a user-friendly implementation of the IHW-BH/Storey procedures with conditional distributions estimated with the Grenander estimator.

## Pre-rendered vignettes

First, below we provide links to pre-rendered vignettes that reproduce figures.

* hQTL data analysis example (Figs. 1 and 7) [[Vignette]](http://bioconductor.org/packages/devel/data/experiment/vignettes/IHWpaper/inst/doc/hqtl_IHW_BY.html)
* Simulations on grouped multiple testing (Figs. 2 and 3) [[Vignette]](http://htmlpreview.github.io/?https://github.com/nignatiadis/IHWStatsPaper/blob/master/vignettes/grouped_multiple_testing.html)
* Simulations on multiple testing with continuous covariates (Figs. 4 and 5) [[Vignette]](http://htmlpreview.github.io/?https://github.com/nignatiadis/IHWStatsPaper/blob/master/vignettes/betamix_simulations.html)
* Simulations on simultaneous two-sample testing (Fig. 6) [[Vignette]](http://htmlpreview.github.io/?https://github.com/nignatiadis/IHWStatsPaper/blob/master/vignettes/two_sample_testing.html)


## Subdirectories of this repository

* **code/**: Main *.R* files that implement benchmarked methods and run simulations. See README.md of subdirectory for description of files within.
* **precomputed_results/**: Simulation results from scripts in the *code* directory that have been pre-computed/cached.
* **vignettes/**: R markdown documents that produce the pre-rendered vignettes above upon being compiled (using files saved in the *precomputed_results* folder). The only exception is the vignette for the hQTL data analysis example, which we include as part of the Bioconductor package [IHWpaper](http://bioconductor.org/packages/devel/IHWpaper).
* **IHWStatsPaper/**: The main functions in the *code* directory wrapped as a standalone R package. Many of these will be ported to the Bioconductor IHW package.[![R build status](https://github.com/nignatiadis/IHWStatsPaper/workflows/R-CMD-check/badge.svg)](https://github.com/nignatiadis/IHWStatsPaper/actions) 
[![codecov](https://codecov.io/gh/nignatiadis/IHWStatsPaper/branch/master/graph/badge.svg)](https://codecov.io/gh/nignatiadis/IHWStatsPaper)
