# Description of files.


## Implemented multiple testing methods.

* **multiple_testing_methods.R**: Main file; loads all multiple testing methods benchmarked in the simulations studies.
* **beta_uniform_glm.R**: Implements EM algorithm to fit the conditional beta-uniform mixture model with or without censoring, as well as functions for the optimization to derive optimal weights through conditional local fdrs.
* **cars_weighting.R**: Implements a heuristic to convert CARS thresholds to weights.
* **groupwise_sabha_pi0.R**: A reimplementation of the ADMM algorithm to fit SABHA (Structure Adaptive Benjamini Hochberg algorithm) which is faster for a larger number of groups. The starting point for the reimplementation was the [original implementation](https://www.stat.uchicago.edu/~rina/sabha.html) by the authors of SABHA. .


## Code for simulations.

### Grouped multiple testing
* **grouped_simulations_functions.R**
* **grouped_simulations_null.R**
* **grouped_simulations.R**

### Multiple testing with continuous covariates
* **betamix_functions.R**
* **betamix_simulations.R**
* **betamix_onesided_simulations.R**

### Simultaneous two-sample t-tests
* **cars_ttest_functions.R**
* **cars_ttest_simulations.R**

