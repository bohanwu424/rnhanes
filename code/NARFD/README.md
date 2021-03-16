# NARFD

This README describes code that accompanies "Non-negative decomposition of functional count data" by Daniel Backenroth, Russell T. Shinohara, Jennifer A. Schrack, and Jeff Goldsmith.

## Data

The main analysis in the paper uses data from the Baltimore Longitudinal Study on Aging. Agreements with the agency that collected these data prevent us from releasing them publicly, but interested groups can apply to access these data here: https://www.blsa.nih.gov/how-apply. 

This supplement includes two simulated datasets. The first, which we describe as "artificial", generates data under the NARFD model using sine and cosine as functional prototypes; this is similar to the data used in the simualtions in the manuscript. Code to generate this dataset and others is also included in this supplement. 

The second simulated data is based on a study of accelerometer data in children. The primary analysis in "New insights into activity patterns in children, found using functional data analyses" was a function-on-scalar regression. The results of that analysis was used to generate observations that mimic real data: they use realistic covariates, regression coefficients based on a fitted model, and residuals based on those in a full analysis. In particular, these data were not generated under the NARFD model. 

Both datasets are in the `data` subdirectory, and are structured in "long" format, with columns for subject, timepoint, and observed value. 


## Code

The `source` subdirectory contains code used for fitting the NARFD model, as well as our novel implementation of GFPCA. There are three files:

* `NARFD.R` contains the `NARFD` function, and a description of that function. This is the primary user-facing function. 
* `NARFD_Utilities.R` contains various helper functions, including objective functions, functions for each step in our alternating minimization procedure, and a function to generate fitted values
* `Simulate_Data.R` contains two fuctions that are used in generating simulated data under NARFD and GFPCA models, and organizing these data for analysis. 

The code has been tested on R version 3.6.1. To install packages required for NARFD and the code in this supplement, please use:

```{r}
install.packages(
  "nloptr", "splines", "pbs", "refund", "NNLM", "dplyr", "tidyr", "data.table", 
  "glmnet", "mgcv", "tibble")
```


## Instructions for use

The `examples` subdirectory contains two reproducible reports with similar structures. These implement the NARFD model for the datasets described above, plot estimated prototypes, and show fitted values for a single subject. The artificial data example is smaller, and runs in less than a minute, while (on a somewhat under-powered laptop) the example that mimics accelerometer data runs in roughly one hour. 
