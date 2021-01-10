# Code
We provide the full set of analyses code in this GitHub repository that can be used to replicate all the results in the main paper. Specific instructions are provided in the [Reproducibility workflow](https://github.com/trambakbanerjee/crejm-code#reproducibility-workflow) below on how to replicate the tables and figures in the main paper. Additionally, the code rely on specific software and hardware requirements that have been described below in detail. 

## Description

There are 3 folders in this repository:

1. data - holds the data files that were used for all the analyses in this paper. See the [README file](https://github.com/trambakbanerjee/crejm-code/tree/main/data#data) inside this folder for more information.
2. spcov - this folder has the R functions from the package `spcov` (Bien, J., and Tibshirani, R. (2011)). 
3. library - holds the main R functions. Please see the [README file](https://github.com/trambakbanerjee/crejm-code/tree/main/library#description) inside this folder for more information. 

Other than these 3 folders, there are 8 scripts in this repository. We describe them below:

1. `dataprocessing.R` - use this script to read the raw data files and output `out.RData` which is a list that holds the processed training and prediction data. 

2. `datasummary.R` - use this script to reproduce table 4 in the main paper. The script writes out two `CSV` files that hold the summary statistics in table 4. 

3. `motivatingfigures.R` - use this script to reproduce figures 2 and 4 in the main paper. 

# Reproducibility workflow
All figures and tables in the paper are reproducible except for figures 1 and 3, and table 3. These two figures are developed in Power Point and do not require any numerical inputs. Table 3 is the data dictionary. [Below](https://github.com/trambakbanerjee/crejm-code/blob/main/README.md#workflow), we provide the steps that must be followed to reproduce the different tables and figures in the paper.

## Supporting software and hardware requirements 
The primary software used is R (>=4.3.0), and the following R packages and their dependencies must be installed for any reproducibility analysis: Rfast (>=2.0.1), msos (>=1.2.0), mvtnorm (>=1.1.1), CVXR (>=1.0.8), igraph (>=1.2.6), [ggb (>=0.10.0)](https://github.com/jacobbien/ggb#ggb), GLMMadaptive (>=0.7.15), lme4 (>=1.1.26), foreach (>=1.5.1), doParallel (>=1.0.16), ggplot2 (>=3.3.2), gridExtra (>=2.3), glmmLasso (>=1.5.1), rpql (>=0.8), readr (>=1.4.0), tidyverse (>=1.3.0), lattice (>=0.20.41), viridisLite (>=0.3.0), reshape2 (>=1.4.4), scales (>=1.1.1), Rcpp (>=1.0.5), RcppArmadillo (>=0.10.1.2.0), RcppProgress (>=0.4.2).

The analyses presented in the paper are based on the following Hardware specifications: Windows 10, 64 bit, with 128GB RAM on an Intel Xeon Gold 6230 CPU. At a minimum, the authors recommend access to 64GB RAM and 10 CPU cores for enabling multi-core parallelization.

## Workflow

We provide a sequence steps for reproducing the different tables and figures in the paper. Expected Run Time for each of these steps will be denoted by ERT.

### Figures 2 and 4
1. Run the script `motivatingfigures.R` to reproduce figures 2 and 4 in the main paper. Please make sure that the `R` working directory is set to `(your folder structure)/crejm-code/data`. `[ERT < 5 minutes]`

### Table 4
1. Run the script `datasummary.R` to reproduce table 4 in the main paper. The script writes out two CSV files that hold the summary statistics in table 4. Before running this script, please make sure that the working directory is set to `<your folder structure>/crejm-code/data`. `[ERT < 5 minutes]`

### Table 1
This is the table that presents the selected fixed / composite effects and the coefficient estimates under the submodels Login Indicator, Duration of Play
and Purchase Propensity. The following steps when executed in the order described below reproduce table 1.

1. Run `dataprocessing.R` to output `out.RData` which is a list that holds the processed training and prediction data. Before running this script, please make sure that the working directory is set to `<your folder structure>/crejm-code/data`. `[ERT < 10 minutes]`
2. Run `crejm_estimation.R` available inside the folder `selection`. This script writes out a list of initial estimates `init.est.RData`. `[ERT ~ 6 hours]`
3. Run `crejm_selection.R` available inside the folder `selection`. This script writes out a list of selected fixed / composite effect predictors in `selection.RData`. `[ERT ~ 4 hours]`
4. To get the coefficient estimates, run `crejm_postselectionestimation.R` available inside the folder `selection`. This script writes out a list `postselection.est.RData` that stores the coefficient estimates, and the estimated covariance matrices of the player and guild specific random effects. `[ERT ~ 6 hours]`  

Please be advised that steps 2, 3 and 4 in the sequence above are both memory and time intensive processes.
  
### Figure 5
1. Figure 5 relies on the output `postselection.est.RData` that is obtained from step (4) above. To reproduce figure 5, run `crejm_randomeffect_network.R` available inside the folder `selection`. `[ERT < 5 minutes]`

### Table 2
Table 2 presents results related to the predictive performance of CREJM and Benchmarks I, II. The following steps when executed in the order described below reproduce table 2.

1. Run `glmmlasso_prediction.R` to output the False Positive (FP) rates, False Negative (FN) rates and Prediction Errors (PE) for Benchmark I.  Please make sure that the working directory and the path to `out.RData` (from step 1 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) are properly set in lines 5 and 8, respectively, of this script. `[ERT < 20 minutes]`
2. Run `rpql_prediction.R` to output the False Positive (FP) rates, False Negative (FN) rates and Prediction Errors (PE) for Benchmark II.  Please make sure that the working directory and the path to `out.RData` (from step 1 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) are properly set in lines 5 and 8, respectively, of this script. Also, please set the appropriate directory in line 83 for sourcing the R scripts in the folder `library`. `[ERT < 30 minutes]`
3. Run `crejm_prediction.R` to output the False Positive (FP) rates, False Negative (FN) rates and Prediction Errors (PE) for CREJM. Please make sure that:
    1. On lines 8 and 9 - the working directory and the path to `selection.RData` (from step 3 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) are properly set.
    2. On lines 11 and 12 - the appropriate directory for sourcing the R scripts in the folder `library` is set.
    3. On line 22 - the appropriate directory for reading `postselection.est.RData` (from step 4 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) is set.
    4. On line 25 - the path to `out.RData` (from step 1 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) is properly set. 
  
 This script uses R packages `foreach` and `doParallel` and expects access to a local cluster of size 10 on line 88. `[ERT < 30 minutes]`    

### Figure 6

1. Run `crejm_guildrandeffs.R` to reproduce figure 6. Please make sure that:
    1. On lines 8 and 9 - the working directory and the path to `selection.RData` (from step 3 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) are properly set.
    2. On lines 11 and 12 - the appropriate directory for sourcing the R scripts in the folder `library` is set.
    3. On line 22 - the appropriate directory for reading `postselection.est.RData` (from step 4 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) is set.
    4. On line 25 - the path to `out.RData` (from step 1 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) is properly set. 
  
 This script uses R packages `foreach` and `doParallel` and expects access to a local cluster of size 10 on line 83. `[ERT < 30 minutes]`    

### Figures 7 and 8
Figure 7 presents three heat-maps, one for each of the three responses, that plot the mean predicted correlation over time of all players that are members of guild k where k = 1,..,50. Please be advised that reproducing figure 7 is both memory and time intensive.

1. Run `crejm_guildcorrelations.R` to reproduce figures 7 and 8. Please make sure that:
    1. On lines 8 and 9 - the working directory and the path to `selection.RData` (from step 3 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) are properly set.
    2. On lines 11 and 12 - the appropriate directory for sourcing the R scripts in the folder `library` is set.
    3. On line 22 - the appropriate directory for reading `postselection.est.RData` (from step 4 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) is set.
    4. On line 25 - the path to `out.RData` (from step 1 of [Table 1](https://github.com/trambakbanerjee/crejm-code#table-1)) is properly set. 

This script uses R packages `foreach` and `doParallel` and expects access to a local cluster of size 20 on lines 87 and 155. `[ERT ~ 48 hours]`   


# References

1. Bien, J., and Tibshirani, R. (2011), "Sparse Estimation of a Covariance Matrix," Biometrika. 98(4). 807–820
