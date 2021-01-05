# Code
We provide the full set of analyses code in this GitHub repository that can be used to replicate all the results in the main paper. The code rely on specific software and hardware requirements that have been described below in detail. Specific instructions are provided in the Reproducibility workflow below on how to replicate the tables and figures in the main paper.

## Description

There are 5 folders in this repository:

1. data - holds the data files that were used for all the analyses in this paper.
2. processing - holds R scripts that process the data and generates figures 2, 4 and table 4 in the paper. 
3. selection - R scripts that reproduce results in table 1 and figure 5 in the paper.
4. prediction - holds R scripts that reproduce table 2 and figures 6,7 and 8 in the paper.
5. library - holds the main R functions.

Please navigate into each of these folders and consult the readme files in there. 

## Supporting software and hardware requirements 
The primary software used is R (>=4.3.0), and the following R packages and their dependencies must be installed for any reproducibility analysis: Rfast (>=2.0.1), msos (>=1.2.0), mvtnorm (>=1.1.1), CVXR (>=1.0.8), igraph (>=1.2.6), ggb (>=0.10.0), GLMMadaptive (>=0.7.15), lme4 (>=1.1.26), foreach (>=1.5.1), doParallel (>=1.0.16), ggplot2 (>=3.3.2), gridExtra (>=2.3), glmmLasso (>=1.5.1), rpql (>=0.8), readr (>=1.4.0), tidyverse (>=1.3.0), lattice (>=0.20.41), viridisLite (>=0.3.0), reshape2 (>=1.4.4), scales (>=1.1.1), Rcpp (>=1.0.5), RcppArmadillo (>=0.10.1.2.0), RcppProgress (>=0.4.2).

The analyses presented in the paper are based on the following Hardware specifications: Windows 10, 64 bit, with 128GB RAM on an Intel Xeon Gold 6230 CPU. At a minimum, the authors recommend access to 64GB RAM and 10 CPU cores for enabling multi-core parallelization .


# Reproducibility workflow


## Scope


## Workflow details


## Expected run-time
