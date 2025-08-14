# constucting_ipcw

This repository contains SAS and R code for the manuscript [xxxx] in which we describe and illustrate how to contruct inverse probability of censoring weights to address informative right-censoring due to loss to follow-up. This repository was created as an R project, and can easily be cloned in R Studio for replication. We use the `renv` package to manage package versions. The data and SAS code can be downloaded and used in SAS statistical software.

Two datasets are present in the folder and described below.
- `actg320.23nov16.dat` -- Original ACTG 320 data, as downloaded from https://github.com/edwardsjk/semiparametric_gcomp.
- `actg320_sim_censor.csv` -- The analytic dataset containing simulated loss to follow-up and outcome times for individauls in the control arm of the trial. These data were generated in R using the file: `create_analytic_dataset.R`.

We provide two codes for implementing inverse probability of censoring weights in R and SAS:
- `r_informative-censoring.R`
- `sas_informative-censoring.sas`
