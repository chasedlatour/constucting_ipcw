# constructing_ipcw

This repository contains SAS and R code for the manuscript [xxxx] in which we describe and illustrate how to construct inverse probability of censoring weights to address informative right-censoring due to loss to follow-up. This repository was created to interface with both R and SAS. This repository can be cloned to your computer; and you may interface with the local repo through either software.

A CSV file with the data for this project is stored in the `data/` folder.

R code is saved in the `R/` folder and consists of 3 primary files:
- `ipcw_example.R` -- R script where the analysis is done. Calls in the other two files.
- `helper functions.R`, `ipcw_functions.R` -- R scripts with functions, called in by `ipcw_example.R` to conduct the analysis.

SAS code is saved in the `SAS/` folder and contains only 1 file where all analytic code is contained.
