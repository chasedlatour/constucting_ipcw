# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROGRAM: create_trial_dataset.R
# PROGRAMMER: Chase Latour
# PURPOSE: To create a version of the ACTG 320 dataset with informative
# censoring due to loss to follow-up (LTFU).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up environment ------------------------------------------------------

## Use the renv package to manage R packages in the project
## If not yet installed, install the renv package on your machine
# install.packages("renv")

## Now, install all packages recorded in the renv lockfile. 
renv::restore()

## Used libraries ----
library(data.table)
library(survival)
library(ggplot2)

## Set the seed for the data generation
set.seed(2394587)

## Read in the ACTG 320 trial data and add column names ----
actg <- fread("data/actg320.23nov16.dat") |> 
  setnames(c("id", "male", "black", "hispanic", "idu", "art", "delta", "drop", 
             "r", "age", "karnof", "days", "cd4", "stop"))

## Subset to 548 people in control arm with complete follow-up ----
### This will serve as the population in which we aim to estimate
### risks with IPCW (actg_cc = actg complete case)
actg_cc <- actg[art == 0 & drop == 0]




# Dichotomize CD4 count ---------------------------------------------------

## We are interested in dichotomizing CD4 cell count to create the variable `Z`. 
## CD4 cell count has a strong relationship with the outcome. 

## First, inspect the distribution ----
hist(actg_cc$cd4)
summary(actg_cc$cd4)

# The median is 68.5 -- discretize at 68. 

## Dichotomize CD4 count at 68.5 ----
actg_cc[, cutoff := ifelse(cd4 < 68.5, 0, 1)]

## Confirm that we see different risk of the outcome across the cutoff strata for CD4
## Risk difference (simple linear model) with Wald Confidence Intervals ----

lmodel <- lm(delta ~ cutoff, data = actg_cc)

sprintf("Risk difference of %3.2f (95%% CI: %3.2f, %3.2f)", 
        coef(lmodel)[2], confint(lmodel)[2,1], confint(lmodel)[2,2])

### See MUCH lower risk (RD = -15%) among those whose CD4 is 68.5 cells/mm3 or greater



# Induce differential censoring -------------------------------------------

## According to the value of cd4 cutoff
## Create new censoring time for each person ----
### where the average time to censoring is differential by CD4.
### We will use a Weibull distribution to create the censoring times.
actg_cc[, days_c := rweibull(.N, shape = 1.5 + 2.5 * cutoff, 
                             scale = ifelse(cutoff == 1, 600, 850))]

## Rename true value variables ----
setnames(actg_cc, c("cutoff", "days", "delta"), c("z", "t_true", "delta_true")) 

## Create new time to event, outcome, and LTFU times ----
### New time to outcome or censoring
actg_cc[, t := as.integer(ifelse(days_c < t_true, days_c, t_true))] 
### New event indicator
actg_cc[, delta := ifelse((days_c >= t_true) & (delta_true == 1), 1, 0)]
### New LTFU indicator
actg_cc[, ltfu := ifelse(days_c < t_true, 1, 0)] 


# Descriptive statistics for censored -------------------------------

mprop <- prop.table(table(actg_cc$ltfu))
sprintf("%3.1f%% of people were LTFU", 100*mprop[2])

## Contingency table for CD4 cutoff and LTFU ----
tab <- table(actg_cc$z, actg_cc$ltfu)
rownames(tab) <- c("<68.5", "68.5+")
colnames(tab) <- c("Not LTFU", "LTFU")

## Proportions (row percents of LTFU by CD4 count) ----
rprop <- prop.table(tab, margin = 1)
sprintf("Among people with CD4 < 68.5 cells/mm, %3.1f%% were LTFU",
        100*rprop[1,2])
sprintf("Among people with CD4 >=68.5 cells/mm, %3.1f%% were LTFU",
        100*rprop[2,2])


# Output new csv ----------------------------------------------------------

fwrite(actg_cc[, .(id, t, ltfu, delta, z, t_true, delta_true)],
       "data/actg320_sim_censor.csv")

