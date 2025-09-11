#########################################################################
# PROGRAM: create_trial_dataset.R
# PROGRAMMER: Chase Latour
# PURPOSE: To create a version of the ACTG 320 dataset with informative
# censoring due to loss to follow-up (LTFU).
#########################################################################

## Use the renv package to manage R packages in the project
# If not yet installed, install the renv package on your machine
# install.packages("renv")
# Now, call in the renv library
library(renv)

# Now, install all packages in the renv package
renv::restore()

# Call in libraries that are needed
library(data.table)
library(survival)
library(ggplot2)

# Set the seed for the data generation
set.seed(2394587)

# Load in the ACTG 320 trial data
df <- fread("actg320.23nov16.dat") 
# Add column names
setnames(df, c("id", "male", "black", "hispanic", "idu", "art", "delta", "drop", "r",
               "age", "karnof", "days", "cd4", "stop"))


###################################################################
# Subset to those individuals in the control arm who were not LTFU
# This will serve as the population in which we aim to estimate
# risks with IPCW.

dfs <- df[art == 0 & drop == 0]



###################################################################
# We are interested in dichotomizing CD4 cell count to create
# the variable `Z`. CD4 cell count is a good variable because it 
# has a strong relationship with the outcome. 

# Look at the distribution

hist(dfs$cd4)
summary(dfs$cd4)

# The median is 68.5, so we discretize at that value

# Dichotomize CD4 count as <60 and more than 60

dfs[, cutoff := ifelse(cd4 < 68, 0, 1)]

# Confirm that we see different risk of the outcome across the cutoff strata for CD4
# Risk difference (using riskratio.wald from epitools for example)
library(epitools)
tab <- table(dfs$cutoff, dfs$delta)
rr <- riskratio.wald(tab)  # includes RD in output
print(rr)

# See MUCH lower risk among those where CD4 is greater than or equal to 68




###################################################################################
# Induce differential censoring according to the value of cutoff (dichotomized CD4)

# Create a new value for censoring time for each person, where the average time to censoring is
# differential by CD4.
dfs[, days_c := rweibull(.N, shape = 1.5 + 2.5 * cutoff, scale = ifelse(cutoff == 1, 600, 850))]

# Now create new values for tieme to event, outcome, and drop
dfs[, t := ifelse(days_c < days, days_c, days)] # New time to outcome or censoring
dfs[, d := ifelse((days_c >= days) & (delta == 1), 1, 0)] # New event indicator
dfs[, ltfu := ifelse(days_c < days, 1, 0)] # New LTFU indicator

setnames(dfs, c("cutoff", "days", "delta"), c("z", "t_true", "delta_true")) # Rename the true variables
setnames(dfs, "d", "delta") # Rename the new outcome indicator
dfs[, t := as.integer(t)]



##################################################################################
# Some descriptive statistics on the number censored

# Table (X, Y)
table(dfs$z, dfs$ltfu)
prop.table(table(dfs$z, dfs$ltfu))
prop.table(table(dfs$z, dfs$ltfu), margin = 1)

prop.table(table(dfs$ltfu))


###################################################################################
# Output a CSV with the new data

fwrite(dfs[, .(id, t, ltfu, delta, z, t_true, delta_true)],
       "actg320_sim_censor.csv")

