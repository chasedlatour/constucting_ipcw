#######################################################################################################################
# PROGRAM: r_informative-censoring.R
# PURPOSE: Illustrate how to contruct inverse probability of censoring to address informative censoring
# conditional on measured variables.
# PROGRAMMER: initial draft by Paul Zivich. Modified by Chase/Catie.
#
# Variable key -
#   id:         identifier for each participant
#   t:          last observed follow-up time. Includes simulated loss-to-follow-up
#   ltfu:       indicator of whether the participant was lost-to-follow-up. Based on simulated data
#   delta:      event indicator. 1=AIDs or death, 0=no event
#   x:          measured variable related to loss to follow-up and the event. Simulated
#   t_true:     observed follow-up time without simulated loss-to-follow-up. Used to create true risk function
#   delta_true: event indicator had no loss-to-follow-up occurred. Used to create true risk function
#
# Last edit: Chase Latour 2025/08/14
#######################################################################################################################


## Use the renv package to manage R packages in the project
# If not yet installed, install the renv package on your machine
# install.packages("renv")
# Now, call in the renv library
library(renv)

# Now, install all packages in the renv package
renv::restore()

# Call in libraries that are needed
library(survival)
library(dplyr)
library(ggplot2)

# Global parameters
end_of_follow_up <- 365
#setwd("C:/Users/zivic/Documents/Research Projects/Steve Cole/IPCW/supplement")

# Read in data
df <- read.csv("actg320_sim_censor.csv", header=TRUE, sep=',')








####################################################
# 0: Prepare and describe the cohort
####################################################

# Prepare the dataset for analysis
df <- df %>% 
  mutate(
    # Create an indicator for LTFU before 365 days of follow-up
    ltfu_efup = ifelse(ltfu == 1 & t < end_of_follow_up, 1, 0),
    
    # If there are any 0 follow-up times, add a small constant
    t = ifelse(t == 0, 0.0001, t)
  )

# Investigate number censored due to LTFU prior to end_of_follow-up

# Number LTFU 
table(df$ltfu_efup)

# Proportion LTFU
prop.table(table(df$ltfu_efup))

## Same but stratified by Z

# Number LTFU 
table(df$z, df$ltfu_efup)

# Proportion LTFU
prop.table(df$z, table(df$ltfu_efup))












####################################################
# 1: Calculate the true risk functions (i.e.,
# ignoring the simulated LTFU)
####################################################

# Data if no censoring had occurred
kmt <- survfit(Surv(t_true, delta_true) ~ 1, data=df)
true <- 1 - kmt$surv







####################################################
# 2: Construct the IPCW using a pooled logistic 
# model
####################################################

##################
#  STEP 1: From the person-level dataset, create a “long” dataset where each row 
#  corresponds to one-unit (e.g., 1 day, week or month) of follow-up time for each person 

# Create a new variable for the observed event or censoring time
df$t_obs <- df$t

# Calculate the minimum t in the data
min <- min(df$t)

# Converting to long data set using 1-day intervals of follow-up
dfl <- survSplit(Surv(t, delta)~., data=df, event="delta", cut=seq(min, end_of_follow_up))

# Create an indicator for not being censored during that interval
dfl$not_censored <- ifelse((dfl$ltfu==1) & (dfl$t_obs == dfl$t), 0, 1)

##################
#  STEP 2: Create a new dataset for the pooled logistic regression model 
#  that excludes intervals where events occur. 

# Remove those intervals where an event occurs
dfl_2 <- subset(dfl, delta == 0)



##################
#  STEP 2

# Fitting pooled logistic model
logit_den <- glm(not_censored ~ x + tstart + I(tstart^2) + I(tstart^3)+ I(x*tstart) + I(x*tstart^2) + I(x*tstart^3), data=dfl, family='binomial')
dfl$pr_c_den <- predict(logit_den, type='response') # predicted probabilities
dfl <- dfl %>% 
  group_by(id) %>%
  mutate(pr_c_den = cumprod(pr_c_den)) # cumulative product by ID

dfl$ipcw <- 1 / dfl$pr_c_den  # unstabilized IPCW

# STABILIZED IPCW ALTERNATIVE
# logit_num <- glm(not_censored ~ tstart + I(tstart^2) + I(tstart^3), data=dfl, family='binomial')
# dfl$pr_c_num <- predict(logit_num, type='response') # predicted probabilities
# dfl <- dfl %>%
#  group_by(id) %>%
#  mutate(pr_c_num = cumprod(pr_c_num)) # cumulative product by ID
# dfl$ipcw <- dfl$pr_c_num / dfl$pr_c_den  # stabilized IPCW

# Estimating a IPC-weighted Kaplan-Meier
kmw <- survfit(Surv(tstart, t, delta) ~ 1, data=dfl, weights=dfl$ipcw)

####################################################
# 3: naive Kaplan-Meier
####################################################

# Observed data only
kmo <- survfit(Surv(t, delta) ~ 1, data=df)
obs <- 1 - kmo$surv

####################################################
# Creating Plot
####################################################

ggplot() + theme_light() + 
  geom_rect(aes(xmin = bounds$t, xmax = lead(bounds$t), 
                ymin = bounds$r_lower, ymax = bounds$r_upper), 
            fill = "gray", alpha = 0.25) +
  geom_step(aes(x=kmt$time, y=true), color="gray") + 
  geom_step(aes(x=kmo$time, y=obs), color="black") +
  geom_step(aes(x=kmw$time, y=ipcw), color="black", linetype="dashed") +
  xlab("t (days)") +
  ylab("Risk at t")
