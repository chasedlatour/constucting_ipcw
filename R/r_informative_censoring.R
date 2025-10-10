# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROGRAM: r_informative-censoring.R
# PURPOSE: Illustrate how to construct inverse probability of censoring to 
#          address informative censoring conditional on measured variables.
# PROGRAMMER: initial draft by Paul Zivich. Modified by Chase/Catie.
#
# VARIABLE KEY -
#   id:         identifier for each participant
#   t:          last observed follow-up time. Includes simulated loss-to-follow-up
#   ltfu:       simulated indicator of whether the participant was lost-to-follow-up
#   delta:      event indicator. 1=AIDs or death, 0=no event
#   z:          Simulated measured variable related to loss to follow-up and the event
#   t_true:     observed follow-up time without simulated loss-to-follow-up
#   delta_true: event indicator had no loss-to-follow-up occurred
#
# LAST EDIT: Chase Latour 2025/08/14
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up environment ------------------------------------------------------

## Use the renv package to manage R packages in the project ----
## If not yet installed, install the renv package on your machine
# install.packages("renv")

## Install all packages recorded in the renv lockfile. 

renv::restore()

## Package libraries ----
library(survival)
library(dplyr)
library(ggplot2)

## Global parameters ----
end_of_follow_up <- 365  # end of follow-up time for risk estimation

## Read in data ----
actg_cc <- read.csv("data/actg320_sim_censor.csv", header=TRUE, sep=',')


# 0: Prepare and describe cohort ------------------------------------------

## Indicator if LTFU prior to end_of_follow-up
actg_cc$ltfu_efup <- with(actg_cc, as.numeric(ltfu == 1 & t < end_of_follow_up))

## Investigate number censored due to LTFU prior to end_of_follow-up ----

### Number (Percent) LTFU ----
sprintf("%s (%3.1f%%) were LTFU", 
        table(actg_cc$ltfu_efup)[2],
        100*prop.table(table(actg_cc$ltfu_efup))[2])


### LTFU Stratified by Z ----

rprop <- prop.table(table(actg_cc$z, actg_cc$ltfu_efup), margin = 1)
sprintf("Among people with CD4 < 68.5 cells/mm, %3.1f%% were LTFU",
        100*rprop[1,2])
sprintf("Among people with CD4 >=68.5 cells/mm, %3.1f%% were LTFU",
        100*rprop[2,2])



# 1: Calculate true risk functions ignoring simulated LTFU ---------------

## KM estimator for Risk function had no censoring had occurred ----
kmt <- survfit(Surv(t_true, delta_true) ~ 1, data=actg_cc)
true <- 1 - kmt$surv



# 2: Construct IPCW using pooled logistic model ---------------------------

## STEP 1: Create long dataset ----
### From the person-level dataset, create a “long” dataset where each row 
###   corresponds to one-unit (e.g., 1 day, week or month) of follow-up time for 
###   each person 

### 1a. Create a placeholder for the full observed event or censoring time ----
actg_cc$t_obs <- actg_cc$t

### 1b. Converting to long data set with 1-day intervals of follow-up ----
actg_long <- survSplit(Surv(t, delta)~., 
                       data=actg_cc, 
                       # Variable with event indicator
                       event="delta", 
                       # Create one day intervals from first time to end of FUP
                       cut=seq(min(actg_cc$t), end_of_follow_up))

actg_long$tstart_f <- as.factor(actg_long$tstart)
### 1c. Create time-varying censor indicator for each 1-day interval ----
actg_long$not_censored <- with(actg_long, as.numeric(ltfu!=1 | (t_obs != t)))

### Inspect resulting data
head(actg_long)


## STEP 2: Create dataset that excludes events --------
### For the pooled logistic regression censoring model need dataset that excludes 
###   intervals where events occur. 

### Remove those intervals where an event occurs ----
actg_long_2 <- subset(actg_long, delta == 0)



## STEP 3: Estimate interval specific censoring probabilities ----
actg_long_3 <- actg_long_2 |> 
  # Remove intervals where everyone is uncensored (not needed for model fitting)  
  dplyr::filter(!all(not_censored == 1), .by = tstart)

### Fit pooled logistic model ----
logit_den <- glm(not_censored ~ z*tstart_f, 
                 data=actg_long_3, family='binomial')

### Predicted probabilities of remaining uncensored in each time interval ----
actg_long_3$pr_c_den <- predict(logit_den, 
                              type='response') # predicted probabilities


## STEP 4: Calculate cumulative probability of remaining uncensored ----

actg_unstab <- actg_long |> 
  ## Merge the predicted probabilities onto the original dataset with 
  ##   all event times
  dplyr::left_join(actg_long_3 |> 
                     dplyr::select(id, tstart, pr_c_den),
                   by = dplyr::join_by(id, tstart)) |> 
  dplyr::mutate(
    ## Replace missing predicted probabilities with 1
    pr_c_den = replace(pr_c_den, is.na(pr_c_den), 1),
    ## Cumulative probability by person ID
    cum_pr_c_den = cumprod(pr_c_den), default = 1,
    .by = id
  )

## STEP 5: Calculate interval specific weights ----
actg_unstab$ipcw <- 1 / actg_unstab$cum_pr_c_den  # unstabilized IPCW




# 3: Estimate Kaplan-Meier functions -----------------------------------------

## IPC-weighted Kaplan-Meier ----
kmw <- survfit(Surv(tstart, t, delta) ~ 1, 
               data=actg_unstab, 
               weights = ipcw)

ipcw <- 1-kmw$surv

## Naive Kaplan-Meier on "observed" data ----
kmo <- survfit(Surv(t, delta) ~ 1, data=actg_cc)
obs <- 1 - kmo$surv



# Plot naive, true, and IPCW risk functions -------------------------------

dplyr::bind_rows(
  data.frame(time=kmt$time, risk=true, type="Truth"),
  data.frame(time=kmo$time, risk=obs, type="Naive"),
  data.frame(time=kmw$time, risk=ipcw, type="IPCW")
) |> 
  dplyr::mutate(
    type = factor(type, levels = c("Truth", "Naive", "IPCW"))
  ) |> 
  ggplot(aes(x = time, y = risk, group = type, linetype = type,
             color = type)) + 
  theme_classic() + 
  geom_step() + 
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365)) + 
  scale_y_continuous(limits = c(0, 0.15)) + 
  xlab("Days") +
  ylab("Cumulative incidence")






# ALTERNATIVE: Stabilized IPCW ----

### Estimate numerator used pooled logit ----

actg_long_3$pr_c_num <- predict(
  glm(not_censored ~ tstart_f, 
      data=actg_long_3, 
      family='binomial'),
  type = "response"
)

### Calculate stabilized weights ----
actg_stab <- actg_long |> 
  ## Merge the predicted probabilities onto the original dataset with 
  ##   all event times
  dplyr::left_join(actg_long_3 |> 
                     dplyr::select(id, tstart, pr_c_den, pr_c_num),
                   by = dplyr::join_by(id, tstart)) |> 
  dplyr::mutate(
    ## Replace missing predicted probabilities with 1
    pr_c_den = replace(pr_c_den, is.na(pr_c_den), 1),
    pr_c_num = replace(pr_c_num, is.na(pr_c_num), 1),
    
    ### Cumulative probability by person ID ----
    cum_pr_c_den = cumprod(pr_c_den),
    cum_pr_c_num = cumprod(pr_c_num),
    
    ### Stabilized IPCW ----
    icpw = cum_pr_c_num/cum_pr_c_den,
    
    .by = id
  )

