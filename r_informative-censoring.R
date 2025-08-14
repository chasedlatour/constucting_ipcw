#######################################################################################################################
# Some title here...
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
# Last edit: Paul Zivich 2020/01/30
#######################################################################################################################

library(survival)
library(dplyr)
library(ggplot2)

# Global parameters
end_of_follow_up <- 365
setwd("C:/Users/zivic/Documents/Research Projects/Steve Cole/IPCW/supplement")

# Reading in data
df <- read.csv("actg320_sim-censor.csv", header=TRUE, sep=',')

####################################################
# 0: True risk functions
####################################################

# Data if no censoring had occurred
kmt <- survfit(Surv(t_true, delta_true) ~ 1, data=df)
true <- 1 - kmt$surv

####################################################
# 1: Bounds on informative censoring
####################################################

# Upper bound: all censored have the event at their lost to follow-up time
df$delta_upper <- ifelse(df$ltfu==1, 1, df$delta)
kmu <- survfit(Surv(t, delta_upper) ~ 1, data=df)
upper <- data.frame(r_upper=1 - kmu$surv)
upper$t <- kmu$time

# Lower bound: all lost to follow-up don't have the event for full study period
df$t_lower <- ifelse(df$ltfu==1, end_of_follow_up, df$t)
kml <- survfit(Surv(t_lower, delta) ~ 1, data=df)
lower <- data.frame(r_lower=1 - kml$surv)
lower$t <- kml$time

bounds <- merge(lower, upper, by='t')

####################################################
# 2: IPCW - pooled logistic model
####################################################

# Converting to long data set
df$t_obs <- df$t
dfl <- survSplit(Surv(t, delta)~., data=df, event="delta", cut=seq(1, end_of_follow_up))
dfl$not_censored <- ifelse((dfl$ltfu==1) & (dfl$t_obs == dfl$t), 0, 1)

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
