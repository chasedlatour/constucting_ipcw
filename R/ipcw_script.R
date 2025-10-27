library(dplyr)
library(survival)
library(ggplot2)

source("R/create_analytic_dataset_lau.R")
source("R/r_informative_censoring_lau.R")

set.seed(1992)

lau_cc_sim <- simulate_censoring(3)

## Quick inspect of HRs
survival::coxph(Surv(t, dthev) ~ cd4_b, data = lau_cc_sim)
survival::coxph(Surv(true_t, true_delta) ~ cd4_b, data = lau_cc_sim)
survival::coxph(Surv(t, artev) ~ cd4_b, data = lau_cc_sim)

weighted_df <- ipcw(lau_cc_sim, "t_f*cd4_b", event_var = dthev, censor_var = artev)

results <- compare_cumul_inc(lau_cc_sim, weighted_df, dthev)

results$cum_inc_plot
