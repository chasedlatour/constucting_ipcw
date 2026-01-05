# PROGRAM:    ipcw_script.R                                                   
# PURPOSE:    Illustrate how to construct inverse probability of censoring to 
#             address informative censoring conditional on measured variables.
# PROGRAMMER: CD Latour, C Wiener, PN Zivich
#
# REQUIRES:   Lau data, ipcw_functions.R, helper_functions.R


# Load packages, source scripts, read in data -----------------------------

library(tidyr)
library(dplyr)
library(survival)
library(ggplot2)

source("R/helper_functions.R")
source("R/ipcw_functions.R")

end_time <- 3 # end time in years

# read in lau data and cut-off follow-up at 3 years
lau <- cutoff_followup(end_of_fup = end_time,
                       path = "data/lau.csv")


lau$t <- as.integer(lau$t)
# combine admin censoring and art censoring
lau$cens <- as.numeric(lau$eventtype == 1)

# descriptive statistics on censoring -----------------------

# Get quintiles of cd4 counts for descriptive purposes
quintiles <- quantile(lau$cd4nadir, probs = c(0.2, 0.4, 0.6, 0.8, 1))
quintiles

descrip <- lau %>% 
  mutate(cd4_quintiles = case_when(cd4nadir <= quintiles[[1]] ~ 1,
                                   cd4nadir <= quintiles[[2]] ~ 2,
                                   cd4nadir <= quintiles[[3]] ~ 3,
                                   cd4nadir <= quintiles[[4]] ~ 4,
                                   cd4nadir <= quintiles[[5]] ~ 5)) 

table(descrip$cd4_quintiles, descrip$eventtype)
round(prop.table(table(descrip$cd4_quintiles, descrip$eventtype), margin = 2)*100, 0)


# Create quadratic restricted splines for CD4 count -----------------------

# Create splines for CD4 -- knots at the 33rd and 67th percentiles for CD4 count
cd4_splines <- qrspline(lau$cd4nadir, 
                        knots = quantile(lau$cd4nadir, probs = c(0.33, 0.67)))
cd4_colnames <- paste0("cd4_spline_", seq_len(ncol(cd4_splines)))
colnames(cd4_splines) <- cd4_colnames

# Concatenate CD4 splines onto the lau data
lau_cc <- cbind(lau, cd4_splines)


# Convert data to long ----------------------------------------------------

# Step 1
lau_long <- convert_to_long(lau_cc,
                            event_var = dthev, 
                            censor_var = cens,
                            end_of_fup = ceiling(365.25*end_time),
                            interval_length = 1)


# Create quadratic restricted splines for time ----------------------------

# Create splines for time on the interval data frame using knots from the original
# times - knots at teh 25th, 50th, and 75th percentiles of follow-up

time_splines <- qrspline(lau_long$t, # t from long
                         knots = quantile(lau$t, probs = c(0.25, 0.5, 0.75)))

time_colnames <- paste0("time_spline_", seq_len(ncol(time_splines)))
colnames(time_splines) <- time_colnames

# bind splines with the long dataframe
lau_long_cc <- cbind(lau_long, time_splines)



# Estimate IPCW & add to data frame -------------------------------

# Steps 2 - 5
weighted_df <- ipcw(lau_long_cc,
                    # Modeling time and CD4 as restricted quadratic splines. 
                    # The spline basis for time is interacted with the CD4 spline
                    # basic, an indicator for race, and indicator for baseline 
                    # injection drug use. 
                    model_form = sprintf("(%s)*(%s)", 
                                         paste0(time_colnames, collapse = "+"),
                                         paste0(cd4_colnames, collapse = "+")))

# Step 6 -- Evaluate the distribution of weights at each time point
summary(weighted_df$ipcw)



# Estimate weighted and naive incidence functions and compare -------------

# Step 7
# compare IP cumulative incidence function + naive cumulative incidence. 
results <- compare_cumul_inc(lau, weighted_df, dthev)
results$cum_inc_plot 

results$cum_inc_fns |> 
  dplyr::filter(time >= 365*end_time)

