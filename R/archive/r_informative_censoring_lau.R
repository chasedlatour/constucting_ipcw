#' Convert data to long -- data preparation for IPCW
#'
#' @param df         dataframe containing time to event data
#' @param event_var  event indicator variable name
#' @param censor_var censor indicator variable name
#' @param end_of_fup final day of follow-up to estimate weights
#' @param interval_length how to coarsen time
#'
#' @returns dataframe with tstart and t variables for intervals and time-varying
#'          event and censoring indicators. 

convert_to_long <- function(df, 
                            model_form, 
                            event_var, 
                            censor_var,
                            end_of_fup = 1099,
                            interval_length = 1){
  
  df <- dplyr::mutate(
    df, 
    # Rename event and censoring variables 
    event_var = {{event_var}},
    censor_var = {{censor_var}},
    # Maintain original t before converting to long
    t_obs = t
  )

  # Convert df to long for pooled logistic regression
  df_long <- survSplit(Surv(t, event_var)~.,
                       data=df, 
                       # Variable with event indicator
                       event="delta", 
                       # Create one day intervals from first time to end of FUP
                       cut=seq(0, end_of_fup, by = interval_length))
  
  df_long$t_f <- as.factor(df_long$tstart)
  
  df_long$not_censored <- with(df_long, as.numeric(!(censor_var == 1 & (t_obs == t))))
  
  df_long
  
}

#' IPCW - inverse probability of censoring weights function
#'
#' @param df_long    dataframe containing a row for each time interval and time-
#'                   varying event and censoring indicators. 
#' @param model_form character string indicating censoring model form
#' 
#' @returns dataframe with variable named ipcw with censoring weights

ipcw <- function(df_long, 
                 model_form){
  
  df_long_no_events <- subset(df_long, delta == 0)
  
  if (grepl("t_f*", model_form)){
    df_long_no_events <- df_long_no_events |> 
      # Remove intervals where everyone is uncensored (not needed for model fitting)  
      dplyr::filter(!all(not_censored == 1), .by = tstart)
  }
  
  logit_den <- glm(as.formula(paste0("not_censored ~ ", model_form)),
                   data=df_long_no_events, family='binomial')
  
  ### Predicted probabilities of remaining uncensored in each time interval ----
  df_long_no_events$pr_c_den <- predict(logit_den, 
                                        type='response') # predicted probabilities
  
  
  ## STEP 4: Calculate cumulative probability of remaining uncensored ----
  
  df_unstab <- df_long |> 
    ## Merge the predicted probabilities onto the original dataset with 
    ##   all event times
    dplyr::left_join(df_long_no_events |> 
                       dplyr::select(id, t_f, pr_c_den),
                     by = dplyr::join_by(id, t_f)) |> 
    dplyr::mutate(
      ## Replace missing predicted probabilities with 1
      pr_c_den = replace(pr_c_den, is.na(pr_c_den), 1),
      ## Cumulative probability by person ID
      cum_pr_c_den = cumprod(pr_c_den),
      lag_cum_pr_c_den = lag(cum_pr_c_den, default = 1),
      .by = id
    )
  
  ## STEP 5: Calculate interval specific weights ----
  df_unstab$ipcw <- 1 / df_unstab$lag_cum_pr_c_den  # unstabilized IPCW
  
  df_unstab
}

#' compare_cumul_inc
#' @description Function to create Kaplan Meier survival functions and plot for 
#' IPCW, Naive, and True cumulative incidence functions
#' @param df traditional tte data with t indicating time to event and delta 
#' indicating whether an event occurred (delta = 1) or censoring (delta = 0)
#' @param df_long dataframe returned by the ipcw function
#' @param event_var name of event indicator variable in df
#' @returns list containing the data frame with the survival functions and the 
#' stacked cumulative incidence plots

compare_cumul_inc <- function(df,
                              df_long, 
                              event_var){
  
  if (!c("t") %in% colnames(df)){
    stop("df must contain column t")
  }
  
  if (any(!c("t", "delta") %in% colnames(df_long))){
    stop("df must contain columns t and delta")
  }
  
  df <- dplyr::mutate(
    df, 
    # Rename event variable
    event_var = {{event_var}}
  )

  ## Fit weighted KM
  kmw <- survfit(Surv(tstart, t, delta) ~ 1, 
                 data=df_long, 
                 weights = ipcw)
  
  ## Fit naive KM
  kmn <- survfit(Surv(t, event_var) ~ 1, 
                 data=df)
  
  ## Stack Naive and Weighted survival functions
  h <- bind_rows(
    dplyr::tibble(
      type = "Naive",
      time = kmn$time,
      risk = 1-kmn$surv
    ),
    dplyr::tibble(
      type = "IPCW",
      time = kmw$time,
      risk = 1-kmw$surv
    )
  )
  
  ## If the truth is in the data, fit KM for truth
  
  if (any(!c("true_t", "true_delta") %in% colnames(df))){
    warning("df missing true_t and true_delta - truth won't be plotted")
  }
  else{
    kmt <- survfit(Surv(true_t, true_delta) ~ 1, 
                   data=df)
    h <- bind_rows(
      dplyr::tibble(
        type = "Truth",
        time = kmt$time,
        risk = 1-kmt$surv
      ),
      h
    )
  }
  
  ## Create grouped plot
  p <- h |> 
    ggplot(aes(x = time, color = type, y = risk)) +
    geom_step()
  
  ## Return survival function and plot
  list(
    cum_inc_fns = h,
    cum_inc_plot = p
  )
  
}