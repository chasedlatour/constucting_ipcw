# PROGRAM:    helper_functions.R                                                   
# PURPOSE:    Functions to help create data for the IPCW example. Includes
#             function for cutting off follow-up and deriving splines
# PROGRAMMER: C Wiener

#' Cut-off follow-up and rederived censoring/event indicators based on new
#' cut-off time. 
#' @description Limits Lau data to 3 years (can change)
#'
#' @param end_of_fup when to end follow-up (in years)
#' @param path path to the lau dataset used in the example
#' 
cutoff_followup <- function(end_of_fup, path = "data/lau.csv"){
  
  lau <- read.csv(path) |> 
    dplyr::as_tibble()
  
  lau |> 
    dplyr::mutate(
      dthev = dplyr::if_else(t > 3, 0, dthev),
      artev = dplyr::if_else(t > 3, 0, artev),
      t = pmin(t, 3)*365.25 # convert time to days
    )
}

#' Quadratic restricted splines
#' @description Derived quadratic restricted splines based on user-supplied knots.
#'
#' @param x vector of values to create spline basis for
#' @param knots where the knots of the spline should be
#' 
#' @returns data frame of the original vector and the splines. 

qrspline <- function(x, 
                     knots = quantile(x, probs = c(0.25, 0.5, 0.75))){
  
  upper_tail_restriction <- ifelse(x > max(knots), (x - max(knots))^2, 0)
  # create restricted quadratic splines
  
  X <- matrix(nrow = length(x), ncol = length(knots))
  X[,1] <- x
  
  for (i in seq_along(knots[-1])){
    X[,i+1] <- ifelse(x > knots[i], 
                      (x - knots[i])^2 - upper_tail_restriction,
                      0)
  }
  
  colnames(X) <- c("x", paste0('qrs', knots[-length(knots)]))
  
  as.data.frame(X)
}

