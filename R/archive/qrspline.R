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

