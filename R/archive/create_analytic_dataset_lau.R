
#' simulate_censoring
#' @description Limits Lau data to 3 years and simulates censoring due to CD4 count
#'
#' @param cd4_cutoff Binary cutoff for CD4 to enforce relationship with censoring. 
#' 
simulate_censoring <- function(cd4_cutoff){
  
  lau <- read.csv("data/lau.csv") |> 
    dplyr::as_tibble()
  
  lau |> 
    dplyr::filter(eventtype == 2 | t >= 3) |> 
    dplyr::mutate(
      cd4_b = cd4nadir <= cd4_cutoff,
      true_delta = if_else(t > 3, 0, dthev),
      true_t = pmin(t,3),
      true_artev = artev,
    ) |> 
    dplyr::mutate(
      art_t = rweibull(nrow(lau_cc), shape = 0.2*cd4_b + 1.3*(1-cd4_b), 
                       scale = (1-cd4_b)*8 + (cd4_b)*5),
      dthev = dplyr::if_else(t > 3 | art_t < t, 0, dthev),
      artev = dplyr::if_else(art_t < t & art_t < 3, 1, 0),
      t = ceiling(pmin(t, art_t, 3)*365.25),
      true_t = true_t*365.25
    )
}



