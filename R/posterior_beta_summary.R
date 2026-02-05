#' Summarise posterior beta(t) draws into quantiles
#' @param beta_t   matrix [timepoints x n_draws] from posterior_beta()
#' @param q3       a quantile function that returns c("2.5%","50%","97.5%")
#' @keywords internal
summarise_beta_posterior <- function(beta_t, q3) {
  beta_t_q <- apply(beta_t, 1, q3)  # expect 3 x T matrix

  if (is.null(dim(beta_t_q))) {
    beta_t_q <- rbind(`2.5%` = beta_t_q["2.5%"],
                      `50%`  = beta_t_q["50%"],
                      `97.5%`= beta_t_q["97.5%"])
  }

  beta_median_series <- as.numeric(beta_t_q["50%",  ])
  beta_lower_series  <- as.numeric(beta_t_q["2.5%", ])
  beta_upper_series  <- as.numeric(beta_t_q["97.5%",])

  list(
    beta_t_q            = beta_t_q,
    beta_median_series  = beta_median_series,
    beta_lower_series   = beta_lower_series,
    beta_upper_series   = beta_upper_series
  )
}
