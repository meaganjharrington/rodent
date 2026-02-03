#' Summarise posterior beta(t) draws into quantiles
#' @param beta_t   matrix [timepoints x n_draws] from posterior_beta()
#' @param q3fun    a quantile function that returns c("2.5%","50%","97.5%")
#' @keywords internal
summarise_beta_posterior <- function(beta_t, q3) {
  # apply over time (rows), get a 3 x timepoints matrix
  beta_t_q <- apply(beta_t, 1, q3)

  # ensure orientation is 3 x timepoints (some q3 implementations return named numeric)
  if (is.null(dim(beta_t_q))) {
    beta_t_q <- rbind(`2.5%` = beta_t_q["2.5%"],
                      `50%`  = beta_t_q["50%"],
                      `97.5%`= beta_t_q["97.5%"])
  }

  beta_median_series <- beta_t_q["50%",  ]
  beta_lower_series  <- beta_t_q["2.5%", ]
  beta_upper_series  <- beta_t_q["97.5%", ]

  list(
    beta_t_q            = beta_t_q,
    beta_median_series  = beta_median_series,
    beta_lower_series   = beta_lower_series,
    beta_upper_series   = beta_upper_series
  )
}
