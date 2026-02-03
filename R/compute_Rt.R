#' Build deterministic Rt(t) posterior summary object
#' @param beta_q  list with elements:
#'        beta_median_series, beta_lower_series, beta_upper_series
#' @param S_t     either a numeric vector of S(t) or a list with element $S_t
#' @param N       population size (scalar > 0)
#' @param gamma   recovery rate (scalar > 0)
#' @keywords internal
compute_Rt <- function(beta_q, S_t, N, gamma) {
  # Accept S_t as list or vector
  S_t_vec <- if (is.list(S_t) && !is.null(S_t$S_t)) S_t$S_t else S_t
  if (!is.numeric(S_t_vec)) {
    stop("compute_Rt: 'S_t' must be a numeric vector or a list with element $S_t (numeric).")
  }

  # Pull beta series from the summary list
  beta_median_series <- beta_q$beta_median_series
  beta_lower_series  <- beta_q$beta_lower_series
  beta_upper_series  <- beta_q$beta_upper_series

  # Checks
  stopifnot(is.numeric(beta_median_series),
            is.numeric(beta_lower_series),
            is.numeric(beta_upper_series))
  if (!all(lengths(list(beta_median_series, beta_lower_series, beta_upper_series)) ==
           length(S_t_vec))) {
    stop("compute_Rt: beta series and S_t must have the same length.")
  }

  # Rt(t) = (beta(t) / gamma) * S(t) / N
  Rt_median <- (beta_median_series / gamma) * (S_t_vec / N)
  Rt_lower  <- (beta_lower_series  / gamma) * (S_t_vec / N)
  Rt_upper  <- (beta_upper_series  / gamma) * (S_t_vec / N)

  list(
    Rt_median = Rt_median,
    Rt_lower  = Rt_lower,
    Rt_upper  = Rt_upper
  )
}
