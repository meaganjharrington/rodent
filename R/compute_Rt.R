#' Build deterministic Rt(t) posterior summary object
#' @keywords internal
compute_Rt <- function() {

  # Rt(t) posterior median and 95% CI
  Rt_median <- (beta_median_series / gamma) * (S_t / N)
  Rt_lower  <- (beta_lower_series  / gamma) * (S_t / N)
  Rt_upper  <- (beta_higher_series / gamma) * (S_t / N)

  ## Return objects needed downstream
  list(
    Rt_median = Rt_median,
    Rt_lower  = Rt_lower,
    Rt_upper  = Rt_upper
  )
}
