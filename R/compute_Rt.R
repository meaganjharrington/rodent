#' Compute Rt(t) summaries from beta_q and S(t)
#' @keywords internal
compute_Rt <- function(beta_q, S_t, N, gamma) {
  beta_med <- beta_q$beta_median_series
  beta_lo  <- beta_q$beta_lower_series
  beta_hi  <- beta_q$beta_upper_series

  # Align all lengths (use shortest)
  L <- min(length(beta_med), length(beta_lo), length(beta_hi), length(S_t))
  beta_med <- beta_med[seq_len(L)]
  beta_lo  <- beta_lo[seq_len(L)]
  beta_hi  <- beta_hi[seq_len(L)]
  S_t      <- S_t[seq_len(L)]

  Rt_median <- (beta_med / gamma) * (S_t / N)
  Rt_lower  <- (beta_lo  / gamma) * (S_t / N)
  Rt_upper  <- (beta_hi  / gamma) * (S_t / N)

  list(
    Rt_median = Rt_median,
    Rt_lower  = Rt_lower,
    Rt_upper  = Rt_upper
  )
}
