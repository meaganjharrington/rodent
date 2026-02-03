#' Expand block betas to daily beta(t)
#' @keywords internal
likelihood_beta_from_theta <- function(theta, K, starts, ends, timepoints, map_blocks_exp) {
  log_b <- theta[seq_len(K)]
  map_blocks_exp(log_b, starts, ends, timepoints)  # returns daily beta(t)
}

#' Build parameter list for dust2 likelihood run
#' @keywords internal
likelihood_pack_pars <- function(I0, N, gamma, time_vec, time0, beta_series, dt = 1) {
  # Ensure beta is defined at time0 by prepending the first value
  beta_times_aug  <- c(time0, time_vec)
  beta_values_aug <- c(beta_series[1], beta_series)

  list(
    N          = N,
    I0         = I0,
    gamma      = gamma,
    dt         = dt,
    n_beta     = length(beta_times_aug),
    beta_values = beta_values_aug,
    beta_times  = beta_times_aug
  )
}
