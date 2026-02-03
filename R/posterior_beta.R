#' Expand posterior beta-block draws onto the time grid
#' @keywords internal
posterior_beta <- function(beta_blocks_samp, starts, ends, timepoints) {
  K        <- nrow(beta_blocks_samp)
  n_draws  <- ncol(beta_blocks_samp)

  beta_t_samp <- matrix(NA_real_, nrow = timepoints, ncol = n_draws)

  for (j in seq_len(n_draws)) {
    out <- numeric(timepoints)
    for (k in seq_len(K)) {
      out[starts[k]:ends[k]] <- beta_blocks_samp[k, j]
    }
    beta_t_samp[, j] <- out
  }

  beta_t_samp
}
