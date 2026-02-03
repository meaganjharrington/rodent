#' Build deterministic prior object
#' @keywords internal

build_prior <- function() {

## Priors (log-scale)
#
logprior <- function(theta) { # Defines log prior density.
  log_b <- theta[seq_len(K)]; log_I <- theta[K + 1L]
  p_b_ind <- sum(stats::dnorm(log_b, mean = mean_beta, sd = sd_beta, log = TRUE)) #Independent priors on log β blocks
  p_b_rw  <- if (K > 1) sum(stats::dnorm(diff(log_b), mean = 0, sd = rw_sd_beta, log = TRUE)) else 0 # random-walk smoothing prior between adjacent β blocks (i.e. don't accept wild jumps between blocks)
  p_I     <- stats::dnorm(log_I, mean = mean_I0, sd = sd_I0, log = TRUE)
  p_b_ind + p_b_rw + p_I
}

## Posterior
density_theta <- function(theta)
  loglikelihood(theta) + logprior(theta) # posterior log density used by monty.

## Return objects needed downstream
list(
  logprior = logprior,
  density_theta = density_theta
)
}
