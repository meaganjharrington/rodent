#' Create the log-likelihood closure
#'
#' @param filter         dust2 filter object
#' @param N, gamma       fixed scalars
#' @param K              number of beta blocks
#' @param time_vec,time0 vectors/scalar for the observation grid
#' @param starts,ends    integer vectors (block boundaries)
#' @param timepoints     length(time_vec)
#' @param map_blocks_exp function(log_b, starts, ends, timepoints) -> beta(t)
#' @param dt             timestep (default 1)
#' @return function(theta) -> scalar log-likelihood
#' @keywords internal
calculate_loglikelihood <- function(
    filter, N, gamma, K,
    time_vec, time0,
    starts, ends, timepoints,
    map_blocks_exp,
    dt = 1
) {
  force(filter); force(N); force(gamma); force(K)
  force(time_vec); force(time0); force(starts); force(ends); force(timepoints)
  force(map_blocks_exp); force(dt)

  loglikelihood <- function(theta) {
    # Unpack theta = [log_beta_1..log_beta_K, log_I0]
    log_I <- theta[K + 1L]
    I0    <- exp(log_I)
    if (!is.finite(I0) || I0 <= 0) return(-Inf)

    beta_series <- likelihood_beta_from_theta(
      theta = theta, K = K,
      starts = starts, ends = ends,
      timepoints = timepoints,
      map_blocks_exp = map_blocks_exp
    )

    pars <- likelihood_pack_pars(
      I0 = I0, N = N, gamma = gamma,
      time_vec = time_vec, time0 = time0,
      beta_series = beta_series, dt = dt
    )

    # dust2 returns a vector of log-likelihood contributions; sum them
    sum(dust2::dust_likelihood_run(filter, pars))
  }

  loglikelihood
}
