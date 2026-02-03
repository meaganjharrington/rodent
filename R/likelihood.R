#' Build deterministic dust2 likelihood object (orchestrator)
#'
#' Wraps odin2 generator creation, dust2 filter construction,
#' and returns a log-likelihood closure.
#'
#' @keywords internal
likelihood_create <- function(
    time_vec,
    time0,
    cases,
    starts,
    ends,
    timepoints,
    N,
    gamma,
    K,
    map_blocks_exp,
    dt = 1,
    n_particles = 1
) {
  gen    <- likelihood_build_generator()
  filter <- likelihood_build_filter(gen, time_vec, cases, time0, dt, n_particles)

  loglikelihood <- calculate_loglikelihood(
    filter = filter, N = N, gamma = gamma, K = K,
    time_vec = time_vec, time0 = time0,
    starts = starts, ends = ends, timepoints = timepoints,
    map_blocks_exp = map_blocks_exp, dt = dt
  )

  list(
    gen            = gen,
    filter         = filter,
    loglikelihood  = loglikelihood
  )
}
