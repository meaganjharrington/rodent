#' Build deterministic dust2 likelihood object
#'
#' Wraps the odin2 generator creation, dust2 filter construction,
#' and the log-likelihood function exactly as in the original script.
#'
#' @keywords internal

likelihood <- function(
    time_vec,
    time0,
    cases,
    starts,
    ends,
    timepoints,
    N,
    gamma,
    map_blocks_exp)
  {

## odin2 generator from make_sir_timevary_generator
# generator for SIR with time-varying beta
gen <- make_sir_timevary_generator()

## Build a dust2 deterministic likelihood
# prepare data frame columns for the dust filter
dust_data <- data.frame(time = time_vec, cases = cases)

# Create a dust filter object
filter  <- dust2::dust_filter_create(
  gen,
  data = dust_data,
  time_start = time0,  # strictly less than first data time (dust2 requirement)
  dt = 1,       # daily step
  n_particles = 1 # deterministic, could make stochatic in future
)

## Likelihood using dust2
# using test theta (vector of unknown model parameters)
loglikelihood <- function(theta) {
  log_b <- theta[seq_len(K)]
  log_I <- theta[K + 1L]
  I0    <- exp(log_I); if (!is.finite(I0) || I0 <= 0) return(-Inf)
  beta_series <- map_blocks_exp(log_b, starts, ends, timepoints) # construct beta(t) (?)

  ## Augment beta_times/beta_values to include time_start (copy first beta)
  beta_times_aug  <- c(time0, time_vec)
  beta_values_aug <- c(beta_series[1], beta_series) # augment beta so it is defined at time_start

  pars <- list(
    N = N, I0 = I0, gamma = gamma, dt = 1,
    n_beta = length(beta_times_aug),
    beta_values = beta_values_aug,
    beta_times  = beta_times_aug
  )

  sum(dust2::dust_likelihood_run(filter, pars)) }

## Return objects needed downstream
list(
  gen = gen,
  filter = filter,
  loglik = loglikelihood)
}
