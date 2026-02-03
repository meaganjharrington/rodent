#' Validate main inputs
#'
#' @param incidence A data.frame with columns `time` (consecutive integers)
#'   and `cases` (non-negative integers).
#' @param N Numeric scalar > 0, population size.
#' @param gamma Numeric scalar > 0, fixed recovery rate.
#' @param mcmc List with optional elements:
#'   `n_steps`, `burnin`, `proposal`, `seed`, `n_rt_draws`.
#' @param priors List with optional elements:
#'   `mean_beta`, `sd_beta`, `rw_sd_beta`, `mean_I0`, `sd_I0`.
#' @param inits List with optional elements: `beta`, `I0`.
#'
#' @return
#' @keywords internal
validate_inputs <- function(incidence, N, gamma, mcmc, priors, inits) {

## Validate main inputs - error messages
# ensure data frame with expected columns
stopifnot(is.data.frame(incidence), all(c("time","cases") %in% names(incidence))) # incidence should include time and case columns
stopifnot(all(incidence$cases %% 1 == 0)) # cases as integer!
cases <- incidence$cases
# time must be consecutive integers, but can start anywhere
time_vec <- as.integer(incidence$time)
if (!all(diff(time_vec) == 1L))
  stop("`incidence$time` must be consecutive integers.") # time series
stopifnot(all(is.finite(time_vec)), all(is.finite(cases)), all(cases >= 0)) # all time and cases are finite (and cases are non-neg)
stopifnot(is.numeric(N) && N > 0, is.numeric(gamma) && gamma > 0) # N and gamme +ve

# total_cases  <- sum(cases)
# attack_rate <- total_cases / N # calculate attack rate as logic check, only accepts >1% or <50%
#
# # current epidemic guardrails with attack rate, calculated from N and total cases in dataset
# if (attack_rate < 0.01)
#   stop("attack rate is expected to be >1%, consider changing attack rate to be greater or N to be smaller to represent a true epidemic")
# if (attack_rate > 0.5)
#   stop("attack rate is expected to be <50%, consider making attack rate smaller or N larger to represent a true epidemic")

timepoints    <- length(cases) # number of time points in data
time1 <- time_vec[1L] # first observed time point
time0 <- time1 - 1 # model start time (dust requirement (?), start time < first obs time)

## Pull MCMC settings from inputs with backup defaults
get <- function(x, n, d) if (!is.null(x[[n]])) x[[n]] else d # helper function, pull MCMC settings, with defaults

# WARNING value for X isn't provided, Y used as default

# MCMC control
n_steps     <- get(mcmc, "n_steps", 6000)
burnin      <- get(mcmc, "burnin", 0.5)
proposal    <- get(mcmc, "proposal", NULL)
seed        <- get(mcmc, "seed", 4)
n_rt_draws  <- get(mcmc, "n_rt_draws", 300)

# Priors (on log scale)
mean_beta   <- get(priors, "mean_beta", log(0.3))
sd_beta     <- get(priors, "sd_beta", 0.5)
rw_sd_beta  <- get(priors, "rw_sd_beta", 0.2)
mean_I0     <- get(priors, "mean_I0", log(10))
sd_I0       <- get(priors, "sd_I0", 0.5)

# Initial state (not on log scale)
init_beta   <- get(inits, "beta", 0.3)
init_I0     <- get(inits, "I0", 10)

I0 <- init_I0
if (I0 >= N) return(-Inf) # sanity check: I0 < N
# if (R >= N) return(-Inf) # sanity check: R < N # add back in once

return(list(
  cases = cases,
  time_vec = time_vec,
  timepoints = timepoints,
  time0 = time0,
  n_steps = n_steps,
  burnin = burnin,
  proposal = proposal,
  seed = seed,
  n_rt_draws = n_rt_draws,
  mean_beta = mean_beta,
  sd_beta = sd_beta,
  rw_sd_beta = rw_sd_beta,
  mean_I0 = mean_I0,
  sd_I0 = sd_I0,
  init_beta = init_beta,
  init_I0 = init_I0
))

}
