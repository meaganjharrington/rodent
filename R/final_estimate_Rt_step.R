
#' Rt(t) with stepwise beta(t), fixed gamma; infer I0 and beta blocks via Metropolis-Hastings random walk MCMC
#' Using full odin2 -> dust2 -> monty workflow (deterministic likelihood)
#'
#' @param incidence data.frame with columns: time (consecutive integers), cases >= 0.
#'        Time can start at any integer; we keep it and set time_start to "first time - 1".
#' @param N Numeric scalar > 0, population size
#' @param gamma Numeric scalar > 0, fixed recovery rate
#' @param R Immunity existing in population, default = 0
#' @param beta_breaks Integer vector in ORIGINAL time units; mapped to indices internally
#' @param mcmc List: n_steps (6000), burnin (0.5), proposal (NULL -> diag(0.02^2)), seed (4),
#'              n_rt_draws (300)
#' @param priors List: mean_beta=log(0.3), sd_beta=0.5, rw_sd_beta=0.2,
#'               mean_I0=log(10), sd_I0=0.5
#' @param inits  List: beta=0.3, I0=10
#' @return list(samples, estimates, Rt_series, fit, diagnostics, model_used, blocks, fixed)
#' @export
final_estimate_Rt_step <- function(
    incidence,
    N,
    gamma,
    # R,
    beta_breaks = NULL,
    mcmc  = list(n_steps = 6000, burnin = 0.5, proposal = NULL, seed = 4, n_rt_draws = 300),
    # MCMC control settings (definitions): n_steps = , burnin = , proposal = , seed = , n_rt_draws =
    priors = list(mean_beta = log(0.3), sd_beta = 0.5, rw_sd_beta = 0.2,
                  mean_I0 = log(10), sd_I0 = 0.5), # priors, defined on log scale
    inits  = list(beta = 0.3, I0 = 10) # initial values for MCMC
) {

  ## Validate main inputs
  stopifnot(is.data.frame(incidence), all(c("time","cases") %in% names(incidence))) # incidence should include time and case columns
  stopifnot(all(incidence$cases %% 1 == 0)) # cases as integer!
  cases <- incidence$cases
  time_vec <- as.integer(incidence$time)
  if (!all(diff(time_vec) == 1L))
    stop("`incidence$time` must be consecutive integers.") # time series
  stopifnot(all(is.finite(time_vec)), all(is.finite(cases)), all(cases >= 0)) # all time and cases are finite (and cases are non-neg)
  stopifnot(is.numeric(N) && N > 0, is.numeric(gamma) && gamma > 0) # N and gamme +ve

  timepoints    <- length(cases) # number of time points in data
  time1 <- time_vec[1L] # first observed time point
  time0 <- time1 - 1 # model start time (dust requirement (?), start time < first obs time)(so we prepend the first beta value)

  ## Pull MCMC settings with defaults
  get <- function(x, n, d) if (!is.null(x[[n]])) x[[n]] else d # helper function, pull MCMC settings, with defaults
  n_steps     <- get(mcmc, "n_steps", 6000)
  burnin      <- get(mcmc, "burnin", 0.5)
  proposal    <- get(mcmc, "proposal", NULL)
  seed        <- get(mcmc, "seed", 4)
  n_rt_draws  <- get(mcmc, "n_rt_draws", 300)

  mean_beta   <- get(priors, "mean_beta", log(0.3))
  sd_beta     <- get(priors, "sd_beta", 0.5)
  rw_sd_beta  <- get(priors, "rw_sd_beta", 0.2)
  mean_I0     <- get(priors, "mean_I0", log(10))
  sd_I0       <- get(priors, "sd_I0", 0.5)

  init_beta   <- get(inits, "beta", 0.3)
  init_I0     <- get(inits, "I0", 10)

  I0 <- init_I0
  if (I0 >= N) return(-Inf) # sanity check: I0 < N
  # if (R >= N) return(-Inf) # sanity check: R < N

  ## Beta blocks - time periods with constant transmission rate (piecewise constant)
  ## beta_breaks are supplied in ORIGINAL time units; map to 1-based positions
  make_breaks <- function(breaks, time_vec) { # convert calender times into array indicies
    if (!is.null(breaks)) {
      breaks <- sort(unique(as.integer(breaks)))
      indices_of_betabreaks <- match(breaks, time_vec) # indicies of beta breaks within time vector (?)
      if (any(is.na(indices_of_betabreaks)))
        stop("`beta_breaks` must align to observed times in `incidence$time`.")
      if (indices_of_betabreaks[1] != 1L)
        stop("First beta block must start at the first observation.")
      return(indices_of_betabreaks)
    }
    1L # if no breaks provided, single beta block starting at 1
  }

  starts <- make_breaks(beta_breaks, time_vec)
  ends   <- c(starts[-1] - 1L, timepoints) # block k ends just before k+1, etc
  K      <- length(starts) # number of beta blocks

  ## Parameter vector on log-scale: [log_beta_1....log_beta_K, log_I0]
  par_names <- c(paste0("log_beta_", seq_len(K)), "log_I0") # parameters in MCMC state vector
  init      <- c(rep(log(init_beta), K), log(init_I0)) # initial MCMC position
  if (is.null(proposal)) proposal <- diag(rep(0.05^2, K + 1))

  ## Helpers
  map_blocks_exp <- function(vals_log, starts, ends, timepoints) {
    vals <- exp(vals_log);
    out <- numeric(timepoints)
    for (k in seq_along(starts)) out[starts[k]:ends[k]] <- vals[k]
    out
  }
  q3 <- function(x) stats::quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE)

  ## odin2 generator
  gen <- make_sir_timevary_generator()

  ## Build a dust2 deterministic likelihood
  dust_data <- data.frame(time = time_vec, cases = cases)
  unfilter  <- dust2::dust_unfilter_create(
    gen,
    data = dust_data,
    time_start = time0,  # strictly less than first data time
    dt = 1
  )

  ## Likelihood using dust2
  loglik <- function(theta) {
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

    sum(dust2::dust_likelihood_run(unfilter, pars))
  }

  ## Priors (log-scale)
  logprior <- function(theta) { # Defines log prior density.
    log_b <- theta[seq_len(K)]; log_I <- theta[K + 1L]
    p_b_ind <- sum(stats::dnorm(log_b, mean = mean_beta, sd = sd_beta, log = TRUE)) #Independent priors on log β blocks
    p_b_rw  <- if (K > 1) sum(stats::dnorm(diff(log_b), mean = 0, sd = rw_sd_beta, log = TRUE)) else 0 # random-walk smoothing prior between adjacent β blocks.
    p_I     <- stats::dnorm(log_I, mean = mean_I0, sd = sd_I0, log = TRUE)
    p_b_ind + p_b_rw + p_I
  }

  ## Posterior
  density_theta <- function(theta)
    loglik(theta) + logprior(theta) # posterior log density used by monty.

  ## MCMC via monty (random-walk MH)
  mod <- monty::monty_model(list(density = density_theta, parameters = par_names)) # monty model object (w/ density and pars)(density = log posterior density to sample from)
  set.seed(seed) # fixed random number gen
  sampler <- monty::monty_sampler_random_walk(proposal) # cteaye rw MH sampler
  n_burn  <- floor(burnin * n_steps) # compute burn-in no (burnin 0.5 = first 1/2 of samples dropped)(?)
  smp     <- monty::monty_sample(mod, sampler, n_steps, initial = init, n_chains = 1) # runs MCMC!
  # runs n_steps iterations, produces monty_samples w/ parameter values each time

  ## Extract post-burnin samples
  pars <- smp$pars; # extract raw MCMC parameter array (dimensions vary by chain no)
    dims <- dim(pars); # stores array dimensions
      indices <- (n_burn + 1):dims[2] # indicies of post-burnin iterations
  extract_post <- if (length(dims) == 3) as.matrix(pars[, indices, 1, drop = FALSE]) else as.matrix(pars[, indices, drop = FALSE])
  # extract post burn-in samples (3D parameters x iterations x chains AND 2D parameters x iterations)
  # After extraction
  dims <- dim(pars)
  indices <- (n_burn + 1):dims[2]

  extract_post <- if (length(dims) == 3) {
    matrix(pars[, indices, 1], nrow = dims[1], ncol = length(indices))
  } else {
    as.matrix(pars[, indices, drop = FALSE])
  }

  # Extract beta blocks and I0 using indices
  beta_blocks_samp <- exp(extract_post[1:K, , drop = FALSE])
  I0_samp          <- exp(extract_post[K + 1, , drop = FALSE])

  beta_q <- t(apply(beta_blocks_samp, 1, q3)) # posterior summaries for each beta block (? - is this needed)
  I0_q <- q3(I0_samp) # posterior median/95% interval for I0 (? - is this necessary)

  ## Medians and Rt(t)
  beta_med_blocks <- beta_q[, 1] # extract post. median beta for each beta block
  beta_med_series <- map_blocks_exp(log(beta_med_blocks), starts, ends, timepoints) # expand into full time series

  ## Augment beta vectors for system creation
  beta_times_med  <- c(time0, time_vec)
  beta_values_med <- c(beta_med_series[1], beta_med_series)

  sys_med <- dust2::dust_system_create(
    gen,
    pars = list(
      N = N, I0 = I0_q[1], gamma = gamma, dt = 1,
      n_beta = length(beta_times_med),
      beta_values = beta_values_med,
      beta_times  = beta_times_med
    ),
    n_particles = 1, deterministic = TRUE
  ) # generate determinstic system using posterior median parameters on I0 and beta(t)
  dust2::dust_system_set_state_initial(sys_med)
  res_med <- dust2::dust_system_simulate(sys_med, time_vec)  # simulate over data times only

  ## Extract S(t) and incidence
  rn_med <- rownames(res_med)
  if (!is.null(rn_med)) {
    S_row   <- which(rn_med == "S")
    inc_row <- which(rn_med == "incidence")
  } else {
    # fallback if rownames missing
    S_row   <- 1L
    inc_row <- 3L  # incidence is the 3rd row in your generator
  }

  S_t        <- as.numeric(res_med[S_row, ])
  lambda_med <- pmax(as.numeric(res_med[inc_row, ]), 1e-12)
  Rt_med     <- (beta_med_series * S_t) / (gamma * N)

  ## Return
  list(
    samples = data.frame(t(rbind(beta_blocks_samp, I0_samp))), # posterior samples
    estimates = list(
      beta_blocks = data.frame(
        block = seq_len(K), start = starts, end = ends,
        beta_median = beta_q[, 1], beta_lower = beta_q[, 2], beta_upper = beta_q[, 3]
      ),
      I0 = data.frame(I0_median = I0_q[1], I0_lower = I0_q[2], I0_upper = I0_q[3])
    ),
    Rt_series = data.frame(time = time_vec, Rt_median = Rt_med),
    model_used = "SIR_deterministic_step_beta_gamma_fixed_dust2_monty",
    blocks = list(beta_starts = starts, beta_ends = ends),
    fixed = list(N = N, gamma = gamma)
  )
}
