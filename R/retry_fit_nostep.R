
#' Fit Constant-Beta SIR Model to Incidence Data (Deterministic Likelihood)
#'
#' Estimates beta (and optionally gamma) via MCMC using a deterministic SIR and
#' a Negative Binomial (default) or Poisson observation model. Returns posterior
#' samples, summaries, R0, R(t) with credible intervals, and posterior predictive incidence.
#'
#' @param incidence_data Data frame with columns 'time' and 'cases'
#' @param N Total population size
#' @param I0 Initial infections (default: inferred from early cases)
#' @param gamma Fixed recovery rate (ignored if estimate_gamma=TRUE)
#' @param estimate_gamma Logical: estimate gamma? (default FALSE)
#' @param beta_prior_mean Prior mean for beta (Normal)
#' @param beta_prior_sd Prior SD for beta (Normal)
#' @param gamma_prior_mean Prior mean for gamma (Normal, if estimated)
#' @param gamma_prior_sd Prior SD for gamma (Normal, if estimated)
#' @param use_negbin Logical: use Negative Binomial instead of Poisson? (default TRUE)
#' @param nb_size Dispersion (size) for NB (fixed; larger = less overdispersion). Default 10.
#' @param rho Reporting rate (fraction observed, default 1.0; fixed)
#' @param n_steps MCMC iterations per chain (default 2000)
#' @param n_chains Number of chains (default 3)
#' @param n_burnin Burn-in iterations (default n_steps/2)
#' @param n_ppc Number of posterior predictive draws (default 200)
#' @return A list of class 'epievolve_fit' with: samples, estimates, R_t, fitted_incidence, mcmc_info, fixed_parameters
retry_fit_nostep <- function(incidence_data,
                           N,
                           I0 = NULL,
                           gamma = 0.1,
                           estimate_gamma = FALSE,
                           beta_prior_mean = 0.3,
                           beta_prior_sd = 0.2,
                           gamma_prior_mean = 0.1,
                           gamma_prior_sd = 0.05,
                           use_negbin = TRUE,
                           nb_size = 10,
                           rho = 1.0,
                           n_steps = 2000,
                           n_chains = 3,
                           n_burnin = NULL,
                           n_ppc = 200) {

  # --- 0) Validate and prepare data ---
  if (!all(c("time", "cases") %in% names(incidence_data))) {
    stop("incidence_data must have 'time' and 'cases' columns")
  }
  cases <- as.numeric(incidence_data$cases)

  # Force times to 0,1,2,... to match dust step indexing
  times <- seq.int(0L, length(cases) - 1L)

  # Infer I0 from early cases if not provided
  if (is.null(I0)) {
    I0 <- max(1, round(mean(cases[seq_len(min(5, length(cases)))])))
  }

  if (is.null(n_burnin)) {
    n_burnin <- floor(n_steps / 2)
  }

  # --- 1) Define the SIR model via odin2 ---
  sir <- odin2::odin({
    deriv(S) <- -beta * S * I / N
    deriv(I) <- beta * S * I / N - gamma * I
    deriv(R) <- gamma * I
    deriv(Inc) <- beta * S * I / N  # cumulative incidence derivative

    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0
    initial(Inc) <- 0

    N <- parameter()
    I0 <- parameter()
    beta <- parameter()
    gamma <- parameter()
  })

  # --- 2) Which parameters to estimate ---
  if (estimate_gamma) {
    packer <- monty::monty_packer(c("beta", "gamma"))
    param_names <- c("beta", "gamma")
  } else {
    packer <- monty::monty_packer("beta", fixed = list(gamma = gamma))
    param_names <- "beta"
  }

  # --- 3) Likelihood (Poisson or NegBin observation model) ---
  likelihood_fn <- function(pars) {
    # Unpack parameters into full list (includes fixed gamma when not estimated)
    p <- packer$unpack(pars)
    beta_val  <- p$beta
    gamma_val <- p$gamma

    # Bounds: strictly positive rates
    if (beta_val <= 0 || gamma_val <= 0 || nb_size <= 0 || rho <= 0 || rho > 1) {
      return(-Inf)
    }

    # Run deterministic SIR
    sys <- dust2::dust_system_create(
      sir,
      list(N = N, I0 = I0, beta = beta_val, gamma = gamma_val),
      n_particles = 1,
      deterministic = TRUE
    )
    dust2::dust_system_set_state_initial(sys)
    result <- dust2::dust_system_simulate(sys, times)

    # States: 1=S, 2=I, 3=R, 4=Inc (cumulative)
    cum_inc <- result[4, ]
    daily_inc_true <- c(cum_inc[1], diff(cum_inc))

    # Reporting fraction
    lambda <- pmax(daily_inc_true * rho, 1e-6)

    # Observation density
    if (use_negbin) {
      # NB parameterization with size and mu; Var = mu + mu^2/size
      ll <- sum(stats::dnbinom(x = cases, size = nb_size, mu = lambda, log = TRUE))
    } else {
      ll <- sum(stats::dpois(x = cases, lambda = lambda, log = TRUE))
    }
    ll
  }

  # --- 4) Wrap the likelihood function as a monty model ---

  likelihood <- monty::monty_model(list(
    density    = likelihood_fn,    # function(pars) -> log-density
    parameters = param_names,      # e.g., "beta" or c("beta","gamma")
    domain     = NULL              # full real line
  ))


  # --- 5) Priors (Normal on beta/gamma) ---

  # --- 5) Priors (Normal on beta/gamma) ---
  if (estimate_gamma) {
    prior <- monty::monty_dsl({
      beta_mu <- !!beta_prior_mean
      beta_sd <- !!beta_prior_sd
      gamma_mu <- !!gamma_prior_mean
      gamma_sd <- !!gamma_prior_sd
      beta  ~ Normal(mean = beta_mu,  sd = beta_sd)
      gamma ~ Normal(mean = gamma_mu, sd = gamma_sd)
    })
  } else {
    prior <- monty::monty_dsl({
      beta_mu <- !!beta_prior_mean
      beta_sd <- !!beta_prior_sd
      beta ~ Normal(mean = beta_mu, sd = beta_sd)
    })
  }


  # Combine into posterior
  posterior <- likelihood + prior

  # --- 6) Sampler setup (random-walk proposal) ---
  vcv <- diag(length(param_names))
  diag(vcv) <- 0.01  # tune if needed
  sampler <- monty::monty_sampler_random_walk(vcv = vcv)

  # Initial vector in parameter order
  initial <- if (estimate_gamma) c(beta_prior_mean, gamma_prior_mean) else beta_prior_mean

  # --- 7) Run MCMC ---
  samples_raw <- monty::monty_sample(
    posterior,
    sampler,
    n_steps,
    initial = initial,
    n_chains = n_chains
  )

  # --- 8) Extract post-burnin samples as vectors ---
  keep <- (n_burnin + 1):n_steps
  arr <- samples_raw$pars[, , keep, drop = FALSE] # [param, chain, iter]
  if (estimate_gamma) {
    beta_samples  <- as.vector(arr[1, , ])
    gamma_samples <- as.vector(arr[2, , ])
  } else {
    beta_samples  <- as.vector(arr[1, , ])
    gamma_samples <- rep(gamma, length(beta_samples))
  }
  R0_samples <- beta_samples / gamma_samples

  # --- 9) Posterior summaries ---
  if (estimate_gamma) {
    estimates <- data.frame(
      parameter = c("beta", "gamma", "R0"),
      median = c(stats::median(beta_samples), stats::median(gamma_samples), stats::median(R0_samples)),
      mean   = c(mean(beta_samples), mean(gamma_samples), mean(R0_samples)),
      sd     = c(stats::sd(beta_samples), stats::sd(gamma_samples), stats::sd(R0_samples)),
      lower_95 = c(stats::quantile(beta_samples, 0.025),
                   stats::quantile(gamma_samples, 0.025),
                   stats::quantile(R0_samples, 0.025)),
      upper_95 = c(stats::quantile(beta_samples, 0.975),
                   stats::quantile(gamma_samples, 0.975),
                   stats::quantile(R0_samples, 0.975))
    )
  } else {
    estimates <- data.frame(
      parameter = c("beta", "R0"),
      median = c(stats::median(beta_samples), stats::median(R0_samples)),
      mean   = c(mean(beta_samples), mean(R0_samples)),
      sd     = c(stats::sd(beta_samples), stats::sd(R0_samples)),
      lower_95 = c(stats::quantile(beta_samples, 0.025), stats::quantile(R0_samples, 0.025)),
      upper_95 = c(stats::quantile(beta_samples, 0.975), stats::quantile(R0_samples, 0.975))
    )
  }

  # --- 10) Deterministic fit using posterior medians ---
  beta_med  <- stats::median(beta_samples)
  gamma_med <- stats::median(gamma_samples)

  sys_fit <- dust2::dust_system_create(
    sir,
    list(N = N, I0 = I0, beta = beta_med, gamma = gamma_med),
    n_particles = 1,
    deterministic = TRUE
  )
  dust2::dust_system_set_state_initial(sys_fit)
  result_fit <- dust2::dust_system_simulate(sys_fit, times)

  S_t <- result_fit[1, ]
  cum_inc_fit <- result_fit[4, ]
  daily_inc_fit <- c(cum_inc_fit[1], diff(cum_inc_fit))

  # --- 11) R(t) credible intervals via posterior-draw simulations ---
  set.seed(1)
  draws <- data.frame(beta = beta_samples, gamma = gamma_samples)
  sel_idx <- sample(seq_len(nrow(draws)), size = min(n_ppc, nrow(draws)))
  draws_sel <- draws[sel_idx, , drop = FALSE]

  Rt_mat <- matrix(NA_real_, nrow = nrow(draws_sel), ncol = length(times))
  St_mat <- matrix(NA_real_, nrow = nrow(draws_sel), ncol = length(times))

  for (i in seq_len(nrow(draws_sel))) {
    sys_i <- dust2::dust_system_create(
      sir,
      list(N = N, I0 = I0, beta = draws_sel$beta[i], gamma = draws_sel$gamma[i]),
      n_particles = 1,
      deterministic = TRUE
    )
    dust2::dust_system_set_state_initial(sys_i)
    res_i <- dust2::dust_system_simulate(sys_i, times)
    S_i <- res_i[1, ]
    St_mat[i, ] <- S_i
    Rt_mat[i, ] <- (draws_sel$beta[i] * S_i) / (draws_sel$gamma[i] * N)
  }

  R_effective_mean <- colMeans(Rt_mat)
  R_lower <- apply(Rt_mat, 2, stats::quantile, 0.025)
  R_upper <- apply(Rt_mat, 2, stats::quantile, 0.975)
  S_t_mean <- colMeans(St_mat)

  # --- 12) Posterior predictive incidence bands (PPC) ---
  pp_inc_mat <- matrix(NA_real_, nrow = length(times), ncol = nrow(draws_sel))
  for (i in seq_len(nrow(draws_sel))) {
    sys_i <- dust2::dust_system_create(
      sir,
      list(N = N, I0 = I0, beta = draws_sel$beta[i], gamma = draws_sel$gamma[i]),
      n_particles = 1,
      deterministic = TRUE
    )
    dust2::dust_system_set_state_initial(sys_i)
    res_i <- dust2::dust_system_simulate(sys_i, times)
    cum_inc_i <- res_i[4, ]
    daily_true <- c(cum_inc_i[1], diff(cum_inc_i))
    lambda_i <- pmax(daily_true * rho, 1e-6)
    pp_inc_mat[, i] <- if (use_negbin) {
      stats::rnbinom(n = length(lambda_i), mu = lambda_i, size = nb_size)
    } else {
      stats::rpois(n = length(lambda_i), lambda = lambda_i)
    }
  }
  ppc_mean  <- rowMeans(pp_inc_mat)
  ppc_lower <- apply(pp_inc_mat, 1, stats::quantile, 0.025)
  ppc_upper <- apply(pp_inc_mat, 1, stats::quantile, 0.975)

  # --- 13) Return assembled output ---
  output <- list(
    samples = if (estimate_gamma) {
      data.frame(beta = beta_samples, gamma = gamma_samples, R0 = R0_samples)
    } else {
      data.frame(beta = beta_samples, R0 = R0_samples)
    },

    estimates = estimates,

    R_t = data.frame(
      time = times,
      S = S_t_mean,                 # posterior-averaged S(t)
      R_effective = R_effective_mean,
      R_lower = R_lower,
      R_upper = R_upper
    ),

    fitted_incidence = data.frame(
      time = times,
      observed = cases,
      fitted = daily_inc_fit,       # deterministic fit at posterior medians
      ppc_mean = ppc_mean,
      ppc_lower = ppc_lower,
      ppc_upper = ppc_upper
    ),

    mcmc_info = list(
      n_steps = n_steps,
      n_chains = n_chains,
      n_burnin = n_burnin,
      estimated_params = param_names
    ),

    fixed_parameters = list(
      N = N,
      I0 = I0,
      gamma = if (!estimate_gamma) gamma else "estimated",
      use_negbin = use_negbin,
      nb_size = nb_size,
      rho = rho
    )
  )

  class(output) <- c("epievolve_fit", "list")
  return(output)
}
