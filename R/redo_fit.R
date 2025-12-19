#' Estimate Rt(t) from daily cases via odin2/dust2/monty (SIR)
#' Maximum flexibility for technical user to make adjustments ?
#'
#'
#' @param incidence data.frame with columns time (integer days) and cases (counts)
#' @param N population size
#' @param I0 initial infections (default: inferred from first 5 days)
#' @param gamma fixed recovery rate per day (default 1/7)
#' @param estimate_gamma logical (default FALSE)
#' @param use_negbin logical; if TRUE, use NB observation in model (default FALSE here)
#' @param nb_size NB size if use_negbin (default 10)
#' @param rho reporting fraction if use_negbin (default 1)
#' @param n_particles PF particles (default 500)
#' @param dt discrete step width (default 1)
#' @param n_steps MCMC iterations (default 4000)
#' @param burnin burn-in iterations (default n_steps/2)
#' @param chains number of chains (default 3)
#' @return list: samples, estimates, R_t (mean/lower/upper), fitted_incidence (PPC), diagnostics
redo_fit <- function(incidence,
                            N,
                            I0 = NULL,
                            gamma = 1/7,
                            estimate_gamma = FALSE,
                            use_negbin = FALSE,
                            nb_size = 10,
                            rho = 1,
                            n_particles = 500,
                            dt = 1,
                            n_steps = 4000,
                            burnin = NULL,
                            chains = 3) {


  # 1) Build odin2 generator (discrete-time SIR)
  gen <- sir_model()  # from R/odin_sir.R

  # 2) Create dust2 particle filter
  pars <- list(N = N, I0 = I0, beta = 0.3, gamma = gamma)
  if (use_negbin) {
    pars$nb_size <- nb_size
    pars$rho     <- rho
  }

  filt <- dust2::dust_filter_create(
    gen(), data = incidence, time_start = min(time),
    n_particles = n_particles, dt = dt
  )

  # 3) Likelihood for monty (bridged from dust2 filter)
  packer_names <- if (estimate_gamma) c("beta", "gamma") else "beta"
  packer <- monty::monty_packer(packer_names)
  like   <- dust2::dust_likelihood_monty(filt, packer)

  # 4) Priors (Exponential on beta, gamma)
  prior <- if (estimate_gamma) {
    monty::monty_dsl({
      beta  ~ Exponential(mean = 0.5)    # tune to pathogen context
      gamma ~ Exponential(mean = 0.2)
    })
  } else {
    monty::monty_dsl({
      beta ~ Exponential(mean = 0.5)
    })
  }

  post <- like + prior

  # 5) Sampler (RW with small VCV; you can switch to monty::sampler_adapt())
  if (is.null(burnin)) burnin <- floor(n_steps/2)
  vcv <- diag(length(packer_names)) * 0.05
  sampler <- monty::monty_sampler_random_walk(vcv)

  # Initial values close to prior means
  init <- if (estimate_gamma) c(beta = 0.3, gamma = gamma) else c(beta = 0.3)

  # 6) MCMC
  smp <- monty::monty_sample(post, sampler, n_steps, initial = init, n_chains = chains)

  # 7) Collect samples and summarise
  arr <- smp$pars[,, (burnin+1):n_steps, drop = FALSE] # [param, chain, iter]
  beta_draws  <- as.vector(arr[1,,])
  gamma_draws <- if (estimate_gamma) as.vector(arr[2,,]) else rep(gamma, length(beta_draws))
  R0_draws    <- beta_draws / gamma_draws

  # 8) Posterior predictive + Rt(t) bands via simulate per draw
  ndraw <- min(200L, length(beta_draws))
  sel   <- sample(seq_along(beta_draws), ndraw)
  Rt_mat   <- matrix(NA_real_, ndraw, length(time))
  PPC_mat  <- matrix(NA_real_, ndraw, length(time))

  for (i in seq_len(ndraw)) {
    # Update parameters and simulate a single trajectory deterministically (mean transitions)
    sys <- dust2::dust_system_create(gen(), pars = list(N = N, I0 = I0,
                                                        beta = beta_draws[sel[i]],
                                                        gamma = gamma_draws[sel[i]]),
                                     n_particles = 1, deterministic = TRUE)
    dust2::dust_system_set_state_initial(sys)
    y <- dust2::dust_system_simulate(sys, time)
    S_t <- y[1, ]                     # S state
    inc <- c(y[4,1], diff(y[4,]))     # daily incidence from cumulative
    Rt_mat[i, ] <- (beta_draws[sel[i]] * S_t) / (gamma_draws[sel[i]] * N)

    # generate posterior predictive counts
    mu <- pmax(inc * rho, 1e-6)
    PPC_mat[i, ] <- if (use_negbin) stats::rnbinom(n = length(mu), mu = mu, size = nb_size)
    else             stats::rpois(n = length(mu), lambda = mu)
  }

  Rt_mean  <- colMeans(Rt_mat)
  Rt_lower <- apply(Rt_mat, 2, stats::quantile, 0.025)
  Rt_upper <- apply(Rt_mat, 2, stats::quantile, 0.975)

  ppc_mean  <- colMeans(PPC_mat)
  ppc_lower <- apply(PPC_mat, 2, stats::quantile, 0.025)
  ppc_upper <- apply(PPC_mat, 2, stats::quantile, 0.975)

  # 9) Deterministic fit at posterior medians (for a clean line)
  beta_med  <- stats::median(beta_draws)
  gamma_med <- stats::median(gamma_draws)
  sys_med <- dust2::dust_system_create(gen(), pars = list(N = N, I0 = I0,
                                                          beta = beta_med, gamma = gamma_med),
                                       n_particles = 1, deterministic = TRUE)
  dust2::dust_system_set_state_initial(sys_med)
  y_med <- dust2::dust_system_simulate(sys_med, time)
  daily_fit <- c(y_med[4,1], diff(y_med[4,]))

  list(
    samples = data.frame(beta = beta_draws,
                         gamma = gamma_draws,
                         R0    = R0_draws),
    estimates = dplyr::bind_rows(
      tibble::tibble(parameter = "beta",
                     median = stats::median(beta_draws),
                     mean   = mean(beta_draws),
                     sd     = stats::sd(beta_draws),
                     lower_95 = stats::quantile(beta_draws, 0.025),
                     upper_95 = stats::quantile(beta_draws, 0.975)),
      tibble::tibble(parameter = "gamma",
                     median = stats::median(gamma_draws),
                     mean   = mean(gamma_draws),
                     sd     = stats::sd(gamma_draws),
                     lower_95 = stats::quantile(gamma_draws, 0.025),
                     upper_95 = stats::quantile(gamma_draws, 0.975)),
      tibble::tibble(parameter = "R0",
                     median = stats::median(R0_draws),
                     mean   = mean(R0_draws),
                     sd     = stats::sd(R0_draws),
                     lower_95 = stats::quantile(R0_draws, 0.025),
                     upper_95 = stats::quantile(R0_draws, 0.975))
    ),
    R_t = data.frame(time = time,
                     R_effective = Rt_mean,
                     R_lower = Rt_lower,
                     R_upper = Rt_upper),
    fitted_incidence = data.frame(
      time = time,
      observed = cases,
      fitted   = daily_fit,
      ppc_mean  = ppc_mean,
      ppc_lower = ppc_lower,
      ppc_upper = ppc_upper
    ),
    diagnostics = list(
      n_steps     = n_steps,
      burnin      = burnin,
      chains      = chains,
      n_particles = n_particles,
      observation = if (use_negbin) "NegBinom" else "Poisson"
    )
  )

}
