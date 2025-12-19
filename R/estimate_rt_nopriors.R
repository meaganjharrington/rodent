
#' Estimate Rt(t) from daily cases via odin2/dust2/monty (constant beta, Poisson) — NO PRIOR
#'
#' @param incidence data.frame with columns 'time' (integer days) and 'cases' (counts)
#' @param N population size (scalar)
#' @param gamma fixed recovery rate per day (default 1/7)
#' @param n_steps MCMC iterations (default 2000)
#' @param burnin fraction of iterations to discard (default 0.5)
#' @param proposal_sd random-walk proposal SD for beta (default 0.05)
#' @param I0 optional initial infections; if NULL, inferred from first ~5 days
#' @return list with: samples (beta, R0), estimates (medians & CIs),
#'         Rt_series (time, Rt at posterior-median beta)
#' @export
estimate_rt_constant_ml <- function(incidence, N, gamma = 1/7,
                                    n_steps = 2000, burnin = 0.5,
                                    proposal_sd = 0.05,
                                    I0 = NULL) {
  stopifnot(all(c("time", "cases") %in% names(incidence)))

  incidence <- as.data.frame(incidence)

  incidence$time <- seq_len(nrow(incidence))
  cases <- as.numeric(incidence$cases)
  Tn    <- length(cases)

  # Default I0 from early cases if not provided
  if (is.null(I0)) {
    I0 <- max(1, round(mean(head(cases, min(5, Tn))) * (1 / gamma)))
  }

  # odin2 generator (must define inc_step; make_sir_generator() in your package)
  gen <- make_sir_generator()

  # ---- Pure Poisson likelihood (NO PRIOR) ----
  loglik <- function(beta) {
    if (!is.finite(beta) || beta <= 0) return(-Inf)

    # Deterministic trajectory with proposed beta
    sys <- dust2::dust_system_create(
      gen(),
      pars = list(N = N, I0 = I0, gamma = gamma, beta = beta),
      n_particles = 1, deterministic = TRUE
    )
    dust2::dust_system_set_state_initial(sys)
    res <- dust2::dust_system_simulate(sys, incidence$time)

    inc_t  <- as.numeric(res[4, ])               # inc_step index (S=1, I=2, R=3, inc_step=4)
    lambda <- pmax(inc_t, 1e-12)                 # guard against log(0)

    sum(stats::dpois(x = cases, lambda = lambda, log = TRUE))
  }

  # Wrap into a monty model with likelihood-only density
  mod <- monty::monty_model(list(
    density    = function(theta_vec) loglik(theta_vec[1]),
    parameters = "beta"
  ))

  set.seed(4)  # reproducible
  sampler <- monty::monty_sampler_random_walk(matrix(proposal_sd^2, 1, 1))
  init    <- c(beta = 0.3)                       # neutral starting value; adjust if needed
  n_burn  <- floor(burnin * n_steps)

  smp <- monty::monty_sample(mod, sampler, n_steps, initial = init, n_chains = 1)

  # Extract post-burnin beta samples (dims: [param, samples, chains])
  pars         <- smp$pars
  beta_samples <- as.vector(pars[1, (n_burn + 1):dim(pars)[2], 1])
  R0_samples   <- beta_samples / gamma

  # Posterior (likelihood-driven)
  beta_med <- stats::median(beta_samples)
  R0_med   <- stats::median(R0_samples)

  # ---- Rt(t) curve at posterior-median beta ----
  sys_med <- dust2::dust_system_create(
    gen(),
    pars = list(N = N, I0 = I0, gamma = gamma, beta = beta_med),
    n_particles = 1, deterministic = TRUE
  )
  dust2::dust_system_set_state_initial(sys_med)
  res_med <- dust2::dust_system_simulate(sys_med, incidence$time)

  S_t <- as.numeric(res_med[1, ])                 # susceptible trajectory
  Rt  <- (beta_med * S_t) / (gamma * N)           # Rt(t) = beta * S / (gamma * N)

  # Return
  list(
    samples = data.frame(beta = beta_samples, R0 = R0_samples),
    estimates = data.frame(
      beta_median = unname(beta_med),
      beta_lower  = unname(stats::quantile(beta_samples, 0.025)),
      beta_upper  = unname(stats::quantile(beta_samples, 0.975)),
      R0_median   = unname(R0_med),
      R0_lower    = unname(stats::quantile(R0_samples, 0.025)),
      R0_upper    = unname(stats::quantile(R0_samples, 0.975))
    ),
    Rt_series = data.frame(
      time = incidence$time,
      Rt   = Rt
    )
  )

}
