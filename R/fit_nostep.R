#' Fit Constant-Beta SIR Model to Incidence Data
#'
#' Estimates transmission rate (beta) and initial infections (I0) using
#' Bayesian inference. Returns mechanistic R(t) = beta * S(t) / (gamma * N).
#'
#' @param incidence_data Data frame with 'time' and 'cases' columns
#' @param N Total population size
#' @param gamma Recovery rate (1/infectious_period)
#' @param beta_prior_mean Prior mean for beta (default: 0.3)
#' @param beta_prior_sd Prior SD for beta (default: 0.15)
#' @param I0_prior_mean Prior mean for I0 (default: mean of first 5 cases)
#' @param I0_prior_sd Prior SD for I0 (default: 5)
#' @param n_steps MCMC steps (default: 10000)
#' @param n_chains MCMC chains (default: 3)
#'
#' @return List with: samples, estimates, R_t, fitted_incidence
#' @export
fit_sir_constant <- function(incidence_data,
                             N,
                             gamma,
                             beta_prior_mean = 0.3,
                             beta_prior_sd = 0.15,
                             I0_prior_mean = NULL,
                             I0_prior_sd = 5,
                             n_steps = 10000,
                             n_chains = 3) {

  # Validate and extract data
  if (!all(c("time", "cases") %in% names(incidence_data))) {
    stop("incidence_data must have 'time' and 'cases' columns")
  }
  times <- incidence_data$time
  cases <- incidence_data$cases

  if (is.null(I0_prior_mean)) {
    I0_prior_mean <- mean(cases[1:min(5, length(cases))])
  } #where is this from????

  # Define SIR model
  sir <- odin2::odin({
    deriv(S) <- -beta * S * I / N
    deriv(I) <- beta * S * I / N - gamma * I
    deriv(R) <- gamma * I
    deriv(Inc) <- beta * S * I / N

    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0
    initial(Inc) <- 0

    N <- parameter()
    I0 <- parameter()
    beta <- parameter()
    gamma <- parameter()
  })

  # Posterior density function
  posterior <- function(pars) {
    beta <- pars[1]
    I0 <- pars[2]

    # Bounds
    if (beta <= 0 || I0 <= 0 || I0 >= N) return(-Inf)

    # Prior
    lp <- dnorm(beta, beta_prior_mean, beta_prior_sd, log = TRUE) +
      dnorm(I0, I0_prior_mean, I0_prior_sd, log = TRUE)

    # Likelihood
    sys <- dust2::dust_system_create(sir,
                                     list(N = N, I0 = I0, beta = beta, gamma = gamma),
                                     n_particles = 1, deterministic = TRUE)
    dust2::dust_system_set_state_initial(sys)
    result <- dust2::dust_system_simulate(sys, times)

    daily_inc <- c(result[4, 1], diff(result[4, ]))
    ll <- sum(dpois(cases, lambda = pmax(daily_inc, 0.01), log = TRUE))

    return(lp + ll)
  }

  # Run MCMC
  cat("Running MCMC:", n_chains, "chains x", n_steps, "steps\n")
  model <- monty::monty_model(posterior, c("beta", "I0"))
  sampler <- monty::monty_sampler_random_walk(vcv = diag(c(0.01, 1)))

  mcmc <- monty::monty_sample(model, sampler, n_steps,
                              initial = c(beta_prior_mean, I0_prior_mean),
                              n_chains = n_chains)

  # Extract post-burnin samples
  burnin <- floor(n_steps / 2)
  beta_samples <- as.vector(mcmc$pars[1, , (burnin + 1):n_steps])
  I0_samples <- as.vector(mcmc$pars[2, , (burnin + 1):n_steps])

  # Get posterior medians
  beta_med <- median(beta_samples)
  I0_med <- median(I0_samples)

  # Run model with posterior median to get S(t)
  sys <- dust2::dust_system_create(sir,
                                   list(N = N, I0 = I0_med, beta = beta_med, gamma = gamma),
                                   n_particles = 1, deterministic = TRUE)
  dust2::dust_system_set_state_initial(sys)
  result <- dust2::dust_system_simulate(sys, times)

  S_t <- result[1, ]
  daily_inc <- c(result[4, 1], diff(result[4, ]))

  # Calculate mechanistic R(t) = beta * S(t) / (gamma * N)
  R_t <- beta_med * S_t / (gamma * N)

  # Calculate credible intervals for R(t) using posterior samples
  # For simplicity, use beta uncertainty (S(t) treated as fixed at median)
  R_samples <- beta_samples / gamma * median(S_t) / N  # Rough approximation

  # Return results
  list(
    samples = data.frame(beta = beta_samples, I0 = I0_samples),

    estimates = data.frame(
      parameter = c("beta", "I0"),
      median = c(beta_med, I0_med),
      lower_95 = c(quantile(beta_samples, 0.025), quantile(I0_samples, 0.025)),
      upper_95 = c(quantile(beta_samples, 0.975), quantile(I0_samples, 0.975))
    ),

    R_t = data.frame(
      time = times,
      S = S_t,
      R_effective = R_t,
      R_lower = quantile(beta_samples, 0.025) * S_t / (gamma * N),
      R_upper = quantile(beta_samples, 0.975) * S_t / (gamma * N)
    ),

    fitted_incidence = data.frame(
      time = times,
      observed = cases,
      fitted = daily_inc
    )
  )
}
