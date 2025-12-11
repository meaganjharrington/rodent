

#' Fit One-Step SIR Model to Incidence Data
#'
#' Estimates two transmission rates (beta1, beta2) with one change point.
#' Returns mechanistic R(t) = beta(t) * S(t) / (gamma * N).
#'
#' @param incidence_data Data frame with 'time' and 'cases' columns
#' @param N Total population size
#' @param gamma Recovery rate
#' @param beta1_prior_mean Prior mean for beta1 (default: 0.4)
#' @param beta2_prior_mean Prior mean for beta2 (default: 0.2)
#' @param change_time_prior_mean Prior mean for change time (default: midpoint)
#' @param n_steps MCMC steps (default: 10000)
#' @param n_chains MCMC chains (default: 3)
#'
#' @return List with: samples, estimates, R_t (time-varying!), fitted_incidence
#' @export
fit_sir_onestep <- function(incidence_data,
                            N,
                            gamma,
                            beta1_prior_mean = 0.4,
                            beta2_prior_mean = 0.2,
                            change_time_prior_mean = NULL,
                            I0_prior_mean = NULL,
                            n_steps = 10000,
                            n_chains = 3) {

  times <- incidence_data$time
  cases <- incidence_data$cases

  if (is.null(I0_prior_mean)) I0_prior_mean <- mean(cases[1:min(5, length(cases))])
  if (is.null(change_time_prior_mean)) change_time_prior_mean <- median(times)

  # SIR model with step change
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
    beta <- if (t < change_time) beta1 else beta2
    beta1 <- parameter()
    beta2 <- parameter()
    change_time <- parameter()
    gamma <- parameter()
  })

  # Posterior
  posterior <- function(pars) {
    beta1 <- pars[1]
    beta2 <- pars[2]
    change_time <- pars[3]
    I0 <- pars[4]

    if (beta1 <= 0 || beta2 <= 0 || I0 <= 0 || I0 >= N ||
        change_time < min(times) || change_time > max(times)) {
      return(-Inf)
    }

    # Prior
    lp <- dnorm(beta1, beta1_prior_mean, 0.2, log = TRUE) +
      dnorm(beta2, beta2_prior_mean, 0.15, log = TRUE) +
      dnorm(change_time, change_time_prior_mean, 10, log = TRUE) +
      dnorm(I0, I0_prior_mean, 5, log = TRUE)

    # Likelihood
    sys <- dust2::dust_system_create(
      sir,
      list(N = N, I0 = I0, beta1 = beta1, beta2 = beta2,
           change_time = change_time, gamma = gamma),
      n_particles = 1, deterministic = TRUE
    )
    dust2::dust_system_set_state_initial(sys)
    result <- dust2::dust_system_simulate(sys, times)

    daily_inc <- c(result[4, 1], diff(result[4, ]))
    ll <- sum(dpois(cases, lambda = pmax(daily_inc, 0.01), log = TRUE))

    return(lp + ll)
  }

  # MCMC
  cat("Running MCMC for one-step model...\n")
  model <- monty::monty_model(posterior, c("beta1", "beta2", "change_time", "I0"))
  sampler <- monty::monty_sampler_random_walk(vcv = diag(c(0.01, 0.01, 2, 1)))

  mcmc <- monty::monty_sample(
    model, sampler, n_steps,
    initial = c(beta1_prior_mean, beta2_prior_mean, change_time_prior_mean, I0_prior_mean),
    n_chains = n_chains
  )

  # Extract samples
  burnin <- floor(n_steps / 2)
  beta1_samp <- as.vector(mcmc$pars[1, , (burnin + 1):n_steps])
  beta2_samp <- as.vector(mcmc$pars[2, , (burnin + 1):n_steps])
  ct_samp <- as.vector(mcmc$pars[3, , (burnin + 1):n_steps])
  I0_samp <- as.vector(mcmc$pars[4, , (burnin + 1):n_steps])

  # Posterior medians
  beta1_med <- median(beta1_samp)
  beta2_med <- median(beta2_samp)
  ct_med <- median(ct_samp)
  I0_med <- median(I0_samp)

  # Run model to get S(t)
  sys <- dust2::dust_system_create(
    sir,
    list(N = N, I0 = I0_med, beta1 = beta1_med, beta2 = beta2_med,
         change_time = ct_med, gamma = gamma),
    n_particles = 1, deterministic = TRUE
  )
  dust2::dust_system_set_state_initial(sys)
  result <- dust2::dust_system_simulate(sys, times)

  S_t <- result[1, ]
  daily_inc <- c(result[4, 1], diff(result[4, ]))

  # Time-varying beta(t) and R(t)
  beta_t <- ifelse(times < ct_med, beta1_med, beta2_med)
  R_t <- beta_t * S_t / (gamma * N)

  # Credible intervals
  beta_lower <- ifelse(times < ct_med,
                       quantile(beta1_samp, 0.025),
                       quantile(beta2_samp, 0.025))
  beta_upper <- ifelse(times < ct_med,
                       quantile(beta1_samp, 0.975),
                       quantile(beta2_samp, 0.975))

  list(
    samples = data.frame(
      beta1 = beta1_samp,
      beta2 = beta2_samp,
      change_time = ct_samp,
      I0 = I0_samp
    ),

    estimates = data.frame(
      parameter = c("beta1", "beta2", "change_time", "I0"),
      median = c(beta1_med, beta2_med, ct_med, I0_med),
      lower_95 = c(quantile(beta1_samp, 0.025), quantile(beta2_samp, 0.025),
                   quantile(ct_samp, 0.025), quantile(I0_samp, 0.025)),
      upper_95 = c(quantile(beta1_samp, 0.975), quantile(beta2_samp, 0.975),
                   quantile(ct_samp, 0.975), quantile(I0_samp, 0.975))
    ),

    R_t = data.frame(
      time = times,
      S = S_t,
      beta = beta_t,
      R_effective = R_t,
      R_lower = beta_lower * S_t / (gamma * N),
      R_upper = beta_upper * S_t / (gamma * N)
    ),

    fitted_incidence = data.frame(
      time = times,
      observed = cases,
      fitted = daily_inc
    )
  )
}


#' Print Fit Summary
#' @export
print.epievolve_fit <- function(x, ...) {
  cat("SIR Bayesian Fit\n")
  cat("================\n\n")
  cat("Parameters:\n")
  print(x$estimates, row.names = FALSE)
  cat("\nR(t) range:", round(range(x$R_t$R_effective), 2), "\n")
}
