#' Simple one-chain MH for constant-beta SIR: estimate beta and I0
#'
#' Component-wise Metropolis–Hastings random walk on log(beta) and log(I0).
#' Gamma is fixed. Uses a deterministic discrete-time SIR (dt = 1).
#' Observation model: Negative Binomial (default) or Poisson.
#'
#' @param cases Numeric vector of non-negative counts (length T).
#' @param N Total population size (scalar).
#' @param gamma Recovery rate per time step (e.g., ~0.14 per day).
#' @param I0 Optional initial infections (prevalence). If NULL (default),
#'   we set a prior mean I0_prior_mean = mean(early cases) * (1/gamma).
#' @param use_negbin Logical; if TRUE (default) use NegBin(size=nb_size), else Poisson.
#' @param nb_size NegBin size (dispersion); larger => less overdispersion (default 10).
#' @param n_steps Total MH iterations (default 4000).
#' @param burnin_frac Fraction to discard as burn-in (default 0.5).
#' @param proposal_sd_logbeta RW proposal SD on log(beta) (default 0.03–0.05 often good).
#' @param proposal_sd_logI0 RW proposal SD on log(I0) (default 0.10).
#' @param prior_meanlog_beta Prior mean of log(beta) (default log(0.30)).
#' @param prior_sdlog_beta   Prior SD   of log(beta) (default 0.50).
#' @param I0_prior_mean Optional prior mean for I0 (scalar). If NULL, computed from data.
#' @param I0_prior_sdlog Prior SD of log(I0) (default 0.70; fairly diffuse).
#' @param n_ppc Posterior draws for PPC/R(t) bands (default 200).
#' @param seed Optional RNG seed for reproducibility (default NULL).
#' @return List with samples (beta, I0, R0), estimates, R_t bands, fitted_incidence, diagnostics.
#' @export
retry_fit_nostep <- function(cases,
                             N,
                             gamma = 0.14,
                             I0 = NULL,
                             use_negbin = TRUE,
                             nb_size = 10,
                             n_steps = 4000,
                             burnin_frac = 0.5,
                             proposal_sd_logbeta = 0.03,
                             proposal_sd_logI0 = 0.10,
                             prior_meanlog_beta = log(0.30),
                             prior_sdlog_beta = 0.50,
                             I0_prior_mean = NULL,
                             I0_prior_sdlog = 0.70,
                             n_ppc = 200,
                             seed = NULL) {

  # ---- Validate inputs ----
  stopifnot(is.numeric(cases), all(is.finite(cases)), all(cases >= 0),
            length(N) == 1, N > 0, length(gamma) == 1, gamma > 0,
            length(nb_size) == 1, nb_size > 0,
            burnin_frac > 0, burnin_frac < 1,
            n_steps >= 1000)

  if (!is.null(seed)) set.seed(seed)

  Tn <- length(cases)
  times <- seq_len(Tn) - 1L

  # ---- Prior mean for I0 if not provided ----
  if (is.null(I0_prior_mean)) {
    early_mean <- mean(cases[seq_len(min(5, Tn))])
    I0_prior_mean <- max(1, round(early_mean * (1 / gamma)))
  }
  # If user gave I0 (fixed), initialise near it; otherwise we estimate it.
  # We will *always* estimate I0 here (as requested); I0 is only used to set initialisation.
  if (is.null(I0)) I0 <- I0_prior_mean

  # Clamp to sensible bounds
  I0 <- min(max(I0, 1), N - 1)

  S0_from_I0 <- function(I0_val) {
    s <- N - I0_val
    if (!is.finite(s) || s <= 0) s <- 1e-6
    s
  }

  # ---- Deterministic SIR simulator (dt = 1) ----
  sim_sir_daily <- function(beta, I0_val) {
    S <- numeric(Tn); I <- numeric(Tn); R <- numeric(Tn); inc <- numeric(Tn)
    S[1] <- S0_from_I0(I0_val); I[1] <- I0_val; R[1] <- 0
    inc[1] <- beta * S[1] * I[1] / N
    for (t in 2:Tn) {
      new_inf <- beta * S[t-1] * I[t-1] / N
      recov   <- gamma * I[t-1]
      S[t] <- S[t-1] - new_inf
      I[t] <- I[t-1] + new_inf - recov
      R[t] <- R[t-1] + recov
      inc[t] <- max(new_inf, 0)
    }
    list(S = S, I = I, R = R, inc = pmax(inc, 1e-8))
  }

  # ---- Log-likelihood ----
  loglik <- function(beta, I0_val) {
    # basic bounds
    if (!is.finite(beta) || beta <= 0) return(-Inf)
    if (!is.finite(I0_val) || I0_val <= 0 || I0_val >= N) return(-Inf)
    sim <- sim_sir_daily(beta, I0_val)
    lambda <- sim$inc
    if (use_negbin) {
      sum(stats::dnbinom(cases, size = nb_size, mu = lambda, log = TRUE))
    } else {
      sum(stats::dpois(cases, lambda = lambda, log = TRUE))
    }
  }

  # ---- Log-priors ----
  logprior_beta  <- function(logbeta) stats::dnorm(logbeta, mean = prior_meanlog_beta, sd = prior_sdlog_beta, log = TRUE)
  logprior_logI0 <- function(logI0)   stats::dnorm(logI0,   mean = log(I0_prior_mean),  sd = I0_prior_sdlog,    log = TRUE)

  # ---- Initial values (log-scale) ----
  logbeta <- numeric(n_steps); logI0 <- numeric(n_steps)
  logbeta[1] <- prior_meanlog_beta
  logI0[1]   <- log(I0)

  # Compute initial log-posterior
  curr_beta <- exp(logbeta[1]); curr_I0 <- exp(logI0[1])
  curr_logpost <- loglik(curr_beta, curr_I0) + logprior_beta(logbeta[1]) + logprior_logI0(logI0[1])

  # ---- MH loop (component-wise: beta then I0 each iteration) ----
  acc_beta <- logical(n_steps - 1L)
  acc_I0   <- logical(n_steps - 1L)

  for (i in 2:n_steps) {
    # 1) Propose log(beta)
    prop_logbeta <- logbeta[i - 1] + stats::rnorm(1, 0, proposal_sd_logbeta)
    prop_beta    <- exp(prop_logbeta)
    prop_logpost_beta <- loglik(prop_beta, curr_I0) + logprior_beta(prop_logbeta) + logprior_logI0(logI0[i - 1])
    a_beta <- prop_logpost_beta - curr_logpost
    if (log(stats::runif(1)) < a_beta) {
      logbeta[i]   <- prop_logbeta
      curr_beta    <- prop_beta
      curr_logpost <- prop_logpost_beta
      acc_beta[i - 1] <- TRUE
    } else {
      logbeta[i] <- logbeta[i - 1]
      acc_beta[i - 1] <- FALSE
    }

    # 2) Propose log(I0)
    prop_logI0 <- logI0[i - 1] + stats::rnorm(1, 0, proposal_sd_logI0)
    prop_I0    <- exp(prop_logI0)
    prop_logpost_I0 <- loglik(curr_beta, prop_I0) + logprior_beta(logbeta[i]) + logprior_logI0(prop_logI0)
    a_I0 <- prop_logpost_I0 - curr_logpost
    if (log(stats::runif(1)) < a_I0) {
      logI0[i]     <- prop_logI0
      curr_I0      <- prop_I0
      curr_logpost <- prop_logpost_I0
      acc_I0[i - 1] <- TRUE
    } else {
      logI0[i] <- logI0[i - 1]
      acc_I0[i - 1] <- FALSE
    }
  }

  # ---- Post-burn-in ----
  n_burnin <- floor(burnin_frac * n_steps)
  keep <- (n_burnin + 1):n_steps
  logbeta_post <- logbeta[keep]
  logI0_post   <- logI0[keep]

  beta_post <- exp(logbeta_post)
  I0_post   <- exp(logI0_post)
  R0_post   <- beta_post / gamma

  # ---- Summaries ----
  q <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
  summ <- rbind(
    data.frame(parameter = "beta", lower = q(beta_post)[1], median = q(beta_post)[2], upper = q(beta_post)[3]),
    data.frame(parameter = "I0",   lower = q(I0_post)[1],   median = q(I0_post)[2],   upper = q(I0_post)[3]),
    data.frame(parameter = "R0",   lower = q(R0_post)[1],   median = q(R0_post)[2],   upper = q(R0_post)[3])
  )

  # ---- Fitted trajectory at posterior medians ----
  beta_med <- stats::median(beta_post)
  I0_med   <- stats::median(I0_post)
  sim_med  <- sim_sir_daily(beta_med, I0_med)
  Rt_med   <- beta_med * sim_med$S / (gamma * N)

  # ---- R(t) and PPC bands from posterior draws ----
  nd <- min(n_ppc, length(beta_post))
  idx <- sample.int(length(beta_post), nd)

  Rt_mat <- matrix(NA_real_, nd, Tn)
  inc_mat_det <- matrix(NA_real_, nd, Tn) # deterministic means
  for (j in seq_len(nd)) {
    b <- beta_post[idx[j]]
    i0 <- I0_post[idx[j]]
    simj <- sim_sir_daily(b, i0)
    Rt_mat[j, ]      <- (b * simj$S) / (gamma * N)
    inc_mat_det[j, ] <- simj$inc
  }

  Rt_mean  <- colMeans(Rt_mat)
  Rt_lower <- apply(Rt_mat, 2, stats::quantile, 0.025)
  Rt_upper <- apply(Rt_mat, 2, stats::quantile, 0.975)

  # Posterior predictive counts (NB or Poisson)
  pp_mat <- matrix(NA_real_, nd, Tn)
  for (j in seq_len(nd)) {
    mu <- pmax(inc_mat_det[j, ], 1e-8)
    pp_mat[j, ] <- if (use_negbin) stats::rnbinom(Tn, mu = mu, size = nb_size) else stats::rpois(Tn, lambda = mu)
  }
  ppc_mean  <- colMeans(pp_mat)
  ppc_lower <- apply(pp_mat, 1, stats::quantile, 0.025)
  ppc_upper <- apply(pp_mat, 1, stats::quantile, 0.975)

  # ---- Output ----
  list(
    samples = data.frame(beta = beta_post, I0 = I0_post, R0 = R0_post),
    estimates = summ,
    R_t = data.frame(
      time = times,
      R_effective = Rt_mean, R_lower = Rt_lower, R_upper = Rt_upper,
      R_effective_median_betaI0 = Rt_med
    ),
    fitted_incidence = data.frame(
      time = times,
      observed = cases,
      fitted = as.numeric(sim_med$inc),
      ppc_mean = ppc_mean, ppc_lower = ppc_lower, ppc_upper = ppc_upper
    ),
    diagnostics = list(
      acceptance_beta = mean(acc_beta),
      acceptance_I0   = mean(acc_I0),
      n_steps = n_steps, burnin = n_burnin,
      proposal_sd_logbeta = proposal_sd_logbeta,
      proposal_sd_logI0   = proposal_sd_logI0,
      use_negbin = use_negbin, nb_size = nb_size,
      N = N, gamma = gamma,
      I0_prior_mean = I0_prior_mean, I0_prior_sdlog = I0_prior_sdlog
    )
  )
}
