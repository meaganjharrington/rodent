
#' Estimate Rt(t) from daily new cases (deterministic SIR, Poisson) — NO PRIOR
#'   Estimating BOTH beta and gamma, with user-specified fixed I0
#'
#' - Deterministic odin2/dust2 trajectory
#' - Observed NEW incidence via Poisson likelihood on inc_step
#' - Two parameters (beta, gamma) inferred with likelihood-only MCMC on log-scale
#' - I0 is FIXED and must be supplied by the user
#' - Rt(t) = beta * S(t) / (gamma * N), reported at posterior-median beta,gamma
#'
#' @param incidence data.frame with columns:
#'   - `time`: consecutive daily integers (e.g., 1,2,...,T)
#'   - `cases`: non-negative counts (NEW infections per day)
#' @param N Population size (positive scalar)
#' @param I0 Initial infected (positive scalar, REQUIRED)
#' @param n_steps MCMC iterations (default 4000)
#' @param burnin Fraction of iterations to discard (default 0.5)
#' @param proposal_vcv 2x2 proposal variance-covariance on (log_beta, log_gamma).
#'        Default diag(c(0.05, 0.05)^2).
#' @param init_beta Initial guess for beta (default 0.3)
#' @param init_gamma Initial guess for gamma (default 1/7)
#' @param seed RNG seed (default 4)
#' @param n_rt_draws Number of posterior draws to use for Rt(t) intervals
#'        (default 300; set 0 for median-only)
#'
#' @return list:
#'   - `samples`: data.frame(beta, gamma, R0)
#'   - `estimates`: medians & 95% CIs for beta, gamma, R0
#'   - `Rt_series`: data.frame(time, Rt_median, Rt_lower, Rt_upper)
#'   - `model_used`: "SIR_deterministic_beta_gamma"
#' @export
final_estimate_Rt <- function(incidence, N, I0,
                              n_steps = 4000, burnin = 0.5,
                              proposal_vcv = diag(c(0.02, 0.02)^2),
                              init_beta = 0.3, init_gamma = 1/7,
                              seed = 4, n_rt_draws = 300) {
  ## ---- Validate input ----
  stopifnot(is.data.frame(incidence),
            all(c("time", "cases") %in% names(incidence)))
  time  <- incidence$time
  cases <- as.numeric(incidence$cases)
  stopifnot(all(is.finite(time)), all(is.finite(cases)), all(cases >= 0))
  if (!all(diff(time) == 1)) {
    stop("`incidence$time` must be consecutive daily integers (e.g., 1,2,...,T).")
  }
  stopifnot(is.numeric(N) && length(N) == 1 && is.finite(N) && N > 0)
  stopifnot(is.numeric(I0) && length(I0) == 1 && is.finite(I0) && I0 > 0)

  Tn    <- length(cases)
  t_idx <- seq_len(Tn)

  ## ---- odin2 generator ----
  gen <- make_sir_generator()

  ## ---- Deterministic Poisson likelihood on NEW infections (inc_step) ----
  loglik <- function(beta, gamma) {
    if (!is.finite(beta) || beta <= 0 || !is.finite(gamma) || gamma <= 0) return(-Inf)

    sys <- dust2::dust_system_create(
      gen(),
      pars = list(N = N, I0 = I0, gamma = gamma, beta = beta),
      n_particles = 1, deterministic = TRUE
    )
    dust2::dust_system_set_state_initial(sys)
    res <- dust2::dust_system_simulate(sys, t_idx)

    # Locate inc_step (by name if available)
    inc_row <- 4L
    rn <- rownames(res)
    if (!is.null(rn)) {
      maybe <- which(rn == "inc_step")
      if (length(maybe) == 1) inc_row <- maybe
    }

    lambda <- pmax(as.numeric(res[inc_row, ]), 1e-12)
    sum(stats::dpois(x = cases, lambda = lambda, log = TRUE))
  }

  ## ---- Likelihood-only MCMC on theta = (log_beta, log_gamma) ----
  density_theta <- function(theta_vec) {
    log_beta  <- theta_vec[1]
    log_gamma <- theta_vec[2]
    beta  <- exp(log_beta)
    gamma <- exp(log_gamma)
    # Jacobian for (beta,gamma) -> (log_beta,log_gamma) is + log_beta + log_gamma
    loglik(beta, gamma) + log_beta + log_gamma
  }

  mod <- monty::monty_model(list(
    density    = density_theta,
    parameters = c("log_beta", "log_gamma")
  ))

  set.seed(seed)
  vcv <- if (is.matrix(proposal_vcv)) proposal_vcv else diag(rep(proposal_vcv, 2))
  sampler <- monty::monty_sampler_random_walk(vcv)
  init    <- c(log(init_beta), log(init_gamma))
  n_burn  <- floor(burnin * n_steps)

  smp <- monty::monty_sample(mod, sampler, n_steps, initial = init, n_chains = 1)

  ## ---- Extract samples (post-burnin) ----
  pars <- smp$pars
  if (is.null(pars)) {
    stop("No parameter samples returned by monty. Check initial values and likelihood.")
  }
  dims <- dim(pars)
  if (length(dims) < 2) {
    stop("Unexpected shape for smp$pars; expected a matrix or 3D array.")
  }

  n_steps_total <- dims[2]
  if (n_burn >= n_steps_total) {
    stop("Burn-in is >= total steps. Decrease `burnin` or increase `n_steps`.")
  }
  idx <- (n_burn + 1):n_steps_total

  # Handle both 2D ([param x steps]) and 3D ([param x steps x chains]) cases
  if (length(dims) == 3) {
    log_beta_samples  <- as.numeric(pars[1, idx, 1])
    log_gamma_samples <- as.numeric(pars[2, idx, 1])
  } else if (length(dims) == 2) {
    log_beta_samples  <- as.numeric(pars[1, idx])
    log_gamma_samples <- as.numeric(pars[2, idx])
  } else {
    stop("Unexpected shape for smp$pars; expected 2D or 3D.")
  }

  if (length(log_beta_samples) == 0 || length(log_gamma_samples) == 0) {
    stop("No samples remaining after burn-in. Reduce `burnin` or increase `n_steps`.")
  }

  # Transform to natural scale and validate
  beta_samples  <- exp(log_beta_samples)
  gamma_samples <- exp(log_gamma_samples)
  ok <- is.finite(beta_samples) & is.finite(gamma_samples) & beta_samples > 0 & gamma_samples > 0
  if (!any(ok)) {
    stop("All post-burn-in samples are invalid. Check `proposal_vcv`, initial values, and data.")
  }
  beta_samples  <- beta_samples[ok]
  gamma_samples <- gamma_samples[ok]
  R0_samples    <- beta_samples / gamma_samples

  # Summaries
  beta_med  <- stats::median(beta_samples)
  gamma_med <- stats::median(gamma_samples)
  R0_med    <- stats::median(R0_samples)

  if (!is.finite(beta_med) || !is.finite(gamma_med) || beta_med <= 0 || gamma_med <= 0) {
    stop("Median beta/gamma are not finite or non-positive; cannot compute Rt(t).")
  }

  ## ---- Rt(t): posterior-median path + (optional) credible bands ----
  # Median path
  sys_med <- dust2::dust_system_create(
    gen(),
    pars = list(N = N, I0 = I0, gamma = gamma_med, beta = beta_med),
    n_particles = 1, deterministic = TRUE
  )
  dust2::dust_system_set_state_initial(sys_med)
  res_med <- dust2::dust_system_simulate(sys_med, t_idx)

  # Get S(t)
  S_row <- 1L
  rn_med <- rownames(res_med)
  if (!is.null(rn_med)) {
    maybeS <- which(rn_med == "S")
    if (length(maybeS) == 1) S_row <- maybeS
  }
  S_t <- as.numeric(res_med[S_row, ])
  Rt_median <- (beta_med * S_t) / (gamma_med * N)

  # Credible bands via thin posterior draws
  Rt_lower <- Rt_upper <- rep(NA_real_, Tn)
  if (n_rt_draws > 0 && length(beta_samples) > 0) {
    draw_idx <- seq_along(beta_samples)
    if (length(draw_idx) > n_rt_draws) draw_idx <- sample(draw_idx, n_rt_draws)

    Rt_mat <- matrix(NA_real_, nrow = Tn, ncol = length(draw_idx))
    for (j in seq_along(draw_idx)) {
      b <- beta_samples[draw_idx[j]]
      g <- gamma_samples[draw_idx[j]]
      sys_j <- dust2::dust_system_create(
        gen(),
        pars = list(N = N, I0 = I0, gamma = g, beta = b),
        n_particles = 1, deterministic = TRUE
      )
      dust2::dust_system_set_state_initial(sys_j)
      res_j <- dust2::dust_system_simulate(sys_j, t_idx)
      S_j <- as.numeric(res_j[S_row, ])
      Rt_mat[, j] <- (b * S_j) / (g * N)
    }
    Rt_lower <- apply(Rt_mat, 1, stats::quantile, 0.025, na.rm = TRUE)
    Rt_upper <- apply(Rt_mat, 1, stats::quantile, 0.975, na.rm = TRUE)
  }

  ## ---- Return ----
  list(
    samples = data.frame(
      beta  = beta_samples,
      gamma = gamma_samples,
      R0    = R0_samples
    ),
    estimates = data.frame(
      beta_median  = unname(beta_med),
      beta_lower   = unname(stats::quantile(beta_samples, 0.025)),
      beta_upper   = unname(stats::quantile(beta_samples, 0.975)),
      gamma_median = unname(gamma_med),
      gamma_lower  = unname(stats::quantile(gamma_samples, 0.025)),
      gamma_upper  = unname(stats::quantile(gamma_samples, 0.975)),
      R0_median    = unname(R0_med),
      R0_lower     = unname(stats::quantile(R0_samples, 0.025)),
      R0_upper     = unname(stats::quantile(R0_samples, 0.975))
    ),
    Rt_series = data.frame(
      time      = incidence$time,
      Rt_median = Rt_median,
      Rt_lower  = Rt_lower,
      Rt_upper  = Rt_upper
    ),
    model_used = "SIR_deterministic_beta_gamma"
  )
}


# > # After running `out <- rt_estimate_constant(...)`
# > samps <- out$samples
# >
#   > par(mfrow = c(1, 2))
# > plot(samps$beta, type = "l", main = expression(beta), ylab = "", xlab = "draw")
# > plot(samps$gamma, type = "l", main = expression(gamma), ylab = "", xlab = "draw")
# >
#   > # Correlation / pair plot
#   > plot(samps$gamma, samps$beta, pch = 16, col = rgb(0,0,0,0.2),
#          +      xlab = expression(gamma), ylab = expression(beta),
#          +      main = "Posterior draws: beta vs gamma")
# > abline(a = 0, b = mean(samps$gamma/samps$beta), col = "red", lty = 2)
# >
#   >
#   > acc_proxy <- mean(diff(out$samples$beta) != 0)
#   > message(sprintf("~Acceptance proxy: %.1f%%", 100 * acc_proxy))
# FLAG !!! ~Acceptance proxy: 1.9%
