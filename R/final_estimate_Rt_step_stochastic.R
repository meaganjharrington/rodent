
# R/final_estimate_Rt_step_pmmh.R
#' Rt(t) with stochastic SIR via PMMH; step-wise beta(t), fixed gamma & N; infer I0 + beta blocks
#'
#' @param incidence data.frame with columns: time=1..T (consecutive), cases>=0
#' @param N Numeric > 0, population size (fixed)
#' @param gamma Numeric > 0, recovery rate (fixed)
#' @param beta_step Integer block length (e.g., 7), OR use `beta_breaks`
#' @param beta_breaks Integer vector of block starts (must include 1)
#' @param mcmc List: n_steps=6000, burnin=0.5, proposal=NULL->diag(0.03^2),
#'              seed=4, n_rt_draws=300, n_particles=400
#' @param priors List: mean_beta=log(0.3), sd_beta=0.6, rw_sd_beta=0.25,
#'               mean_I0=log(10), sd_I0=0.6
#' @param inits  List: beta=0.3, I0=10
#' @return list(samples, estimates, Rt_series, fit, diagnostics, model_used, blocks, fixed)
#' @export
final_estimate_Rt_step_stochastic <- function(
    incidence, N, gamma,
    beta_step = NULL, beta_breaks = NULL,
    mcmc  = list(n_steps = 6000, burnin = 0.5, proposal = NULL,
                 seed = 4, n_rt_draws = 300, n_particles = 400),
    priors = list(mean_beta = log(0.3), sd_beta = 0.6, rw_sd_beta = 0.25,
                  mean_I0 = log(10), sd_I0 = 0.6),
    inits  = list(beta = 0.3, I0 = 10)
) {
  ## ---- Validate ----
  stopifnot(is.data.frame(incidence), all(c("time","cases") %in% names(incidence)))
  time  <- incidence$time
  cases <- as.numeric(incidence$cases)
  if (!all(diff(time) == 1)) stop("`incidence$time` must be consecutive 1..T.")
  stopifnot(is.numeric(N) && N > 0, is.numeric(gamma) && gamma > 0)
  Tn <- length(cases); t_idx <- seq_len(Tn)

  ## ---- Settings ----
  get <- function(x, n, d) if (!is.null(x[[n]])) x[[n]] else d
  n_steps    <- get(mcmc, "n_steps", 6000)
  burnin     <- get(mcmc, "burnin", 0.5)
  proposal   <- get(mcmc, "proposal", NULL)
  seed       <- get(mcmc, "seed", 4)
  n_rt_draws <- get(mcmc, "n_rt_draws", 300)
  n_particles<- get(mcmc, "n_particles", 400)

  mean_beta  <- get(priors, "mean_beta", log(0.3))
  sd_beta    <- get(priors, "sd_beta", 0.6)
  rw_sd_beta <- get(priors, "rw_sd_beta", 0.25)
  mean_I0    <- get(priors, "mean_I0", log(10))
  sd_I0      <- get(priors, "sd_I0", 0.6)

  init_beta  <- get(inits, "beta", 0.3)
  init_I0    <- get(inits, "I0", 10)

  ## ---- Blocks ----
  make_breaks <- function(step, breaks, Tn) {
    if (!is.null(breaks)) {
      breaks <- sort(unique(as.integer(breaks)))
      if (breaks[1] != 1) stop("beta_breaks must include 1.")
      if (any(breaks < 1 | breaks > Tn)) stop("beta_breaks out of range 1..T.")
      return(breaks)
    }
    if (is.null(step)) step <- Tn
    seq(1, Tn, by = as.integer(step))
  }
  starts <- make_breaks(beta_step, beta_breaks, Tn)
  ends   <- c(starts[-1] - 1L, Tn)
  K      <- length(starts)

  ## ---- Parameter vector (log-scale)
  par_names <- c(paste0("log_beta_", seq_len(K)), "log_I0")
  init      <- c(rep(log(init_beta), K), log(init_I0))
  if (is.null(proposal)) proposal <- diag(rep(0.03^2, K + 1))

  ## ---- Helpers ----
  map_blocks_exp <- function(vals_log, starts, ends, Tn) {
    vals <- exp(vals_log); out <- numeric(Tn)
    for (k in seq_along(starts)) out[starts[k]:ends[k]] <- vals[k]
    out
  }
  q3 <- function(x) stats::quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE)

  ## ---- Stochastic generator ----
  gen <- make_sir_step_stochastic_generator()

  ## ---- Monte-Carlo marginal log-likelihood (PMMH core) ----
  ## Unbiased estimator: log(mean exp(loglik_j)) with log-sum-exp stabilisation
  mc_loglik <- function(theta) {
    log_b <- theta[seq_len(K)]
    log_I <- theta[K + 1L]
    I0    <- exp(log_I); if (!is.finite(I0) || I0 <= 0) return(-Inf)
    beta_series <- map_blocks_exp(log_b, starts, ends, Tn)

    loglik_j <- numeric(n_particles)
    for (j in seq_len(n_particles)) {
      sys <- dust2::dust_system_create(
        gen(),
        pars = list(N = N, I0 = I0, gamma = gamma,
                    beta_values = beta_series, beta_times = t_idx),
        n_particles = 1, deterministic = FALSE
      )
      dust2::dust_system_set_state_initial(sys)
      res <- dust2::dust_system_simulate(sys, t_idx)

      inc_row <- 4L; rn <- rownames(res)
      if (!is.null(rn)) { w <- which(rn == "inc_step"); if (length(w) == 1) inc_row <- w }
      lambda <- pmax(as.numeric(res[inc_row, ]), 1e-12)

      loglik_j[j] <- sum(stats::dpois(x = cases, lambda = lambda, log = TRUE))
    }
    m <- max(loglik_j)
    m + log(sum(exp(loglik_j - m))) - log(n_particles)
  }

  ## ---- Priors (log-scale) ----
  logprior <- function(theta) {
    log_b <- theta[seq_len(K)]; log_I <- theta[K + 1L]
    p_b_ind <- sum(stats::dnorm(log_b, mean = mean_beta, sd = sd_beta, log = TRUE))
    p_b_rw  <- if (K > 1) sum(stats::dnorm(diff(log_b), mean = 0, sd = rw_sd_beta, log = TRUE)) else 0
    p_I     <- stats::dnorm(log_I, mean = mean_I0, sd = sd_I0, log = TRUE)
    p_b_ind + p_b_rw + p_I
  }

  ## ---- PMMH target (include Jacobian for exp transform) ----
  density_theta <- function(theta) mc_loglik(theta) + logprior(theta) + sum(theta)

  ## ---- MCMC via monty ----
  mod <- monty::monty_model(list(density = density_theta, parameters = par_names))
  set.seed(seed)
  sampler <- monty::monty_sampler_random_walk(proposal)
  n_burn  <- floor(burnin * n_steps)
  smp     <- monty::monty_sample(mod, sampler, n_steps, initial = init, n_chains = 1)

  ## ---- Extract samples (post-burnin) ----
  pars <- smp$pars; dims <- dim(pars); idx <- (n_burn + 1):dims[2]
  Slog <- if (length(dims) == 3) as.matrix(pars[, idx, 1, drop = FALSE]) else as.matrix(pars[, idx, drop = FALSE])
  rownames(Slog) <- par_names

  beta_blocks_samp <- exp(Slog[paste0("log_beta_", seq_len(K)), , drop = FALSE])
  I0_samp          <- exp(Slog["log_I0", , drop = FALSE])
  beta_q           <- t(apply(beta_blocks_samp, 1, q3))
  I0_q             <- q3(I0_samp)

  ## ---- Median path & Rt(t) (use deterministic trajectory for clarity of Rt)
  ## Optional: you can switch to averaging stochastic S(t) over particles if you prefer
  beta_med_blocks <- beta_q[, 1]
  beta_med_series <- map_blocks_exp(log(beta_med_blocks), starts, ends, Tn)

  gen_det <- make_sir_tvbeta_generator()  # deterministic generator from V2
  sys_med <- dust2::dust_system_create(
    gen_det(),
    pars = list(N = N, I0 = I0_q[1], gamma = gamma,
                beta_values = beta_med_series, beta_times = t_idx),
    n_particles = 1, deterministic = TRUE
  )
  dust2::dust_system_set_state_initial(sys_med)
  res_med <- dust2::dust_system_simulate(sys_med, t_idx)

  S_row <- 1L; inc_row <- 4L; rn_med <- rownames(res_med)
  if (!is.null(rn_med)) {
    mS <- which(rn_med == "S"); if (length(mS) == 1) S_row <- mS
    mI <- which(rn_med == "inc_step"); if (length(mI) == 1) inc_row <- mI
  }
  S_t        <- as.numeric(res_med[S_row, ])
  lambda_med <- pmax(as.numeric(res_med[inc_row, ]), 1e-12)
  Rt_med     <- (beta_med_series * S_t) / (gamma * N)

  ## ---- Credible bands via stochastic draws
  Rt_lo <- Rt_up <- rep(NA_real_, Tn); y_lo <- y_up <- rep(NA_real_, Tn)
  if (n_rt_draws > 0) {
    n_it <- ncol(Slog)
    draw <- if (n_it > n_rt_draws) sample(seq_len(n_it), n_rt_draws) else seq_len(n_it)

    Rt_mat <- matrix(NA_real_, Tn, length(draw))
    Lam    <- matrix(NA_real_, Tn, length(draw))
    for (j in seq_along(draw)) {
      b_log <- Slog[paste0("log_beta_", seq_len(K)), draw[j]]
      I0_j  <- exp(Slog["log_I0", draw[j]])
      b_ser <- map_blocks_exp(b_log, starts, ends, Tn)

      ## One stochastic trajectory per draw (you can average multiple particles if desired)
      sys_j <- dust2::dust_system_create(
        gen(),
        pars = list(N = N, I0 = I0_j, gamma = gamma,
                    beta_values = b_ser, beta_times = t_idx),
        n_particles = 1, deterministic = FALSE
      )
      dust2::dust_system_set_state_initial(sys_j)
      res_j <- dust2::dust_system_simulate(sys_j, t_idx)

      S_j <- as.numeric(res_j[S_row, ])
      lam <- pmax(as.numeric(res_j[inc_row, ]), 1e-12)
      Rt_mat[, j] <- (b_ser * S_j) / (gamma * N)
      Lam[, j]    <- lam
    }
    Rt_lo <- apply(Rt_mat, 1, stats::quantile, 0.025, na.rm = TRUE)
    Rt_up <- apply(Rt_mat, 1, stats::quantile, 0.975, na.rm = TRUE)
    y_lo  <- apply(Lam,    1, stats::quantile, 0.025, na.rm = TRUE)
    y_up  <- apply(Lam,    1, stats::quantile, 0.975, na.rm = TRUE)
  }

  ## ---- Acceptance proxy ----
  any_change <- function(mat) {
    diffs <- apply(mat, 1, function(row) c(NA, diff(row)))
    mean(colSums(!is.na(diffs) & diffs != 0) > 0)
  }
  acc <- any_change(Slog)

  ## ---- Return ----
  list(
    samples = data.frame(t(rbind(beta_blocks_samp, I0_samp))),
    estimates = list(
      beta_blocks = data.frame(
        block = seq_len(K), start = starts, end = ends,
        beta_median = beta_q[, 1], beta_lower = beta_q[, 2], beta_upper = beta_q[, 3]
      ),
      I0 = data.frame(I0_median = I0_q[1], I0_lower = I0_q[2], I0_upper = I0_q[3])
    ),
    Rt_series = data.frame(time = incidence$time, Rt_median = Rt_med,
                           Rt_lower = Rt_lo, Rt_upper = Rt_up),
    fit = data.frame(time = incidence$time, observed = cases,
                     pred_median = lambda_med, pred_lower = y_lo, pred_upper = y_up),
    diagnostics = list(acceptance_proxy = acc, n_particles = n_particles),
    model_used = "SIR_stochastic_step_beta_gamma_fixed_PMMH",
    blocks = list(beta_starts = starts, beta_ends = ends),
    fixed = list(N = N, gamma = gamma)
  )
}
