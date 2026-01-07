
#' Rt(t) with stepwise beta(t), fixed gamma; infer I0 and beta blocks via MH-RW MCMC
#' Canonical odin2 -> dust2 -> monty workflow (deterministic likelihood)
#'
#' @param incidence data.frame with columns: time (consecutive integers), cases >= 0.
#'        Time can start at any integer; we keep it (Option A) and set time_start to "first time - 1".
#' @param N Numeric scalar > 0, population size
#' @param gamma Numeric scalar > 0, fixed recovery rate
#' @param beta_step Integer block length (e.g., 7) OR provide `beta_breaks`
#' @param beta_breaks Integer vector in ORIGINAL time units; mapped to indices internally
#' @param mcmc List: n_steps (6000), burnin (0.5), proposal (NULL -> diag(0.02^2)), seed (4),
#'              n_rt_draws (300)
#' @param priors List: mean_beta=log(0.3), sd_beta=0.5, rw_sd_beta=0.2,
#'               mean_I0=log(10), sd_I0=0.5
#' @param inits  List: beta=0.3, I0=10
#' @return list(samples, estimates, Rt_series, fit, diagnostics, model_used, blocks, fixed)
#' @export
final_estimate_Rt_step <- function(
    incidence, N, gamma,
    beta_step = NULL, beta_breaks = NULL,
    mcmc  = list(n_steps = 6000, burnin = 0.5, proposal = NULL, seed = 4, n_rt_draws = 300),
    priors = list(mean_beta = log(0.3), sd_beta = 0.5, rw_sd_beta = 0.2,
                  mean_I0 = log(10), sd_I0 = 0.5),
    inits  = list(beta = 0.3, I0 = 10)
) {
  ## ---- Validate core inputs ----
  stopifnot(is.data.frame(incidence), all(c("time","cases") %in% names(incidence)))
  time_vec <- as.integer(incidence$time)
  cases    <- as.integer(incidence$cases)
  stopifnot(all(is.finite(time_vec)), all(is.finite(cases)), all(cases >= 0))
  if (!all(diff(time_vec) == 1L)) stop("`incidence$time` must be consecutive integers.")
  stopifnot(is.numeric(N) && N > 0, is.numeric(gamma) && gamma > 0)

  Tn    <- length(cases)
  time0 <- time_vec[1L]
  time_start <- time0 - 1L  # strict requirement: first time > time_start

  ## ---- Pull settings with defaults ----
  get <- function(x, n, d) if (!is.null(x[[n]])) x[[n]] else d
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

  ## ---- Blocks (built over indices 1..Tn; independent of absolute origin) ----
  ## beta_breaks are supplied in ORIGINAL time units; map to 1-based positions
  make_breaks <- function(step, breaks, Tn) {
    if (!is.null(breaks)) {
      breaks <- sort(unique(as.integer(breaks)))
      idx <- match(breaks, time_vec)
      if (any(is.na(idx))) stop("`beta_breaks` must align to observed times in `incidence$time`.")
      if (idx[1] != 1L) stop("First beta block must start at the first observation (time_start+1).")
      return(idx)
    }
    if (is.null(step)) step <- Tn
    seq(1L, Tn, by = as.integer(step))
  }
  starts <- make_breaks(beta_step, beta_breaks, Tn)
  ends   <- c(starts[-1] - 1L, Tn)
  K      <- length(starts)

  ## ---- Parameter vector on log-scale: [log_beta_1..log_beta_K, log_I0]
  par_names <- c(paste0("log_beta_", seq_len(K)), "log_I0")
  init      <- c(rep(log(init_beta), K), log(init_I0))
  if (is.null(proposal)) proposal <- diag(rep(0.02^2, K + 1))

  ## ---- Helpers ----
  map_blocks_exp <- function(vals_log, starts, ends, Tn) {
    vals <- exp(vals_log); out <- numeric(Tn)
    for (k in seq_along(starts)) out[starts[k]:ends[k]] <- vals[k]
    out
  }
  q3 <- function(x) stats::quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE)

  ## ---- odin2 generator ----
  gen <- make_sir_timevary_generator()

  ## ---- Build a dust2 deterministic likelihood (unfilter) with your data ----
  dust_data <- data.frame(time = time_vec, cases = cases)
  unfilter  <- dust2::dust_unfilter_create(
    gen,
    data = dust_data,
    time_start = time_start,  # strictly less than first data time
    dt = 1
  )

  ## ---- Likelihood using dust2 ----
  loglik <- function(theta) {
    log_b <- theta[seq_len(K)]
    log_I <- theta[K + 1L]
    I0    <- exp(log_I); if (!is.finite(I0) || I0 <= 0) return(-Inf)
    beta_series <- map_blocks_exp(log_b, starts, ends, Tn)

    ## Augment beta_times/beta_values to include time_start (copy first beta)
    beta_times_aug  <- c(time_start, time_vec)
    beta_values_aug <- c(beta_series[1], beta_series)

    pars <- list(
      N = N, I0 = I0, gamma = gamma, dt = 1,
      n_beta = length(beta_times_aug),
      beta_values = beta_values_aug,
      beta_times  = beta_times_aug
    )

    dust2::dust_likelihood_run(unfilter, pars)
  }

  ## ---- Priors (log-scale): independent + RW smoothness + I0 ----
  logprior <- function(theta) {
    log_b <- theta[seq_len(K)]; log_I <- theta[K + 1L]
    p_b_ind <- sum(stats::dnorm(log_b, mean = mean_beta, sd = sd_beta, log = TRUE))
    p_b_rw  <- if (K > 1) sum(stats::dnorm(diff(log_b), mean = 0, sd = rw_sd_beta, log = TRUE)) else 0
    p_I     <- stats::dnorm(log_I, mean = mean_I0, sd = sd_I0, log = TRUE)
    p_b_ind + p_b_rw + p_I
  }

  ## ---- Posterior target (+ Jacobian for exp transform) ----
  density_theta <- function(theta) loglik(theta) + logprior(theta) + sum(theta)

  ## ---- MCMC via monty (random-walk MH) ----
  mod <- monty::monty_model(list(density = density_theta, parameters = par_names))
  set.seed(seed)
  sampler <- monty::monty_sampler_random_walk(proposal)
  n_burn  <- floor(burnin * n_steps)
  smp     <- monty::monty_sample(mod, sampler, n_steps, initial = init, n_chains = 1)

  ## ---- Extract post-burnin samples (shape-robust) ----

  pars <- smp$pars; dims <- dim(pars); idx <- (n_burn + 1):dims[2]
  Slog <- if (length(dims) == 3) as.matrix(pars[, idx, 1, drop = FALSE]) else as.matrix(pars[, idx, drop = FALSE])
  rownames(Slog) <- par_names

  beta_blocks_samp <- exp(Slog[paste0("log_beta_", seq_len(K)), , drop = FALSE])
  I0_samp          <- exp(Slog["log_I0", , drop = FALSE])
  beta_q           <- t(apply(beta_blocks_samp, 1, q3))
  I0_q             <- q3(I0_samp)

  ## ---- Median path & Rt(t) ----
  beta_med_blocks <- beta_q[, 1]
  beta_med_series <- map_blocks_exp(log(beta_med_blocks), starts, ends, Tn)

  ## Augment beta vectors for system creation too
  beta_times_med  <- c(time_start, time_vec)
  beta_values_med <- c(beta_med_series[1], beta_med_series)

  sys_med <- dust2::dust_system_create(
    gen(),
    pars = list(
      N = N, I0 = I0_q[1], gamma = gamma, dt = 1,
      n_beta = length(beta_times_med),
      beta_values = beta_values_med,
      beta_times  = beta_times_med
    ),
    n_particles = 1, deterministic = TRUE
  )
  dust2::dust_system_set_state_initial(sys_med)
  res_med <- dust2::dust_system_simulate(sys_med, time_vec)  # simulate over data times only

  ## Extract S(t) and incidence
  S_row <- 1L; inc_row <- 4L
  rn_med <- rownames(res_med)
  if (!is.null(rn_med)) {
    mS <- which(rn_med == "S");         if (length(mS) == 1) S_row  <- mS
    mI <- which(rn_med == "incidence"); if (length(mI) == 1) inc_row <- mI
  }
  S_t        <- as.numeric(res_med[S_row, ])
  lambda_med <- pmax(as.numeric(res_med[inc_row, ]), 1e-12)
  Rt_med     <- (beta_med_series * S_t) / (gamma * N)

  ## ---- Credible bands (thin draws) ----
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

      beta_times_j  <- c(time_start, time_vec)
      beta_values_j <- c(b_ser[1], b_ser)

      sys_j <- dust2::dust_system_create(
        gen(),
        pars = list(
          N = N, I0 = I0_j, gamma = gamma, dt = 1,
          n_beta = length(beta_times_j),
          beta_values = beta_values_j,
          beta_times  = beta_times_j
        ),
        n_particles = 1, deterministic = TRUE
      )
      dust2::dust_system_set_state_initial(sys_j)
      res_j <- dust2::dust_system_simulate(sys_j, time_vec)

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
  any_change <- function(param_by_iter) {
    diffs <- apply(param_by_iter, 1, function(row) c(NA, diff(row)))
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
    Rt_series = data.frame(time = time_vec, Rt_median = Rt_med,
                           Rt_lower = Rt_lo, Rt_upper = Rt_up),
    fit = data.frame(time = time_vec, observed = cases,
                     pred_median = lambda_med, pred_lower = y_lo, pred_upper = y_up),
    diagnostics = list(acceptance_proxy = acc),
    model_used = "SIR_deterministic_step_beta_gamma_fixed_dust2_monty",
    blocks = list(beta_starts = starts, beta_ends = ends),
    fixed = list(N = N, gamma = gamma)
  )
}
