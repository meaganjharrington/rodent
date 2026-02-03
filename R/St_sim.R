#' Build deterministic S(t) series using posterior medians
#' @keywords internal
St_sim <- function(beta_q, I0_q, gen, N, gamma, time0, time_vec) {

  # Pull medians from provided summaries
  beta_median_series <- beta_q$beta_median_series
  if (is.null(beta_median_series)) {
    stop("St_sim: 'beta_q$beta_median_series' not found")
  }

  I0_median <- I0_q$I0_median
  if (is.null(I0_median)) {
    # fallback if not provided
    I0_median <- unname(I0_q$I0_q["50%"])
  }

  # Build beta-times/values (define at time0 by repeating first value)
  beta_times_med  <- c(time0, time_vec)
  beta_values_med <- c(beta_median_series[1], beta_median_series)

  # Create and simulate deterministic system
  sys_med <- dust2::dust_system_create(
    gen,  # use gen() here if your generator requires construction
    pars = list(
      N           = N,
      I0          = I0_median,
      gamma       = gamma,
      dt          = 1,
      n_beta      = length(beta_times_med),
      n_particles = 1,
      beta_values = beta_values_med,
      beta_times  = beta_times_med
    )
  )

  dust2::dust_system_set_state_initial(sys_med)
  res_med <- dust2::dust_system_simulate(sys_med, time_vec)

  # Extract S(t)
  rn_med <- rownames(res_med)
  if (!is.null(rn_med)) {
    S_row <- which(rn_med == "S")
    if (length(S_row) != 1L) S_row <- 1L
  } else {
    S_row <- 1L
  }

  S_t <- as.numeric(res_med[S_row, ])

  list(
    S_t            = S_t,
    res_med        = res_med,
    beta_times_med = beta_times_med,
    beta_values_med= beta_values_med
  )
}
