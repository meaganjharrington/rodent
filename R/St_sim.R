#' Deterministic S(t) simulation using beta median and I0 median
#' @keywords internal
St_sim <- function(beta_q, I0_q, gen, N, gamma, time0, time_vec) {
  beta_med <- as.numeric(beta_q$beta_median_series)  # ensure numeric
  stopifnot(length(beta_med) == length(time_vec))

  beta_times  <- c(time0, time_vec)
  beta_values <- c(beta_med[1], beta_med)

  I0_med <- as.numeric(I0_q["50%"])  # ensure numeric

  sys <- dust2::dust_system_create(
    gen,
    pars = list(
      N = N,
      I0 = I0_med,
      gamma = gamma,
      dt = 1,
      n_beta = length(beta_times),
      n_particles = 1,
      beta_values = beta_values,
      beta_times  = beta_times
    ),
    deterministic = TRUE
  )

  dust2::dust_system_set_state_initial(sys)
  res <- dust2::dust_system_simulate(sys, time_vec)

  rn <- rownames(res)
  S_row <- if (!is.null(rn)) which(rn == "S") else 1L
  S_full <- as.numeric(res[S_row, ])

  # If model returned time0 + time_vec, drop time0 to align
  if (length(S_full) == length(time_vec) + 1L) {
    S <- S_full[-1]
  } else {
    S <- S_full
  }

  S
}
