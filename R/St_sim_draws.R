#' Simulate S(t) for multiple posterior draws and summarise quantiles
#' @param beta_t_samp matrix [timepoints x n_draws] from posterior_beta()
#' @param I0_samp     numeric length n_draws (subset aligned to beta_t_samp)
#' @return list(S_draws = T x D, S_q = T x 3 with cols "2.5%","50%","97.5%")
#' @keywords internal
St_sim_draws <- function(beta_t_samp, I0_samp, gen, N, gamma, time0, time_vec, q3) {
  stopifnot(ncol(beta_t_samp) == length(I0_samp))
  T_len <- nrow(beta_t_samp)
  D     <- ncol(beta_t_samp)

  S_mat <- matrix(NA_real_, nrow = T_len, ncol = D)

  for (j in seq_len(D)) {
    beta_t_j   <- beta_t_samp[, j]
    I0_j       <- I0_samp[j]
    beta_times <- c(time0, time_vec)
    beta_vals  <- c(beta_t_j[1], beta_t_j)

    sys <- dust2::dust_system_create(
      gen,
      pars = list(
        N = N, I0 = I0_j, gamma = gamma, dt = 1,
        n_beta = length(beta_times),
        n_particles = 1,
        beta_values = beta_vals,
        beta_times  = beta_times
      ),
      deterministic = TRUE
    )
    dust2::dust_system_set_state_initial(sys)
    res <- dust2::dust_system_simulate(sys, time_vec)

    rn     <- rownames(res)
    S_row  <- if (!is.null(rn)) which(rn == "S") else 1L
    S_full <- as.numeric(res[S_row, ])
    S_t_j  <- if (length(S_full) == length(time_vec) + 1L) S_full[-1] else S_full

    S_mat[, j] <- S_t_j
  }

  # time-wise quantiles (apply over columns for each time row)
  S_q <- t(apply(S_mat, 1, q3))  # T x 3, columns named by q3

  list(S_draws = S_mat, S_q = S_q)
}
