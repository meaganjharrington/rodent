
#' Build deterministic S(t) series using posterior medians
#' @keywords internal
St_sim <- function() {

## deterministic S(t) using posterior median beta(t) and I0
beta_times_med  <- c(time0, time_vec)
beta_values_med <- c(beta_median_series[1], beta_median_series)

sys_med <- dust2::dust_system_create(
  gen,
  pars = list(
    N = N,
    I0 = I0_q[1],
    gamma = gamma,
    dt = 1,
    n_beta = length(beta_times_med),
    n_particles = 1,
    beta_values = beta_values_med,
    beta_times  = beta_times_med
  ),
  deterministic = TRUE
)

dust2::dust_system_set_state_initial(sys_med)
res_med <- dust2::dust_system_simulate(sys_med, time_vec)

## pull S(t)
rn_med <- rownames(res_med)
if (!is.null(rn_med)) {
  S_row <- which(rn_med == "S")
} else {
  S_row <- 1L
}

S_t <- as.numeric(res_med[S_row, ])


## Return objects needed downstream
list(
  S_t = S_t,
  res_med = res_med,
  beta_times_med = beta_times_med,
  beta_values_med = beta_values_med
)
}

