
#' Simple V1 wrapper: deterministic SIR Rt(t), estimating beta & gamma (no priors)
#' @export
rt_estimate_constant <- function(incidence, N, I0) {
  final_estimate_Rt(
    incidence = incidence,
    N = N,
    I0 = I0,
    n_steps = 4000,
    burnin = 0.5,
    proposal_vcv = diag(c(0.05, 0.05)^2),
    init_beta = 0.3,
    init_gamma = 1/7,
    seed = 4,
    n_rt_draws = 300
  )
}
