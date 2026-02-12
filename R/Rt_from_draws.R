#' Compute Rt from draw matrices, then summarise quantiles
#' @param beta_t_samp T x D matrix
#' @param S_draws     T x D matrix
#' @return list(R_draws = T x D, R_q = T x 3)
#' @keywords internal
Rt_from_draws <- function(beta_t_samp, S_draws, N, gamma, q3) {
  stopifnot(all(dim(beta_t_samp) == dim(S_draws)))
  R_mat <- (beta_t_samp / gamma) * (S_draws / N)  # element-wise
  R_q   <- t(apply(R_mat, 1, q3))                 # T x 3
  list(R_draws = R_mat, R_q = R_q)
}
