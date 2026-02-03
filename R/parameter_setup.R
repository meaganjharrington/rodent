#' Set up MCMC parameter vectors (log-scale)
#'
#' @param K Number of beta blocks.
#' @param init_beta Initial beta value from user input.
#' @param init_I0 Initial I0 value from user input.
#' @param proposal Optional proposal covariance matrix (may be NULL).
#'
#' @return A list containing:
#'   \item{par_names}{Character vector of parameter names.}
#'   \item{init}{Numeric vector of initial values on log-scale.}
#'   \item{proposal}{Proposal covariance matrix (defaulted if NULL).}
#'
#' @keywords internal
parameter_setup <- function(K, init_beta, init_I0, proposal) {

  ## Define parameterisation (log)
  ## Parameter vector on log-scale: [log_beta_1....log_beta_K, log_I0]
  par_names <- c(paste0("log_beta_", seq_len(K)), "log_I0")   # parameters in MCMC state vector
  init      <- c(rep(log(init_beta), K), log(init_I0))         # initial MCMC position on log scale
  if (is.null(proposal)) proposal <- diag(rep(0.05^2, K + 1)) # default diag proposal

  ## Return for downstream use
  list(
    par_names = par_names,
    init = init,
    proposal = proposal
  )
}
