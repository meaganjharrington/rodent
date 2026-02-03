#' Helper functions for posterior calculations
#'
#' @return A list containing the helper functions:
#'   \item{map_blocks_exp}{Expands block-level log-betas to daily beta(t)}
#'   \item{q3}{Computes median and 95% credible interval}
#'
#' @keywords internal
posterior_helpers <- function() {

  ## Helpers
  # build full-length beta(t) from K log-betas
  map_blocks_exp <- function(vals_log, starts, ends, timepoints) {
    vals <- exp(vals_log); # transform to normal
    out <- numeric(timepoints)
    for (k in seq_along(starts)) out[starts[k]:ends[k]] <- vals[k]
    out
  }

  q3 <- function(x) stats::quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE) # extract median, 95% CI

  ## Return these helpers for use in other modules
  list(
    map_blocks_exp = map_blocks_exp,
    q3 = q3
  )
}
