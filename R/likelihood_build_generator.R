#' Create odin2 generator for SIR with time-varying beta
#' @keywords internal
likelihood_build_generator <- function() {
  # Reuses your existing generator factory (which may be cached internally)
  make_sir_timevary_generator()
}
