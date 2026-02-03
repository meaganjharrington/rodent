
#' Build beta block indices
#' @keywords internal
build_beta_blocks <- function(beta_breaks, time_vec, timepoints){

## Beta blocks - time periods with constant transmission rate (piecewise constant)
## beta_breaks are supplied in ORIGINAL time units; map to 1-based indicies
make_breaks <- function(breaks, time_vec) { # convert calender times into array indicies
  if (!is.null(breaks)) {
    breaks <- sort(unique(as.integer(breaks)))
    indices_of_betabreaks <- match(breaks, time_vec) # indicies of beta breaks within time vector (?)
    if (any(is.na(indices_of_betabreaks)))
      stop("`beta_breaks` must align to observed times in `incidence$time`.")
    if (indices_of_betabreaks[1] != 1L)
      stop("First beta block must start at the first observation.")
    return(indices_of_betabreaks)
  }
  1L # if no breaks provided, single beta block starting at 1 (i.e. only one continuous time period, no breaks)
}

starts <- make_breaks(beta_breaks, time_vec)
ends   <- c(starts[-1] - 1L, timepoints) # block k ends just before k+1, etc
K      <- length(starts) # number of beta blocks


list(
  starts = starts,
  ends   = ends,
  K      = K
)


}
