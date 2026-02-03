#' Create dust2 filter object
#'
#' @param gen        odin2 generator object
#' @param time_vec   integer/numeric vector of observation times
#' @param cases      non-negative integer vector (incidence)
#' @param time0      scalar, strictly < first(time_vec)
#' @param dt         timestep (default 1)
#' @param n_particles number of particles (1 => deterministic path)
#' @keywords internal
likelihood_build_filter <- function(gen, time_vec, cases, time0, dt = 1, n_particles = 1) {
  dust_data <- data.frame(time = time_vec, cases = cases)
  dust2::dust_filter_create(
    gen,
    data       = dust_data,
    time_start = time0,
    dt         = dt,
    n_particles = n_particles
  )
}
