#' Simulate SIR Model with Single Time-Varying Transmission Rate
#'
#' @param N Total population size
#' @param I0 Initial number of infected individuals
#' @param gamma Recovery rate (1/infectious period)
#' @param beta Transmission rate before change point
#' @param end_time End time of simulation period
#' @return A list containing:
#'   \item{times}{Time points}
#'   \item{S}{Susceptible counts over time}
#'   \item{I}{Infected counts over time}
#'   \item{R}{Recovered counts over time}
#'   \item{beta}{Transmission rate over time}
#'   \item{parameters}{List of parameters used}
#'
#' @export
sir_nostep <- function(N, I0, gamma, beta, end_time){

  # Create time vector
  times <- seq(0, end_time, by = 1)

  # Validate inputs
  if (N <= 0) stop("N must be positive")
  if (I0 <= 0 || I0 >= N) stop("I0 must be positive and less than N")
  if (gamma <= 0) stop("gamma must be positive")
  if (beta <= 0) stop("beta must be positive")

  # Define the SIR model
  sir <- odin2::odin({
    deriv(S) <- -beta * S * I / N
    deriv(I) <- beta * S * I / N - gamma * I
    deriv(R) <- gamma * I

    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0

    N           <- parameter()
    I0          <- parameter()
    beta        <- parameter()
    gamma       <- parameter()
  })

  pars <- list(
    N = N,
    I0 = I0,
    beta = beta,
    gamma = gamma
  )

  # create dust2 system
  sys <- dust2::dust_system_create(sir, pars)

  dust2::dust_system_set_state_initial(sys)

  # Run the deterministic simulation
  result <- dust2::dust_system_simulate(sys, times)

  # Extract state variables from matrix [state, time]
  S_vals <- result[1, ]
  I_vals <- result[2, ]
  R_vals <- result[3, ]

  # Create output object
  output <- list(
    times = times,
    S = S_vals,
    I = I_vals,
    R = R_vals,
    beta = beta,
    parameters = list(
      N = N,
      I0 = I0,
      gamma = gamma,
      beta = beta
    )
  )

  class(output) <- c("epievolve_sir", "list")
  return(output)
}
