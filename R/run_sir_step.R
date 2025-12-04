

run_sir_step <- function(N, I0, gamma, beta1, beta2, change_time, times=seq(0,100, by=1)){

  # Validate inputs
  if (N <= 0) stop("N must be positive")
  if (I0 <= 0 || I0 >= N) stop("I0 must be positive and less than N")
  if (gamma <= 0) stop("gamma must be positive")
  if (beta1 <= 0 || beta2 <= 0) stop("beta1 and beta2 must be positive")
  if (change_time < min(times)) stop("change_time must be within time range")
  if (!all(diff(times) > 0)) stop("times must be strictly increasing")

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
    beta        <- if (time < change_time) beta1 else beta2
    beta1       <- parameter()
    beta2       <- parameter()
    change_time <- parameter()
    gamma       <- parameter()
  })

  pars <- list(
    N = N,
    I0 = I0,
    beta1 = beta1,
    beta2 = beta2,
    change_time = change_time,
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

  # Calculate beta at each time point for plotting
  beta_vals <- ifelse(times < change_time, beta1, beta2)

  # Create output object
  output <- list(
    times = times,
    S = S_vals,
    I = I_vals,
    R = R_vals,
    beta = beta_vals,
    parameters = list(
      N = N,
      I0 = I0,
      gamma = gamma,
      beta1 = beta1,
      beta2 = beta2,
      change_time = change_time
    )
  )

  class(output) <- c("epievolve_sir", "list")
  return(output)
}
