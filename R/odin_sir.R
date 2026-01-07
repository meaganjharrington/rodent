
#' SIR (odin2) generator for constant-beta model with new infections output
#'
#' Deterministic ODE SIR system:
#'   dS/dt = -beta * S * I / N
#'   dI/dt =  beta * S * I / N - gamma * I
#'   dR/dt =  gamma * I
#'
#' Output:
#'   inc_step = beta * S * I / N   # model incidence (new infections per day)
#'
#' Parameters (user-supplied): N, I0, beta, gamma
#'
#' @return An odin2 generator function compatible with dust2::dust_system_create
#' @export
make_sir_generator <- function() {
  odin2::odin({
    ## Parameters
    N     <- parameter()
    I0    <- parameter()
    beta  <- parameter()
    gamma <- parameter()

    ## States (Deterministic ODE)
    deriv(S) <- -beta * S * I / N
    deriv(I) <-  beta * S * I / N - gamma * I
    deriv(R) <-  gamma * I

    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0

    ## Output needed for the Poisson likelihood
    output(inc_step) <- beta * S * I / N
  })
}

