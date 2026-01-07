
#' Discrete-time deterministic SIR with stepwise beta(t) and Poisson observation
make_sir_timevary_generator <- function() {
  odin2::odin({
    ## States
    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0
    initial(incidence) <- 0  # new infections per step

    ## Hazards per step (dt)
    p_SI <- 1 - exp(-beta_t * I / N * dt)  # infection
    p_IR <- 1 - exp(-gamma * dt)           # recovery

    ## Deterministic expected flows
    new_inf <- S * p_SI
    new_rec <- I * p_IR

    ## Updates
    update(incidence) <- new_inf
    update(S)        <- S - new_inf
    update(I)        <- I + new_inf - new_rec
    update(R)        <- R + new_rec

    ## Time-varying beta(t) via piecewise-constant interpolation
    n_beta <- parameter()
    beta_times  <- parameter()
    beta_values <- parameter()
    dim(beta_times)  <- n_beta
    dim(beta_values) <- n_beta
    beta_t <- interpolate(beta_times, beta_values, "constant")

    ## External scalars
    N     <- parameter()
    I0    <- parameter()
    gamma <- parameter()

    ## Observation model on incidence
    cases <- data()
    cases ~ Poisson(incidence)
  })
}
