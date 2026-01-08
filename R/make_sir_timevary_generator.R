
#' Discrete-time deterministic SIR with stepwise beta(t) and Poisson observation
make_sir_timevary_generator <- function() {
  odin2::odin({
    ## External scalars
    N     <- parameter()
    I0    <- parameter() #estimated by MCMC (?)
    gamma <- parameter()
    n_beta <- parameter() # number of beta transition time points
    # R0 <- parameter(0) # default 0

    ## Initial states
    initial(S) <- N - I0 # - R0 # Initial S, minus I0 (infectious) and R (existing immunity)
    initial(I) <- I0
   # initial(R) <- R0 # allow user to specify existing immunity, default 0
    initial(incidence, zero_every = 1) <- 0  # New infections this time step, reset every step

    ## Hazards per step (dt)
    p_SI <- 1 - exp(-beta_t * I / N * dt)  # Probability a susceptible gets infected in this time step
    p_IR <- 1 - exp(-gamma * dt)           # Probability of recovery in this time step

    ## Deterministic expected flows
    new_inf <- S * p_SI # susceptibles * infection probability
    new_rec <- I * p_IR # recovered = infections * recovery probability

    ## Updates
    update(incidence) <- new_inf #new infections in time step
    update(S)        <- S - new_inf # new infections = dec S
    update(I)        <- I + new_inf - new_rec
    #update(R)        <- R + new_rec

    ## Time-varying beta(t) via piecewise-constant interpolation
    beta_times  <- parameter()
    beta_values <- parameter()
    dim(beta_times)  <- n_beta
    dim(beta_values) <- n_beta
    beta_t <- interpolate(beta_times, beta_values, "constant") # transmission rate at time point

    ## Observation model on incidence
    cases <- data()
    cases ~ Poisson(incidence)
  })
}
