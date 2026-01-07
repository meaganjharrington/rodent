
# SIR with step-wise beta(t) via blocks; stochastic infections & recoveries
make_sir_step_stochastic_generator <- function() {
  odin2::odin({
    ## States
    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0

    ## Per-day hazards
    p_SI  <- 1 - exp(-beta_t * I / N)
    p_rec <- 1 - exp(-gamma)  # daily recovery hazard

    ## Stochastic events (binomial draws)
    new_inf <- rand_binom(S, p_SI)
    new_rec <- rand_binom(I, p_rec)

    ## State updates
    update(S) <- S - new_inf
    update(I) <- I + new_inf - new_rec
    update(R) <- R + new_rec

    ## Observed "new cases" per day (latent)
    inc_step <- new_inf
    initial(inc_step, zero_every = 1) <- 0

    ## Time-varying beta via piecewise-constant interpolation
    beta_t <- interpolate(beta_times, beta_values, "constant")

    ## User inputs
    N <- user(); I0 <- user(); gamma <- user()
    beta_values[] <- user(); beta_times[] <- user()
    dim(beta_values) <- length(beta_times)
  })
}
