
make_sir_generator <- function() {
  odin2::odin({
    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0
    initial(inc_step, zero_every = 1) <- 0

    p_SI <- 1 - exp(-beta * I / N)
    p_IR <- 1 - exp(-gamma)

    n_SI <- Binomial(S, p_SI)
    n_IR <- Binomial(I, p_IR)

    update(S) <- S - n_SI
    update(I) <- I + n_SI - n_IR
    update(R) <- R + n_IR
    update(inc_step) <- inc_step + n_SI

    N     <- parameter()
    I0    <- parameter()
    beta  <- parameter()
    gamma <- parameter()

    cases <- data()
    cases ~ Poisson(inc_step)
  })

}
