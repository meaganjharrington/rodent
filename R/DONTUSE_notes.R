

# library(odin2)
# library(dust2)
# library(monty)
# sir <- odin({
#   update(S) <- S - n_SI
#   update(I) <- I + n_SI - n_IR
#   update(R) <- R + n_IR
#   update(incidence) <- incidence + n_SI
#
#   initial(S) <- N - I0
#   initial(I) <- I0
#   initial(R) <- 0
#   initial(incidence, zero_every = 1) <- 0
#
#   p_SI <- 1 - exp(-beta * I / N * dt)
#   p_IR <- 1 - exp(-gamma * dt)
#   n_SI <- Binomial(S, p_SI)
#   n_IR <- Binomial(I, p_IR)
#
#   N <- parameter(1000)
#   I0 <- parameter(10)
#   beta <- parameter(0.2)
#   gamma <- parameter(0.1)
#
#   cases <- data()
#   cases ~ Poisson(incidence)
# })
# sir
#
# filter <- dust_filter_create(sir, data = data, time_start = 0,
#                              n_particles = 200, dt = 0.25)
# dust_likelihood_run(filter, list(beta = 0.4, gamma = 0.2))
#
# packer <- monty_packer(c("beta", "gamma"))
# packer
#
#
#
# likelihood <- dust_likelihood_monty(filter, packer)
# likelihood
#
# prior <- monty_dsl({
#   beta ~ Exponential(mean = 0.5)
#   gamma ~ Exponential(mean = 0.3)
# })
# prior
#
#
# packer <- monty_packer(c("beta", "gamma"))
# likelihood <- dust_likelihood_monty(filter, packer)
# likelihood
#
#
# posterior <- likelihood + prior
# posterior
#
# vcv <- diag(2) * 0.2
#
# sampler <- monty_sampler_random_walk(vcv)
# sampler
#
# samples <- monty_sample(posterior, sampler, 1000, n_chains = 3)
# samples
#
# samples_df <- posterior::as_draws_df(samples)
# samples_df
#
