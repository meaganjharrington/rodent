#' Build deterministic MCMC object
#' @keywords internal
build_mcmc <- function() {

  ## MCMC via monty (random-walk MH)
  mod <- monty::monty_model(list(density = density_theta, parameters = par_names)) # monty model object (w/ density and pars)(density = log posterior density to sample from)
  set.seed(seed) # fixed random number gen
  sampler <- monty::monty_sampler_random_walk(proposal) # cteaye rw MH sampler
  n_burn  <- floor(burnin * n_steps) # compute burn-in no (burnin 0.5 = first 1/2 of samples dropped)(burn-in iterations as integer count)

  sample <-
    monty::monty_sample(mod, sampler, n_steps,
                        initial = init, n_chains = 1) # runs MCMC! runs n_steps iterations

  ## Extract post-burnin samples
  pars <- sample$pars  # raw MCMC parameter array
  dims <- dim(pars)    # store array dimensions
  indices <- (n_burn + 1):dims[2] # post-burnin iterations

  extract_post <- if (length(dims) == 3) {
    matrix(
      pars[, indices, 1],
      nrow = dims[1],
      ncol = length(indices)
    )
  } else {
    as.matrix(pars[, indices, drop = FALSE])
  } # extract post burn-in samples

  ## Return everything needed downstream
  list(
    mod          = mod,
    sampler      = sampler,
    sample       = sample,
    extract_post = extract_post
  )
}
