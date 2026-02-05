#' Package Workflow Overview
#'
#' High-level idea:
#'   1) Model transmission with an SIR system where β(t) is piecewise-constant
#'      over user-defined time blocks, γ is fixed, and the initial infectious count I0
#'      is unknown. R_t is computed as R_t = (β(t) / γ) * S(t) / N
#'   2) Infer the unknown parameters (β-blocks and I0) by sampling a posterior
#'      using a random-walk MH sampler. The likelihood is evaluated by integrating
#'      the SIR system against the observed daily incidence
#'   3) Propagate posterior uncertainty by simulating S(t) for many posterior
#'      draws and computing pointwise R_t quantiles
#'
#' Model specification
#' - States: S(t), I(t), R(t) with S + I + R = N.
#' - Transmission: β(t) = β_k for t in block k (piecewise constant).
#' - Removal/recovery: fixed γ.
#' - Initial conditions: I0 unknown (inferred), S0 = N - I0 (no pre-existing immunity
#'   in the basic version; can be extended to include R0 > 0 if needed).
#' - Reproduction number:
#'     R_t = (β(t) / γ) * S(t) / N
#'
#' Parameterisation and priors & PRIORS
#' - Parameter vector: θ = (log β_1, ..., log β_K, log I0).
#' - Priors:
#'     * log β_k ~ Normal(mean_beta, sd_beta), independently.
#'     * Random-walk smoothing prior on successive differences: (log β_k - log β_{k-1})
#'       ~ Normal(0, rw_sd_beta) for k > 1.
#'     * log I0 ~ Normal(mean_I0, sd_I0).
#'
#'1. Input Validation
#' 2. Beta‑Block Construction
#' 3. Parameter Vector Setup
#' 4. Posterior Helper Functions (block‑expansion + quantiles)
#' 5. Dust2 Likelihood Construction (generator, filter, loglikelihood)
#' 6. Prior Specification (logprior)
#' 7. Posterior Definition (density = likelihood + prior)
#' 8. MCMC Setup + Sampling (monty)
#' 9. Extract Posterior Samples (post‑burnin)
#' 10. Expand β-blocks into β(t)
#' 11. Summarise Posterior β(t) Quantiles
#' 12. Summarise I0 Posterior Quantiles
#' 13. Deterministic S(t) Simulation (dust2 system)
#' 14. Compute Rt(t) Median + 95% CI
#' 15. Construct Rt Output Objects (summary tables + plot + return list)
#'
#'
#' Workflow
#'1. validate_inputs
#'2. build_beta_timeblocks
#'3. parameter_setup
#'4. posterior_helpers
#'5. likelihood
#'6. build_prior
#'7. mcmc
#'8. extract_beta_I0
#'9. posterior_beta
#'10. I0_quantiles
#'11. St_sim
#'12. compute_Rt
#'13. output_data_plots
#'
#' To-do list:
#' Fix S_t posteriors  - not currently sampling from the posterior? Must be why Rt intervals are so wide
#'
#' Rt(t) with stepwise beta(t), fixed gamma; infer I0 and beta blocks via Metropolis-Hastings random walk MCMC
#' Using full odin2 -> dust2 -> monty workflow (deterministic likelihood)
#' log-posterior = loglikelihood + logprior
#'
#' @param incidence data.frame with columns: time (consecutive integers), cases >= 0.
#' Time can start at any integer, we keep it and set time_start to "first time - 1".
#'   - time : consecutive integers (any starting value allowed)
#'   - cases: non-negative integers (daily incidence)
#' @param N Numeric scalar > 0, population size
#' @param gamma Numeric scalar > 0, fixed recovery rate
#' R Immunity existing in population, default = 0 - add back in once baseline immunity functionality added!
#' @param beta_breaks Integer vector in ORIGINAL time units; mapped to indices internally
#' @param mcmc List with settings:
#'                - n_steps (6000)
#'                - burnin (0.5)
#'                - proposal (NULL -> diag(0.02^2))
#'                - seed (4)
#'                - n_rt_draws (300) posterior draws for Rt uncertainty bands
#' @param priors List: mean_beta=log(0.3), sd_beta=0.5, rw_sd_beta=0.2,
#'               mean_I0=log(10), sd_I0=0.5
#' @param inits  List: beta=0.3, I0=10
#' @return list(samples, estimates, Rt_series, fit, diagnostics, model_used, blocks, fixed)
#' @export
  estimate_Rt <- function(
    incidence,
    N,
    gamma,
    beta_breaks = NULL,
    mcmc  = list(n_steps = 6000, burnin = 0.5, proposal = NULL,
                 seed = 4, n_rt_draws = 300),
    priors = list(mean_beta = log(0.3), sd_beta = 0.5, rw_sd_beta = 0.2,
                  mean_I0 = log(5), sd_I0 = 0.1),
    inits  = list(beta = 0.3, I0 = 10)
  ) {

      ## 1) Validate inputs
      inp <- validate_inputs(incidence, N, gamma, mcmc, priors, inits)

      ## 2) Beta-block construction
      blocks <- build_beta_blocks(beta_breaks, inp$time_vec, inp$timepoints)

      ## 3) Parameter vector setup
      params <- parameter_setup(blocks$K, inp$init_beta, inp$init_I0, inp$proposal)

      ## 4) Posterior helper functions
      helpers <- posterior_helpers()

      ## 5) Likelihood construction (generator + loglik function)
      likelihood <- likelihood_create(
        time_vec      = inp$time_vec,
        time0         = inp$time0,
        cases         = inp$cases,
        starts        = blocks$starts,
        ends          = blocks$ends,
        timepoints    = inp$timepoints,
        N             = N,
        gamma         = gamma,
        K             = blocks$K,
        map_blocks_exp = helpers$map_blocks_exp
      )

      ## 6) Priors + posterior density
      post <- build_prior(
        K           = blocks$K,
        mean_beta   = inp$mean_beta,
        sd_beta     = inp$sd_beta,
        rw_sd_beta  = inp$rw_sd_beta,
        mean_I0     = inp$mean_I0,
        sd_I0       = inp$sd_I0,
        loglik      = likelihood$loglikelihood
      )

      ## 7) MCMC sampling
      mcmc_res <- build_mcmc(
        density_theta = post$density_theta,
        par_names     = params$par_names,
        proposal      = params$proposal,
        n_steps       = inp$n_steps,
        burnin        = inp$burnin,
        seed          = inp$seed,
        init          = params$init
      )

      ## 8) Extract post-burnin samples (beta-blocks and I0)
      extracted <- extract_beta_I0(
        extract_post = mcmc_res$extract_post,  # was mcmc_res$raw
        K           = blocks$K
      )
      ## 9) Expand beta blocks into daily beta(t)
      beta_t <- posterior_beta(
        beta_blocks_samp = extracted$beta_blocks_samp,
        starts           = blocks$starts,
        ends             = blocks$ends,
        timepoints       = inp$timepoints
      )

      ## 10) Summarise posterior beta(t) (time series)
      beta_q <- summarise_beta_posterior(
        beta_t = beta_t,
        q3     = helpers$q3
      )

      ## 11) Summarise I0 posterior
      I0_q <- I0_quantiles(
        I0_samp = extracted$I0_samp,
        q3      = helpers$q3
      )

      ## 12) Deterministic S(t) simulation (aligned to time_vec)
      S_t <- St_sim(
        beta_q   = beta_q,
        I0_q     = I0_q,
        gen      = likelihood$gen,
        N        = N,
        gamma    = gamma,
        time0    = inp$time0,
        time_vec = inp$time_vec
      )

      ## 13) Compute Rt(t) median + CI
      Rt_q <- compute_Rt(
        beta_q = beta_q,
        S_t    = S_t,
        N      = N,
        gamma  = gamma
      )

      ## 14) Rt plot + series
      Rt_outputs <- output_data_plots(
        q3        = helpers$q3,
        time_vec  = inp$time_vec,
        Rt_median = Rt_q$Rt_median,
        Rt_lower  = Rt_q$Rt_lower,
        Rt_upper  = Rt_q$Rt_upper
      )

      ## 15) Final output list
      build_output_list(
        extracted   = extracted,
        blocks      = blocks,
        beta_q      = beta_q,
        I0_q        = I0_q,
        S_t         = S_t,
        Rt_q        = Rt_q,
        Rt_plot     = Rt_outputs$Rt_plot,
        Rt_series   = Rt_outputs$Rt_series,
        inputs      = inp,
        N           = N,
        gamma       = gamma,
        diagnostics = mcmc_res$diagnostics
      )
  }
