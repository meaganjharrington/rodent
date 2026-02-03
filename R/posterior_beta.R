#' Build deterministic posterior summary object
#' @keywords internal
posterior_beta <- function() {

  ## Posterior summaries for beta(t) and Rt(t)

  ## Expand ALL posterior beta draws to time scale
  ## result: matrix [timepoints x n_draws]
  beta_t_samp <- apply(
    beta_blocks_samp,
    2,
    function(b) {
      out <- numeric(timepoints)
      for (k in seq_len(K))
        out[starts[k]:ends[k]] <- b[k]
      out
    }
  )

  ## posterior summaries of beta(t)
  beta_t_q <- apply(beta_t_samp, 1, q3)

  beta_median_series  <- beta_t_q["50%", ]
  beta_lower_series   <- beta_t_q["2.5%", ]
  beta_higher_series  <- beta_t_q["97.5%", ]

  ## Return objects needed downstream
  list(
    beta_t_samp         = beta_t_samp,
    beta_t_q            = beta_t_q,
    beta_median_series  = beta_median_series,
    beta_lower_series   = beta_lower_series,
    beta_higher_series  = beta_higher_series
  )
}
