#' Build outputs to present
#' @keywords internal
build_output_list <- function(
    extracted,
    blocks,
    beta_q,
    I0_q,
    S_t,
    Rt_q,
    Rt_plot,
    Rt_series,
    inputs,
    N,
    gamma,
    diagnostics
) {
  # Metadata
  starts <- blocks$starts
  ends   <- blocks$ends

  # Samples (robust)
  beta_blocks_samp <- as.matrix(extracted$beta_blocks_samp)
  I0_samp          <- as.numeric(extracted$I0_samp)
  n <- min(nrow(beta_blocks_samp), length(I0_samp))
  beta_blocks_samp <- beta_blocks_samp[seq_len(n), , drop = FALSE]
  I0_samp          <- I0_samp[seq_len(n)]
  if (is.null(colnames(beta_blocks_samp))) {
    colnames(beta_blocks_samp) <- paste0("beta_block_", seq_len(ncol(beta_blocks_samp)))
  }
  samples <- data.frame(beta_blocks_samp, I0 = I0_samp, check.names = FALSE)

  # Rt series sanity
  stopifnot(all(c("time","median","lower","upper") %in% names(Rt_series)))
  L_Rt <- min(nrow(Rt_series))
  Rt_series <- Rt_series[seq_len(L_Rt), , drop = FALSE]

  # Beta series on the observed grid
  beta_series <- data.frame(
    time        = inputs$time_vec[seq_len(length(beta_q$beta_median_series))],
    beta_median = beta_q$beta_median_series,
    beta_lower  = beta_q$beta_lower_series,
    beta_upper  = beta_q$beta_upper_series,
    check.names = FALSE
  )

  # S(t) aligned to Rt time if off by one
  L_SRt <- min(length(S_t), nrow(Rt_series))
  St_series <- data.frame(
    time = Rt_series$time[seq_len(L_SRt)],
    St   = S_t[seq_len(L_SRt)]
  )

  list(
    samples   = samples,
    estimates = list(
      beta_series = beta_series,
      I0 = data.frame(
        I0_median = as.numeric(I0_q["50%"]),
        I0_lower  = as.numeric(I0_q["2.5%"]),
        I0_upper  = as.numeric(I0_q["97.5%"])
      )
    ),
    Rt_series = Rt_series[, c("time", "median"), drop = FALSE],
    St_series = St_series,
    Rt_lower  = Rt_series$lower,
    Rt_upper  = Rt_series$upper,
    Rt_plot   = Rt_plot,
    model_used = "SIR_deterministic_step_beta_gamma_fixed_dust2_monty",
    blocks     = list(beta_starts = starts, beta_ends = ends),
    fixed      = list(N = N, gamma = gamma),
    diagnostics = diagnostics
  )
}
