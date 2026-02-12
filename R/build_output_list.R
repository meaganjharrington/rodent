#' Build outputs to present
#' @keywords internal
build_output_list <- function(
    extracted,
    blocks,
    beta_q,
    I0_q,
    S_t,           # now accepts either a numeric vector OR a list with S_median/S_lower/S_upper
    Rt_q,
    Rt_plot,
    Rt_series,
    inputs,
    N,
    gamma,
    diagnostics,
    posterior_paths = NULL   # <--- added to accept per-draw matrices
) {
  # ---- Metadata ----
  starts <- blocks$starts
  ends   <- blocks$ends

  # ---- Samples (rows = draws) ----
  # expected shapes from earlier steps:
  #   extracted$beta_blocks_samp : K x n_draws  (blocks in rows, draws in columns)
  #   extracted$I0_samp          : length n_draws
  beta_blocks_samp <- as.matrix(extracted$beta_blocks_samp)
  I0_samp          <- as.numeric(extracted$I0_samp)

  # Orient to rows = draws
  beta_blocks_draws <- t(beta_blocks_samp)  # n_draws x K
  n_draws <- nrow(beta_blocks_draws)
  if (length(I0_samp) != n_draws) {
    # Truncate to common length just in case
    n <- min(n_draws, length(I0_samp))
    beta_blocks_draws <- beta_blocks_draws[seq_len(n), , drop = FALSE]
    I0_samp           <- I0_samp[seq_len(n)]
    n_draws <- n
  }

  # Name the beta-block columns
  if (is.null(colnames(beta_blocks_draws))) {
    colnames(beta_blocks_draws) <- paste0("beta_block_", seq_len(ncol(beta_blocks_draws)))
  }

  samples <- data.frame(beta_blocks_draws, I0 = I0_samp, check.names = FALSE)

  # ---- Rt series sanity ----
  # Expect Rt_series to contain: time, median, lower, upper (as used by output_data_plots)
  stopifnot(all(c("time", "median", "lower", "upper") %in% names(Rt_series)))
  L_Rt <- nrow(Rt_series)
  Rt_series <- Rt_series[seq_len(L_Rt), c("time", "median", "lower", "upper"), drop = FALSE]

  # ---- Beta series on the observed grid ----
  beta_series <- data.frame(
    time        = inputs$time_vec[seq_len(length(beta_q$beta_median_series))],
    beta_median = as.numeric(beta_q$beta_median_series),
    beta_lower  = as.numeric(beta_q$beta_lower_series),
    beta_upper  = as.numeric(beta_q$beta_upper_series),
    check.names = FALSE
  )

  # ---- S(t) aligned to Rt time ----
  # Backward compatibility:
  # - If S_t is a *numeric* vector, treat it as a single series (legacy median path).
  # - If S_t is a *list* with S_median/S_lower/S_upper, build a richer St_series with CIs.
  if (is.numeric(S_t)) {
    L_SRt <- min(length(S_t), nrow(Rt_series))
    St_series <- data.frame(
      time = Rt_series$time[seq_len(L_SRt)],
      St   = S_t[seq_len(L_SRt)]
    )
  } else if (is.list(S_t) && all(c("S_median", "S_lower", "S_upper") %in% names(S_t))) {
    L_SRt <- min(length(S_t$S_median), nrow(Rt_series))
    St_series <- data.frame(
      time      = Rt_series$time[seq_len(L_SRt)],
      S_median  = as.numeric(S_t$S_median[seq_len(L_SRt)]),
      S_lower   = as.numeric(S_t$S_lower [seq_len(L_SRt)]),
      S_upper   = as.numeric(S_t$S_upper [seq_len(L_SRt)]),
      check.names = FALSE
    )
  } else {
    stop("S_t must be either a numeric vector or a list with S_median, S_lower, S_upper.")
  }

  # ---- I0 quantiles data.frame ----
  I0_df <- data.frame(
    I0_median = as.numeric(I0_q["50%"]),
    I0_lower  = as.numeric(I0_q["2.5%"]),
    I0_upper  = as.numeric(I0_q["97.5%"])
  )

  # ---- Build final list ----
  out <- list(
    samples    = samples,
    estimates  = list(
      beta_series = beta_series,
      I0          = I0_df
    ),
    Rt_series  = Rt_series[, c("time", "median"), drop = FALSE],
    St_series  = St_series,
    Rt_lower   = Rt_series$lower,
    Rt_upper   = Rt_series$upper,
    Rt_plot    = Rt_plot,
    model_used = "SIR_deterministic_step_beta_gamma_fixed_dust2_monty",
    blocks     = list(beta_starts = starts, beta_ends = ends),
    fixed      = list(N = N, gamma = gamma),
    diagnostics = diagnostics,
    posterior_paths = posterior_paths   # <- included if supplied
  )

  class(out) <- c("rodent_rt", class(out))
  out
}
