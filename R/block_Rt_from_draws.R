# Summarise Rt over non-overlapping blocks/windows using posterior draws
# R_draws: matrix [T x D] from your pipeline (out$posterior_paths$R)
# starts, ends: integer vectors (1-based, inclusive) of equal length K
block_Rt_from_draws <- function(R_draws, starts, ends, probs = c(0.025, 0.5, 0.975)) {
  stopifnot(is.matrix(R_draws), length(starts) == length(ends))
  T_len <- nrow(R_draws)
  stopifnot(all(starts >= 1L), all(ends <= T_len), all(ends >= starts))

  K <- length(starts)
  out <- lapply(seq_len(K), function(k) {
    idx <- starts[k]:ends[k]
    # mean within the block per draw (D values)
    per_draw_block_mean <- colMeans(R_draws[idx, , drop = FALSE])
    qs <- stats::quantile(per_draw_block_mean, probs = probs, names = FALSE, type = 8)
    data.frame(
      block   = k,
      t_start = starts[k],
      t_end   = ends[k],
      Mean    = mean(per_draw_block_mean),
      Lower   = qs[1],
      Median  = qs[2],
      Upper   = qs[3]
    )
  })
  do.call(rbind, out)
}
