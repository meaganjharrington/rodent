
#' Build deterministic post-processing object
#' @keywords internal
extract_beta_I0 <- function(extract_post, K) {
  beta_blocks_samp <- exp(extract_post[1:K, , drop = FALSE])
  I0_samp          <- as.numeric(exp(extract_post[K + 1L, ]))
  list(
    beta_blocks_samp = beta_blocks_samp,
    I0_samp          = I0_samp
  )
}
