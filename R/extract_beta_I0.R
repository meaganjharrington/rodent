#' Build deterministic post-processing object
#' @keywords internal
extract_beta_I0 <- function() {

  # Extract beta blocks and I0 using indices
  beta_blocks_samp <- exp(extract_post[1:K, , drop = FALSE])
  I0_samp          <- exp(extract_post[K + 1, , drop = FALSE])

  ## Return objects needed downstream
  list(
    beta_blocks_samp = beta_blocks_samp,
    I0_samp = I0_samp
  )
}
