#' Build deterministic I0 posterior summary object
#' @keywords internal
I0_quantiles <- function(I0_samp, q3) {

  ## posterior summaries of I0
  I0_q <- q3(I0_samp) # 0.5, 0.975, 0.025

  ## Return objects needed downstream
  list(
    I0_q = I0_q
  )
}
