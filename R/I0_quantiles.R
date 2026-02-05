#' Build deterministic I0 posterior summary object
#' @keywords internal
I0_quantiles <- function(I0_samp, q3) {
  q3(I0_samp)  # return the named numeric vector: c("2.5%","50%","97.5%")
}
