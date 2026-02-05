#' Prepare Rt output data and plot
#' @keywords internal
output_data_plots <- function(q3, time_vec, Rt_median, Rt_lower, Rt_upper) {
  # Align defensively
  L <- min(length(time_vec), length(Rt_median), length(Rt_lower), length(Rt_upper))
  time_vec  <- time_vec[seq_len(L)]
  Rt_median <- Rt_median[seq_len(L)]
  Rt_lower  <- Rt_lower[seq_len(L)]
  Rt_upper  <- Rt_upper[seq_len(L)]

  Rt_series <- data.frame(
    time   = time_vec,
    median = Rt_median,
    lower  = Rt_lower,
    upper  = Rt_upper
  )

  library(ggplot2)
  Rt_plot <- ggplot(Rt_series, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = median), color = "steelblue", linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_y_continuous(limits = c(0, 5), oob = scales::squish) +
    labs(x = "Time", y = expression(R[t]), title = "Time-varying Reproduction Number") +
    theme_minimal()

  list(
    Rt_series = Rt_series,
    Rt_plot   = Rt_plot
  )
}
