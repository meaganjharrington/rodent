#' Build block-level beta summaries and Rt plot object
#' @keywords internal
output_data_plots <- function() {

  ## Block-level beta summaries (basically just for return function)
  beta_block_q <- t(apply(beta_blocks_samp, 1, q3))

  # making a basic Rt plot!
  # Combine Rt summaries for plotting
  Rt_series <- data.frame(time = time_vec, Rt_median = Rt_median)
  Rt_lower  <- data.frame(time = time_vec, Rt_lower  = Rt_lower)
  Rt_upper  <- data.frame(time = time_vec, Rt_upper  = Rt_upper)

  rt_df <- merge(Rt_series, Rt_lower, by = "time")
  rt_df <- merge(rt_df, Rt_upper, by = "time")

  Rt_plot <- ggplot2::ggplot(rt_df, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Rt_lower, ymax = Rt_upper),
      fill  = "steelblue",
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Rt_median),
      color     = "steelblue",
      linewidth = 1
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype   = "dashed",
      color      = "grey40"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 5),
      oob    = scales::squish
    ) +
    ggplot2::labs(
      x     = "Time",
      y     = expression(R[t]),
      title = "Time-varying Reproduction Number"
    ) +
    ggplot2::theme_minimal()

  ## Return objects needed downstream
  list(
    beta_block_q = beta_block_q,
    Rt_plot      = Rt_plot,
    rt_df        = rt_df
  )
}
