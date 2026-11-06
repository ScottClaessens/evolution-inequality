#' Plot time window calculations from ancestral states
#'
#' @param time_windows Tibble of time window calculations
#'
#' @returns A patchwork of ggplot objects
#'
plot_time_windows <- function(time_windows) {

  # plot node counts through time
  left_plot <-
    time_windows |>
    rowwise() |>
    mutate(
      n = n / 100,
      n_egalitarian = median(n_egalitarian / 100),
      n_wealth_divided = median(n_wealth_divided / 100),
      n_stratified = median(n_stratified / 100)
    ) |>
    ggplot() +
    geom_ribbon(
      aes(
        x = time,
        ymin = 0,
        ymax = n_egalitarian
      ),
      fill = "#fbe4cc"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = n_egalitarian,
        ymax = n_egalitarian + n_wealth_divided,
      ),
      fill = "#fcad6b"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = n_egalitarian + n_wealth_divided,
        ymax = n,
      ),
      fill = "#e6550e"
    ) +
    geom_line(
      aes(
        x = time,
        y = n
      )
    ) +
    scale_x_continuous(
      name = "Time before present (kya)",
      limits = c(min(time_windows$time), max(time_windows$time))
    ) +
    scale_y_continuous(
      name = "Number of ancestral nodes\nwithin last 1000 years"
    ) +
    theme_classic()

  # plot ratio of stratified to egalitarian nodes
  right_plot <-
    time_windows |>
    rowwise() |>
    mutate(
      max_dist_from_equality_p3 = median(max_dist_from_equality_p3),
      estimate = median(ratio),
      lower50 = quantile(ratio, 0.25),
      upper50 = quantile(ratio, 0.75),
      lower95 = quantile(ratio, 0.025),
      upper95 = quantile(ratio, 0.975)
    ) |>
    ggplot() +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower95,
        ymax = upper95
      ),
      fill = "grey90"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower50,
        ymax = upper50
      ),
      fill = "grey80"
    ) +
    geom_line(
      aes(
        x = time,
        y = estimate
      )
    ) +
    geom_ribbon(
      aes(
        x = ifelse(time == -20, min(time_windows$time) - 0.1, time),
        ymin = lower95 - 0.01,
        ymax = upper95 + 0.01,
        alpha = max_dist_from_equality_p3
      ),
      fill = "white"
    ) +
    scale_x_continuous(
      name = "Time before present (kya)",
      limits = c(min(time_windows$time) - 0.1, max(time_windows$time))
    ) +
    scale_y_continuous(
      name = "Ratio of stratified nodes\nto egalitarian nodes",
      limits = c(0, 1.25),
      breaks = seq(0, 1.25, by = 0.25)
    ) +
    scale_alpha_continuous(
      name = NULL,
      range = c(0.95, 0)
    ) +
    theme_classic() +
    theme(legend.position = "none")

  # put together
  out <- left_plot + right_plot + plot_annotation(tag_levels = "a")

  # save
  ggsave(
    filename = "plots/time.pdf",
    plot = out,
    height = 3.5,
    width = 6,
    device = cairo_pdf
  )

  # return
  out

}
