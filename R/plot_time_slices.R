#' Plot time slice calculations from ancestral states
#'
#' @param time_slices Tibble of time slice calculations
#'
#' @returns A patchwork of ggplot objects
#'
plot_time_slices <- function(time_slices) {

  # plot lineage counts through time
  left_plot <-
    time_slices |>
    unnest(c(n_egalitarian, n_wealth_divided, n_stratified)) |>
    group_by(time) |>
    summarise(
      n_egalitarian = mean(n_egalitarian),
      n_wealth_divided = mean(n_wealth_divided),
      n_stratified = mean(n_stratified)
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
        ymax = n_egalitarian + n_wealth_divided + n_stratified,
      ),
      fill = "#e6550e"
    ) +
    geom_line(
      aes(
        x = time,
        y = n_egalitarian + n_wealth_divided + n_stratified
      )
    ) +
    scale_x_continuous(
      name = "Time before present (kya)",
      limits = c(min(time_slices$time), 0)
    ) +
    scale_y_continuous(
      name = "Number of lineages"
    ) +
    theme_classic()

  # plot ratio of stratified to egalitarian lineages
  right_plot <-
    time_slices |>
    unnest(c(n_egalitarian, n_wealth_divided, n_stratified)) |>
    mutate(ratio = n_stratified / n_egalitarian) |>
    group_by(time) |>
    summarise(
      estimate = median(ratio),
      lower25 = quantile(ratio, 0.375),
      upper25 = quantile(ratio, 0.625),
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
      fill = "grey97"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower50,
        ymax = upper50
      ),
      fill = "grey90"
    ) +
    geom_ribbon(
      aes(
        x = time,
        ymin = lower25,
        ymax = upper25
      ),
      fill = "grey80"
    ) +
    geom_line(
      aes(
        x = time,
        y = estimate
      )
    ) +
    scale_x_continuous(
      name = "Time before present (kya)",
      limits = c(min(time_slices$time), 0)
    ) +
    scale_y_continuous(
      name = "Ratio of stratified lineages\nto egalitarian lineages",
      limits = c(0, 1.02),
      breaks = seq(0, 1, by = 0.25)
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
