#' Combine location maps into a single plot
#'
#' @param plot_locations List of ggplot grobs
#'
#' @returns A patchwork of ggplot objects
#'
combine_location_plots <- function(plot_location) {

  # get list of ggplot objects
  plot_list <- lapply(plot_location, \(x) ggpubr::as_ggplot(x[["grob"]]))

  # combine plots
  out <-
    wrap_plots(
      plot_list,
      guides = "collect"
    ) +
    ggpubr::as_ggplot(plot_location[[1]]$legend) +
    patchwork::plot_layout(
      heights = c(1, 1, 1, 0.5),
      design = '
          AB
          CD
          EF
          GG
        '
    )

  # save
  ggsave(
    plot = out,
    filename = "plots/locations.pdf",
    height = 6,
    width = 6
  )

  # return
  out

}
