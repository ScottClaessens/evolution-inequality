#' Combine phylogeography maps into a single plot
#'
#' @param plot_geo List of ggplot objects
#'
#' @returns A patchwork of ggplot objects
#'
combine_phylogeography_plots <- function(plot_geo) {

  # get list of ggplot objects
  plot_list <- lapply(plot_geo, function(x) ggpubr::as_ggplot(x[["grob"]]))

  # combine plots
  out <-
    wrap_plots(
      plot_list,
      guides = "collect"
    ) +
      ggpubr::as_ggplot(plot_geo[[1]]$legend) +
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
    filename = "plots/geography.pdf",
    height = 6,
    width = 6
  )

  # return
  out

}
