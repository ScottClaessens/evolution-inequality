#' Plot results of phylogeographic modelling
#'
#' @param fit_geo Tibble of posterior samples from the fitted model
#' @param ancestral_states Tibble of ancestral states
#' @param tree Tree object of class multiPhylo
#' @param tree_id Indexes for trees
#' @param time Numeric. Time slice to use in plot.
#' @param breaks Numeric vector. Breaks for the fill gradient.
#' @param limits Numeric vector of length 2. Lower and upper limits for the fill
#'   gradient.
#'
#' @returns A ggplot object
#'
plot_phylogeography <- function(fit_geo, ancestral_states, tree, tree_id,
                                time = -20, breaks = c(0, 10, 100),
                                limits = c(0, 400)) {

  # summarise ancestral states
  ancestral_states <-
    ancestral_states |>
    rowwise() |>
    transmute(
      tree_id, node,
      p1 = median(p1),
      p2 = median(p2),
      p3 = 1 - (p1 + p2),
      egalitarian = p1 > p2 & p1 > p3,
      wealth_divided = p2 > p1 & p2 > p3,
      stratified = p3 > p1 & p3 > p2,
      var_for_plot = list(c(egalitarian, wealth_divided, stratified))
    ) |>
    ungroup()

  # get active lineages over all trees
  active_lineages <- lapply(seq_along(tree_id), function(i) {

    # current tree id
    id <- tree_id[i]

    # get branch ages
    times <- node.depth.edgelength(tree[[id]])
    ages <- times - max(times)

    # get lineages
    edges <- tree[[id]]$edge
    parent_nodes <- edges[, 1]
    child_nodes <- edges[, 2]
    start_times <- ages[parent_nodes]
    end_times <- ages[child_nodes]

    # keep only active lineages
    idx <- which(start_times <= time & end_times >= time)

    # get active lineages
    tibble(
      tree_id = id,
      parent_node = parent_nodes[idx],
      child_node = child_nodes[idx],
      parent_time = start_times[idx],
      child_time = end_times[idx]
    )

  })

  # join ancestral states by parent nodes
  active_lineages <-
    do.call(rbind, active_lineages) |>
    left_join(
      ancestral_states,
      by = c(
        "tree_id" = "tree_id",
        "parent_node" = "node"
      )
    )

  # clean up
  rm(ancestral_states, tree, tree_id)

  # join estimated ancestral locations
  locations <-
    fit_geo |>
    pivot_wider(
      names_from = dimension,
      values_from = post,
      id_cols = c(tree_id, node)
    ) |>
    right_join(
      active_lineages,
      by = c(
        "tree_id" = "tree_id",
        "node" = "parent_node"
      )
    ) |>
    unnest(c(longitude, latitude))

  # clean up
  rm(fit_geo, active_lineages)

  # plot world map
  out <-
    ne_coastline(
      scale = "small",
      returnclass = "sf"
    ) |>
    ggplot() +
    stat_summary_hex(
      data = locations,
      mapping = aes(
        x = longitude,
        y = latitude,
        z = var_for_plot
      ),
      fun = function(x) {
        sum(unlist(lapply(x, \(i) i[3]))) / 500
      },
      colour = NA,
      bins = 40
    ) +
    geom_sf(
      fill = NA,
      colour = "black",
      linewidth = 0.07
    ) +
    scale_fill_gradient(
      name = "Number of stratified lineages",
      low = "white",
      high = "#e6550e",
      na.value = "blue",
      trans = "log1p",
      breaks = breaks,
      limits = limits
    ) +
    coord_sf(
      xlim = c(-160, 160),
      ylim = c(-50, 80)
    ) +
    labs(
      title = paste0(
        ifelse(time == -0.1, 0, abs(time) * 1000),
        " BP"
      )
    ) +
    theme_void() +
    theme(
      panel.background = element_rect(
        fill = "white"
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.4
      ),
      plot.margin = margin(0, 5, 0, 5),
      legend.position = "bottom"
    )

  # clean up
  rm(locations)

  # return as grob and legend
  list(
    grob = ggplot2::ggplotGrob(out + theme(legend.position = "none")),
    legend = ggpubr::get_legend(out)
  )

}
