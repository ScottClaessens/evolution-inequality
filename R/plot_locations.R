#' Plot locations of ancestral lineages from phylogeography
#'
#' @param locations Tibble of estimated longitude-latitude locations of internal
#'   nodes from the global language tree
#' @param ancestral_states Tibble of ancestral states
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param tree_id Indexes for trees
#' @param time Numeric. Time slice to use in plot.
#'
#' @returns List. A grob from a ggplot object and its legend
#'
plot_locations <- function(locations, ancestral_states, data,
                           tree, tree_id, time) {

  # get active lineages at this time over all trees
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
  rm(ancestral_states)

  # get shared taxa
  shared_taxa <- intersect(
    unique(unlist(locations$descendants)),
    unique(data$glottocode)
  )

  # get glottocode descendants from active lineages
  active_lineages <-
    active_lineages |>
    rowwise() |>
    mutate(
      # get posterior summaries for plot
      p1 = median(p1),
      p2 = median(p2),
      p3 = 1 - (p1 + p2),
      probably_egalitarian = p1 > p2 & p1 > p3,
      probably_wealth_divided = p2 > p1 & p2 > p3,
      probably_stratified = p3 > p1 & p3 > p2,
      # variable for plotting
      var_for_plot = list(
        c(probably_egalitarian,
          probably_wealth_divided,
          probably_stratified)
      ),
      # get descendants as integer indexes
      descendants = phangorn::Descendants(
        x = tree[[tree_id]],
        node = parent_node,
        type = "tips"
      ),
      # convert to xd_id tip labels
      descendants = list(
        tree[[tree_id]]$tip.label[descendants]
      ),
      # convert to glottocodes
      descendants = list(
        sort(unique(data$glottocode[match(descendants, data$xd_id)]))
      ),
      # only shared taxa
      descendants = list(
        descendants[descendants %in% shared_taxa]
      )
    ) |>
    ungroup()

  # wrangle locations data to include only descendants that are shared taxa
  locations <-
    locations |>
    rowwise() |>
    mutate(
      descendants = list(
        sort(unique(descendants[descendants %in% shared_taxa]))
      )
    ) |>
    ungroup() |>
    dplyr::select(!c(node)) |>
    # get average locations to avoid duplicating lineages
    group_by(descendants) |>
    mutate(
      latitude = mean(latitude),
      longitude = mean(longitude)
    )

  # join datasets by nodes that share the same glottocode descendants
  locations <-
    left_join(
      x = active_lineages,
      y = locations,
      by = c("tree_id", "descendants")
    ) |>
    dplyr::select(!descendants) |>
    # drop cases without location data
    drop_na(c(latitude, longitude))

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
        x = ((longitude + 180) %% 360) - 180,
        y = latitude,
        z = var_for_plot
      ),
      fun = function(x) {
        num_egalitarian    <- sum(unlist(lapply(x, function(i) i[1])))
        num_wealth_divided <- sum(unlist(lapply(x, function(i) i[2])))
        num_stratified     <- sum(unlist(lapply(x, function(i) i[3])))
        log((num_stratified + 1) / (num_egalitarian + 1)) # log ratio (offset)
      },
      colour = NA,
      binwidth = 10,
      center = c(0, 0)
    ) +
    geom_sf(
      fill = NA,
      colour = "black",
      linewidth = 0.07
    ) +
    scale_fill_gradient2(
      name = "Log ratio of stratified lineages\nto egalitarian lineages",
      low = "#19AAF1",
      mid = "white",
      high = "#E6550E",
      na.value = "green",
      breaks = c(-5, 0, 5),
      limits = c(-7.5, 7.5)
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
        fill = "grey90"
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.4
      ),
      plot.margin = margin(0, 5, 0, 5),
      legend.position = "bottom"
    )

  # convert to grob and legend to save space
  out <-
    list(
      grob = ggplot2::ggplotGrob(out + theme(legend.position = "none")),
      legend = ggpubr::get_legend(out)
    )

  # clean up
  rm(locations, active_lineages, data, tree, shared_taxa,
     time, tree_id, bins, breaks, limits)

  # return
  out

}
