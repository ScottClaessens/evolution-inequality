#' Extract ancestral states from a fitted BayesTraits model
#'
#' @param fit_asr Tibble of posterior samples from the fitted model
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#'
#' @returns Tibble of ancestral states
#'
extract_ancestral_states <- function(fit_asr, data, tree) {

  # internal function to extract the time of a node
  get_node_time <- function(tree_id, node) {
    # get branch times
    times <- node.depth.edgelength(tree[[tree_id]])
    # return time before present (kya)
    times[node] - max(times)
  }

  # internal function to extract the world region of a node
  get_node_region <- function(tree_id, node) {
    # get tip labels
    tip_labels <- tree[[tree_id]]$tip.label
    # get node ids for descendant tips
    descendants <- Descendants(tree[[tree_id]], node, type = "tips")[[1]]
    # convert descendants to xd_id
    xd_ids <- tip_labels[descendants]
    # get world regions for descendants
    world_regions <- data$world_region[data$xd_id %in% xd_ids]
    # return result
    if (length(unique(world_regions)) == 1) {
      unique(world_regions)
    } else {
      NA
    }
  }

  # wrangle bayestraits output
  fit_asr |>
    select(
      tree_id,
      (ends_with("P(2)") | ends_with("P(3)")) & !starts_with("Root")
    ) |>

    # group posterior samples by individual trees
    group_by(tree_id) |>
    summarise(across(everything(), list)) |>

    # get in long format
    pivot_longer(
      cols = !tree_id,
      names_pattern = "N(.*) (.*)",
      names_to = c("node", ".value")
    ) |>

    # modify columns
    rowwise() |>
    transmute(
      tree_id = tree_id,
      node = as.numeric(node),
      time = get_node_time(tree_id, node),
      world_region = get_node_region(tree_id, node),
      # ancestral probabilities
      p1 = list(1 - (`P(2)` + `P(3)`)),
      p2 = list(`P(2)`),
      p3 = list(`P(3)`),
      # which state is most likely?
      probably_egalitarian    = list(p1 > p2 & p1 > p3),
      probably_wealth_divided = list(p2 > p1 & p2 > p3),
      probably_stratified     = list(p3 > p1 & p3 > p2)
    ) |>
    ungroup()

}
