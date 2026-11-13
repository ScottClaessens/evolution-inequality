#' Calculate and summarise time slice calculations over time
#'
#' @param ancestral_states Tibble of ancestral states
#' @param tree Tree object of class multiPhylo
#' @param tree_id Index for tree
#' @param time Numeric. Time for the time slice calculations (in ky before
#'   present).
#'
#' @returns Tibble of time slice calculations
#'
summarise_time_slices <- function(ancestral_states, tree, tree_id, time) {

  # get tree
  tree <- tree[[tree_id]]

  # get branch ages
  times <- node.depth.edgelength(tree)
  ages <- times - max(times)

  # get lineages
  edges <- tree$edge
  parent_nodes <- edges[, 1]
  child_nodes <- edges[, 2]
  start_times <- ages[parent_nodes]
  end_times <- ages[child_nodes]

  # keep only active lineages
  idx <- which(start_times <= time & end_times >= time)

  # get active lineages
  active_lineages <-
    tibble(
      tree_id = tree_id,
      parent_node = parent_nodes[idx],
      child_node = child_nodes[idx],
      parent_time = start_times[idx],
      child_time = end_times[idx]
    ) |>

    # join ancestral states by parent nodes
    left_join(
      dplyr::select(ancestral_states, !c(time, world_region)),
      by = c(
        "tree_id" = "tree_id",
        "parent_node" = "node"
      )
    )

  # clean up to save space
  rm(ancestral_states, tree)

  # unnest posterior draws
  active_lineages |>
    rowid_to_column() |>
    rowwise() |>
    mutate(drawid = list(seq_along(probably_egalitarian))) |>
    unnest(c(drawid, starts_with("probably_"), p1, p2, p3)) |>

    # sum across all nodes, retaining posterior uncertainty
    group_by(drawid) |>
    summarise(
      n_active_lineages = n(),
      n_egalitarian = sum(probably_egalitarian),
      n_wealth_divided = sum(probably_wealth_divided),
      n_stratified = sum(probably_stratified)
    ) |>

    # return posterior samples as lists
    summarise(
      tree_id           = tree_id,
      time              = time,
      n_active_lineages = unique(n_active_lineages),
      n_egalitarian     = list(n_egalitarian),
      n_wealth_divided  = list(n_wealth_divided),
      n_stratified      = list(n_stratified)
    )

}
