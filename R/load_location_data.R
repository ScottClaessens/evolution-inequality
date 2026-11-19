#' Load longitude and latitude data for internal nodes
#'
#' @param locations_file File path to location data
#' @param tree_id Indexes for trees
#'
#' @returns A tibble
#'
load_location_data <- function(locations_file, tree_id) {

  # get all trees with associated location metadata
  trees <-
    treeio::read.beast(
      file = locations_file,
      threads = 4
    )

  # subset to fewer trees to reduce size
  trees <- trees[tree_id]

  # for each tree, return tibble of node locations
  out <-
    parallel::mclapply(
      X = seq_along(tree_id),
      mc.cores = 2,
      FUN = function(i) {

        # extract data and phylo objects
        data <- trees[[i]]@data
        phylo <- trees[[i]]@phylo

        # filter to internal nodes only
        data |>
          filter_out(
            node %in% seq_along(phylo$tip.label)
          ) |>

          # extract longitude and latitude
          unnest_wider(
            col = location,
            names_sep = "_"
          ) |>
          rename(
            latitude = location_1,
            longitude = location_2
          ) |>

          # tidy
          mutate(tree_id = tree_id[i]) |>
          dplyr::select(tree_id, node, everything()) |>

          # get descendant glottologs
          rowwise() |>
          mutate(
            descendants = phangorn::Descendants(
              x = phylo,
              node = node,
              type = "tips"
            ),
            descendants = list(phylo$tip.label[descendants])
          ) |>
          ungroup()

      }
    )

  # return as tibble
  do.call(bind_rows, out)

}
