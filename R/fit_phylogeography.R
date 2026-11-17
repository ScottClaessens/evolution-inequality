#' Fit phylogeography model with BEAST
#'
#' This function fits a phylogeography model to the D-PLACE phylogeny given
#' latitude and longitude locations for societies. The specification for the
#' model is outlined in "beast/phylogeography.xml". We assume geographic
#' diffusion on a sphere with a strict clock diffusion rate. Output files for
#' the MCMC sampler are outputted to the beast/ directory and then cleaned up
#' by the function. The function returns a tibble of post-burn-in posterior
#' draws for internal node coordinates (longitude and latitude values) alongside
#' convergence statistics (Rhat and effective sample sizes).

#' Note: this function requires BEAST 2.7.7 to be installed and has only been
#' tested on a Macbook Pro, so it may need adapting for other operating systems.
#'
#' @param data Tibble of D-PLACE data
#' @param tree Tree object of class multiPhylo
#' @param tree_id Numeric. Index for tree.
#' @param iter Numeric. Number of MCMC sampling iterations.
#' @param burnin Numeric. Number of MCMC burn-in iterations.
#' @param thin Numeric. Number of MCMC iterations to thin in the output.
#' @param seed Numeric. Random seed for the model.
#'
#' @returns A tibble of posterior samples
#'
fit_phylogeography <- function(data, tree, tree_id = 1, iter = 550000,
                               burnin = 50000, thin = 1000, seed = 1) {

  # get tree
  tree <- tree[[tree_id]]

  # reorder data to align with tip labels
  data <- data[match(tree$tip.label, data$xd_id), ]

  # add lon-lat information to tip names for beast
  tree$tip.label <-
    paste(
      tree$tip.label,
      data$latitude,
      data$longitude,
      sep = "_"
    )

  # get newick string for this tree
  newick_string <- write.tree(tree)

  # copy master xml file with added newick string and sampling options
  readLines("beast/phylogeography.xml") |>
    str_replace(
      pattern = 'newick=""',
      replacement = paste0('newick="', newick_string, '"')
    ) |>
    str_replace(
      pattern = 'chainLength=""',
      replacement = paste0('chainLength="', as.integer(iter), '"')
    ) |>
    str_replace(
      pattern = 'logEvery=""',
      replacement = paste0('logEvery="', as.integer(thin), '"')
    ) |>
    writeLines(paste0("beast/phylogeography_", tree_id, ".xml"))

  # run beast in command line
  dir <- paste0(here::here(), "/beast")
  invisible(
    system(
      paste0(
        "cd ", dir, "; ",
        "/Applications/BEAST\\ 2.7.7/bin/beast ",
        "-seed ", seed, " ",
        "-threads 1 ",
        "phylogeography_", tree_id, ".xml"
      ),
      intern = TRUE,
      ignore.stdout = TRUE,
      ignore.stderr = TRUE
    )
  )

  # free up memory
  rm(data, tree)

  # load results
  results <- suppressWarnings(
    treeio::read.beast(
      file = paste0("beast/phylogeography_", tree_id, ".trees")
    )
  )

  # cleanup files in beast directory
  invisible(file.remove(paste0("beast/phylogeography_", tree_id, ".xml")))
  invisible(file.remove(paste0("beast/phylogeography_", tree_id, ".xml.state")))
  invisible(file.remove(paste0("beast/phylogeography_", tree_id, ".log")))
  invisible(file.remove(paste0("beast/phylogeography_", tree_id, ".trees")))

  # get ids for post-burn-in draws
  draw_idx <- seq_along(results)[parse_number(names(results)) > burnin]

  # collect results for internal nodes
  out <-
    lapply(draw_idx, function(i) {
      results[[i]]@data |>
        transmute(
          tree_id = tree_id,
          node,
          locationsgeo
        ) |>
        unnest_wider(
          locationsgeo,
          names_sep = "_"
        ) |>
        rename(
          latitude = locationsgeo_1,
          longitude = locationsgeo_2
        ) |>
        filter_out(node %in% seq_along(results[[i]]@phylo$tip.label))
    })

  # convert to tibble and calculate convergence statistics
  do.call(rbind, out) |>
    mutate(node = as.numeric(node)) |>
    group_by(tree_id, node) |>
    summarise(
      latitude = list(as.numeric(latitude)),
      longitude = list(as.numeric(longitude)),
      .groups = "drop"
    ) |>
    pivot_longer(
      cols = c(latitude, longitude),
      names_to = "dimension",
      values_to = "post"
    ) |>
    arrange(node) |>
    rowwise() |>
    mutate(
      rhat = rstan::Rhat(post),
      ess_bulk = rstan::ess_bulk(post),
      ess_tail = rstan::ess_tail(post)
    ) |>
    ungroup()

}
