options(tidyverse.quiet = TRUE)
library(crew)
library(targets)
library(tarchetypes)
library(tidyverse)

tar_option_set(
  packages = c("ape", "ggtree", "patchwork", "phangorn", "phytools",
               "rnaturalearth", "tidyverse", "withr"),
  controller = crew_controller_local(workers = 8),
  deployment = "main"
)
tar_source()

# pipeline
list(

  # ─────────────────────────────────────────
  # Load D-PLACE data and phylogeny
  # ─────────────────────────────────────────

  # get data urls
  tar_target(
    dplace_data_url,
    paste0(
      "https://raw.githubusercontent.com/D-PLACE/dplace-cldf/",
      "6c2008c187a297d1955b41d8ae80d8e31d404f6c/cldf/data.csv"
    ),
    format = "url"
  ),
  tar_target(
    dplace_societies_url,
    paste0(
      "https://raw.githubusercontent.com/D-PLACE/dplace-cldf/",
      "6c2008c187a297d1955b41d8ae80d8e31d404f6c/cldf/societies.csv"
    ),
    format = "url"
  ),
  tar_target(
    glottolog_languages_url,
    paste0(
      "https://raw.githubusercontent.com/glottolog/glottolog-cldf/",
      "072ca0d0410039fb8b779be8fc165bac575d2cda/cldf/languages.csv"
    ),
    format = "url"
  ),

  # get tree file path
  tar_target(tree_file, "data/tree/dplace.nxs", format = "file"),

  # load tree
  tar_target(tree, read.nexus(tree_file)),

  # compute maximum clade credibility tree
  tar_target(mcc_tree, phangorn::mcc(tree)),

  # load dplace data
  tar_target(
    data,
    load_dplace_data(
      dplace_data_url, dplace_societies_url,
      glottolog_languages_url, mcc_tree
    )
  ),

  # ─────────────────────────────────────────
  # Plot world map
  # ─────────────────────────────────────────

  # plot world map
  tar_target(plot_world_map, plot_map(data)),

  # ─────────────────────────────────────────
  # Calculate phylogenetic signal
  # ─────────────────────────────────────────

  # calculate phylogenetic signal
  tar_target(
    phylogenetic_signal,
    calculate_phylogenetic_signal(data, mcc_tree)
  ),

  # ─────────────────────────────────────────
  # Compare models of evolution for full tree
  # ─────────────────────────────────────────

  # get independent mcmc chains
  tar_target(chain, 1:4),

  # loop over models
  tar_map(

    values = tibble(
      model = c("full", "rectilinear", "unilinear", "relaxed_unilinear")
    ),

    # fit model
    tar_target(
      fit,
      fit_model(data, tree, chain, model),
      pattern = map(chain),
      deployment = "worker",
      storage = "worker",
      retrieval = "worker"
    ),

    # get diagnostics
    tar_target(diagnostics, calculate_model_diagnostics(fit)),

    # plot MCMC trace
    tar_target(plot_trace, plot_mcmc_trace(fit, model))

  ),

  # model comparison table
  tar_target(
    table_model_comparison,
    get_table_model_comparison(
      bind_rows(
        fit_full, fit_rectilinear, fit_unilinear, fit_relaxed_unilinear
      )
    )
  ),

  # ─────────────────────────────────────────────────────
  # Compare models of evolution in specific world regions
  # ─────────────────────────────────────────────────────

  # loop over world regions
  tar_map(

    values = tibble(
      subset_region = c("Africa", "Americas", "Eurasia", "Oceania")
    ),

    # loop over models
    tar_map(

      values = tibble(
        model = c("full", "rectilinear", "unilinear", "relaxed_unilinear")
      ),

      # fit model
      tar_target(
        fit,
        fit_model(data, tree, chain, model,
                  subset_region = subset_region),
        pattern = map(chain),
        deployment = "worker",
        storage = "worker",
        retrieval = "worker"
      ),

      # get diagnostics
      tar_target(diagnostics, calculate_model_diagnostics(fit))

    ),

    # model comparison table
    tar_target(
      table_model_comparison,
      get_table_model_comparison(
        bind_rows(
          fit_full, fit_rectilinear, fit_unilinear, fit_relaxed_unilinear
        )
      )
    )

  ),

  # ─────────────────────────────────────────
  # Estimate ancestral states
  # ─────────────────────────────────────────

  # get tree ids
  tar_target(tree_id, sample(1:length(tree), size = 100, replace = FALSE)),

  # fit ancestral state reconstruction model
  tar_target(
    fit_asr,
    fit_model(data, tree, chain, model = "relaxed_unilinear",
              iter = 550000, burnin = 50000, stones = FALSE,
              asr = TRUE, tree_id = tree_id),
    pattern = cross(chain, tree_id),
    # run in parallel
    deployment = "worker",
    storage = "worker",
    retrieval = "worker"
  ),

  # calculate model diagnostics
  tar_target(diagnostics_asr, calculate_model_diagnostics(fit_asr)),

  # plot tree
  tar_target(
    plot_tree_states,
    plot_tree(data, tree, tree_id, fit_asr)
  ),

  # ─────────────────────────────────────────
  # Produce manuscript
  # ─────────────────────────────────────────

  # knit quarto manuscript
  tar_quarto(manuscript, "quarto/manuscript.qmd", quiet = FALSE),

  # write session info
  tar_target(
    sessionInfo,
    writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  )

)
