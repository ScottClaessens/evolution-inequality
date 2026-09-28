#' Produce table of log Bayes Factors for model comparison
#'
#' @param fits Model outputs combined with \code{bind_rows()}
#'
#' @returns Tibble
#'
get_table_model_comparison <- function(fits) {

  # get table of log marginal likelihoods
  log_liks <-
    fits |>
    filter(chain == 1) |>
    mutate(
      model = factor(model, levels = c("full", "unilinear", "relaxed_unilinear",
                                       "rectilinear"))
    ) |>
    group_by(model) |>
    summarise(log_lik = unique(log_lik))

  # function to calculate log bayes factors
  calculate_log_bfs <- function(model) {
    2 * (log_liks$log_lik - log_liks$log_lik[log_liks$model == model])
  }

  # calculate log bayes factors
  log_liks$full              <- calculate_log_bfs(model = "full")
  log_liks$rectilinear       <- calculate_log_bfs(model = "rectilinear")
  log_liks$unilinear         <- calculate_log_bfs(model = "unilinear")
  log_liks$relaxed_unilinear <- calculate_log_bfs(model = "relaxed_unilinear")

  # function to print log bayes factors
  print_log_bf <- function(x) {
    ifelse(x == 0, "-", format(round(x, 2), trim = TRUE, nsmall = 2))
  }

  # return table
  log_liks |>
    transmute(
      Model = str_to_sentence(str_replace(model, "_", " ")),
      `Log marginal likelihood` = format(round(log_lik, 2), nsmall = 2),
      `Full` = print_log_bf(full),
      `Unilinear` = print_log_bf(unilinear),
      `Relaxed unilinear` = print_log_bf(relaxed_unilinear),
      `Rectilinear` = print_log_bf(rectilinear)
    )

}
