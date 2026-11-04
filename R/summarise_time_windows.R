#' Calculate and summarise time window calculations over time
#'
#' @param ancestral_states Tibble of ancestral states
#' @param start Numeric. Start time for the time window calculations (in ky).
#' @param end Numeric. End time for the time window calculations (in ky).
#' @param by Numeric. Increments from start time to end time (in ky).
#' @param window Numeric. Width of the time window for calculations (in ky).
#'
#' @returns Tibble of time window calculations
#'
summarise_time_windows <- function(ancestral_states, start = -20, end = 0,
                                   by = 0.1, window = 1) {

  # run time step calculations
  out <- lapply(seq(start, end, by = by), function(t) {

    # get nodes in time window
    ancestral_states |>
      filter(time >= (t - window) & time < t) |>

      # sum across all nodes, retaining posterior uncertainty
      rowid_to_column() |>
      rowwise() |>
      mutate(drawid = list(1:length(probably_egalitarian))) |>
      unnest(c(drawid, starts_with("probably_"), p1, p2, p3)) |>
      group_by(drawid) |>
      summarise(
        n = n(),
        max_dist_from_equality_p1 = max(abs(p1 - (1/3))),
        max_dist_from_equality_p2 = max(abs(p2 - (1/3))),
        max_dist_from_equality_p3 = max(abs(p3 - (1/3))),
        n_egalitarian = sum(probably_egalitarian),
        n_wealth_divided = sum(probably_wealth_divided),
        n_stratified = sum(probably_stratified),
        ratio = n_stratified / n_egalitarian
      ) |>

      # return posterior samples as lists
      summarise(
        time                      = t,
        n                         = unique(n),
        max_dist_from_equality_p1 = list(max_dist_from_equality_p1),
        max_dist_from_equality_p2 = list(max_dist_from_equality_p2),
        max_dist_from_equality_p3 = list(max_dist_from_equality_p3),
        n_egalitarian             = list(n_egalitarian),
        n_wealth_divided          = list(n_wealth_divided),
        n_stratified              = list(n_stratified),
        ratio                     = list(ratio)
      )

  })

  # bind rows and return
  do.call(rbind, out)

}
