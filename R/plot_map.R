#' Plot class differentiation data on world map
#'
#' @param data Tibble of D-PLACE data
#'
#' @returns A ggplot object
#'
plot_map <- function(data) {

  # plot
  out <-
    ne_countries(
      scale = "small",
      returnclass = "sf"
    ) |>
    filter(continent != "Antarctica") |>
    ggplot() +
    geom_sf(
      fill = "grey80",
      colour = NA
    ) +
    geom_point(
      data =
        data |>
        mutate(
          class = as.numeric(class_differentiation),
          class = case_when(
            class == 1 ~ "Absence of distinctions",
            class == 2 ~ "Wealth distinctions",
            class %in% 3:5 ~ "Stratification",
            is.na(class) ~ NA
          ),
          class = factor(
            class,
            levels = c(
              "Absence of distinctions",
              "Wealth distinctions",
              "Stratification"
            ),
            ordered = TRUE
          )
        ) |>
        arrange(class),
      mapping = aes(
        x = longitude,
        y = latitude,
        fill = class
      ),
      size = 1.2,
      shape = 21
    ) +
    scale_fill_brewer(
      type = "seq",
      palette = 7,
      guide = guide_legend(
        override.aes = list(
          size = 4,
          shape = "square filled"
        )
      )
    ) +
    theme_void() +
    theme(
      legend.title = element_blank(),
      legend.key.spacing.y = unit(-2, "mm"),
      legend.position = "inside",
      legend.position.inside = c(0.12, 0.04),
      legend.box.background = element_rect(),
      legend.margin = margin(0, 5, 0, 0)
    )

  # save
  ggsave(
    filename = "plots/world_map.pdf",
    plot = out,
    height = 4,
    width = 8
  )

  # return
  out

}
