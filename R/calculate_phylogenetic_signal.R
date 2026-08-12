#' Calculate phylogenetic signal (delta) for the class differentiation trait
#'
#' @param data Tibble of D-PLACE data
#' @param mcc_tree Tree object for maximum clade credibility tree
#'
#' @returns A list of results (delta and associated p-value)
#'
calculate_phylogenetic_signal <- function(data, mcc_tree) {

  # drop missing observations
  data <- drop_na(data, class_differentiation)

  # retain tips with observed data
  mcc_tree <- keep.tip(mcc_tree, tip = data$xd_id)

  # get trait (in line with phylogeny)
  trait <- data$class_differentiation[match(mcc_tree$tip.label, data$xd_id)]

  # convert trait to numeric
  trait <- ifelse(as.numeric(trait) %in% 3:5, 3, as.numeric(trait))

  # calculate delta (uses code in R/helper.R)
  delta_EA066 <- delta(
    trait = trait,   # trait data
    tree = mcc_tree, # maximum clade credibility tree
    lambda0 = 0.1,   # rate parameter of the proposal
    se = 0.5,        # standard deviation of the proposal
    sim = 10000,     # number of iterations
    thin = 10,       # thinning frequency
    burn = 100       # burn-in iterations
  )

  # calculate p-value
  random_delta <- rep(NA, 100)
  for (i in 1:100) {
    rtrait <- sample(trait)
    random_delta[i] <- delta(rtrait, mcc_tree, 0.1, 0.5, 10000, 10, 100)
  }
  p_value <- sum(random_delta > delta_EA066) / length(random_delta)

  # return
  list(
    delta = delta_EA066,
    p_value = p_value
  )

}
