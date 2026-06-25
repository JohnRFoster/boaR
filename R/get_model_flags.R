#' Create boolean values to build the correct model
#' @param constants The constants list the is used to build a nimble model.
#' @return a list of values
#' 	single_property is TRUE of only one property is in the data.
#' 	single_method is TRUE of only one removal method is used.
#' 	use_shooting is TRUE if firearms, fixed-wing aircraft, or helicopters are used.
#' 	use_traps_and_snares is TRUE if both traps and snares are used.
#' 	use_traps_or_snares is TRUE if one of traps or snares are used, but not both.
#'
#' @export

get_model_flags <- function(df) {
  out <- list()
  out$single_property <- length(unique(df$propertyID)) == 1
  out$single_method <- length(unique(df$method)) == 1

  use_fire <- any(df$method == "FIREARMS")
  use_fixed <- any(df$method == "FIXED WING")
  use_heli <- any(df$method == "HELICOPTER")
  use_snares <- any(df$method == "SNARE")
  use_traps <- any(df$method == "TRAPS")

  out$use_shooting <- any(use_fire | use_fixed | use_heli)
  out$use_traps_and_snares <- use_traps & use_snares
  out$use_traps_or_snares <- use_traps + use_snares == 1
  out$no_landcover <- all(is.na(df$c_canopy))
  out$use_beta_p <- !out$no_landcover

  out
}
