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

get_model_flags <- function(constants) {
	out <- list()
	out$single_property <- constants$n_property == 1
	out$single_method <- constants$n_method == 1
	out$use_shooting <- any(constants$shooting == 1)

	use_snares <- any(constants$method == 4)
	use_traps <- any(constants$method == 5)
	out$use_traps_and_snares <- use_traps & use_snares
	out$use_traps_or_snares <- use_traps + use_snares == 1

	out
}
