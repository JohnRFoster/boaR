#' Prepare data for nimble model
#' @param df A data frame containing the survey data.
#' @return A list of data for the nimble model.
#' @export

nimble_data <- function(df) {
	X <- create_X(df)
	y_sum <- removed_in_pp_cumsum(df)

	list(
		y = df$take,
		y_sum = y_sum,
		J = data_ls, # mean litter size year from VerCauteren et al. 2019 pg 63
		X_p = X,
		effort_per = df$effort_per,
		log_effort_per = log(df$effort_per),
		n_trap_m1 = df$trap_count - 1,
		log_survey_area_km2 = log(df$property_area_km2)
	)
}
