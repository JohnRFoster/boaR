#' Read MIS data, filter, and create primary periods
#' @param file The file path to the raw data CSV.
#' @param interval The number of days in each primary period.
#' @param create_new Whether to create new primary periods or use existing ones (from calibration).
#'
#'
#' @export

get_data <- function(file, interval, create_new) {
	all_take <- readr::read_csv(file, show_col_types = FALSE) |>
		dplyr::filter(start.date >= lubridate::ymd("2014-01-01")) |>
		dplyr::mutate(
			cnty_name = if_else(
				grepl("ST ", cnty_name),
				gsub("ST ", "ST. ", cnty_name),
				cnty_name
			),
			cnty_name = if_else(grepl("KERN", cnty_name), "KERN", cnty_name)
		)

	data_mis <- all_take |>
		dplyr::mutate(property_area_km2 = round(property.size * 0.00404686, 2)) |>
		dplyr::filter(property_area_km2 >= 1.8, st_name != "HAWAII") |>
		dplyr::mutate(
			effort = if_else(
				cmp_name %in% c("TRAPS, CAGE", "SNARE"),
				cmp.days,
				cmp.hours
			),
			effort_per = effort / cmp.qty,
			cmp_name = if_else(cmp_name == "TRAPS, CAGE", "TRAPS", cmp_name)
		) |>
		dplyr::rename(method = cmp_name, trap_count = cmp.qty) |>
		dplyr::select(-wt_work_date, -hours, -cmp.hours, -cmp.days) |>
		dplyr::distinct() |>
		dplyr::mutate(propertyID = paste0(agrp_prp_id, "-", alws_agrprop_id)) |>
		dplyr::arrange(propertyID, start.date, end.date)

	# create PP of length [interval]
	data_timestep <- create_primary_periods(data_mis, interval, create_new) |>
		resolve_duplicate() |> # resolve duplicate property areas
		take_filter() |> # remove properties with zero pigs taken
		dynamic_filter() |> # filter out bad events & properties
		condition_first_capture() |> # condition on first positive removal event for each property
		dynamic_filter() |> # filter out bad events & properties
		order_interval() |> # determine midpoints from start/end dates
		order_stochastic() |> # randomly order events
		order_of_events() |> # assign order number, check
		county_codes() # county codes and renaming

	# now we have two columns for time
	# primary_period is how [interval] sequences are aligned across the data set
	# observed_timestep is the sequence of primary periods with
	# removal events within a property
	timestep_df <- create_timestep_df(data_timestep)

	data_mis <- dplyr::left_join(
		data_timestep,
		timestep_df,
		by = join_by(propertyID, primary_period)
	) |>
		dplyr::mutate(primary_period = primary_period - min(primary_period) + 1)

	## observation covariates ----
	# TODO get more (and updated) observation covariates
	data_obs <- get_obs_covars(land_csv)

	## join MIS with observation covariates ----
	data_join <- dplyr::left_join(data_mis, data_obs, by = join_by(county_code))
	data_join |>
		dplyr::filter(!is.na(c_road_den)) |>
		create_ids()
}


# Deal with observation covariates, a small percentage of which are missing
get_obs_covars <- function(df) {
	df |>
		dplyr::mutate(county_code = sprintf("%05d", FIPS)) |>
		dplyr::select(-starts_with('sd'), -NAME, -FIPS) |>
		dplyr::group_by(STATE_NAME) |>
		dplyr::mutate(
			rural.road.density = mean_impute(rural.road.density),
			prop.pub.land = mean_impute(prop.pub.land),
			mean.ruggedness = mean_impute(mean.ruggedness),
			mean.canopy.density = mean_impute(mean.canopy.density)
		) |>
		dplyr::ungroup() |>
		dplyr::mutate(
			c_road_den = center_scale(rural.road.density),
			c_rugged = center_scale(mean.ruggedness),
			c_canopy = center_scale(mean.canopy.density)
		)
}
