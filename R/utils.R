# utility functions for boaR
# mostly used for processing data and MCMC samples, creating nimble lists

# sample for the specified duration 'n_iters'
# drop previous samples (we save them after each time this function is called, so we don't need them)
continue_sampling <- function() {
	Cmcmc$run(niter = n_iters, reset = FALSE, resetMV = TRUE, thin = 10)
	samples <<- as.matrix(Cmcmc$mvSamples) # need to define this in the global environment to use in subset_mcmc
	samples
}

# for taking nodes and splitting them from an mcmc list
# as is, does one node at a time but can be passed to lapply to use on a vector of nodes
subset_mcmc <- function(node) {
	s <- samples[, grep(node, colnames(samples), value = TRUE, fixed = TRUE)] |>
		as.matrix() |>
		tibble::as_tibble()
	if (ncol(s) == 1) {
		colnames(s) <- node
	}
	s
}

# get parameter nodes using subset_mcmc
subset_params <- function() {
	purrr::map_dfc(
		lapply(params_check, subset_mcmc),
		tibble::as_tibble
	) |>
		as.matrix()
}

# get observed abundance nodes using subset_mcmc
subset_N_observed <- function() {
	nodes <- paste0("N[", model_constants$N_full_unique, "]")
	purrr::map_dfc(
		lapply(nodes, subset_mcmc),
		tibble::as_tibble
	) |>
		as.matrix()
}

# get unobserved abundance nodes using subset_mcmc
subset_N_unobserved <- function() {
	nodes <- paste0("N[", model_constants$N_quant_unique, "]")
	N_unobserved <- purrr::map_dfc(
		lapply(nodes, subset_mcmc),
		tibble::as_tibble
	) |>
		as.matrix()
	t(apply(N_unobserved, 2, quantile, c(0.025, 0.5, 0.975)))
}

write_out_p <- function(p, d, dest) {
	f <- "paramSamples.rds"
	readr::write_rds(list(params = p, diagnostic = d), file.path(dest, f))
}

write_abundance <- function(no, nu, dest) {
	f <- "observedAbundanceSamples.rds"
	readr::write_rds(no, file.path(dest, f))

	f <- "unobservedAbundanceQuantiles.rds"
	readr::write_rds(nu, file.path(dest, f))
}

# miscellaneous helper functions

# need all timesteps whether there are observations or not
create_all_primary_periods <- function(df) {
	pp_min_max <- df |>
		dplyr::select(property, primary_period) |>
		dplyr::distinct() |>
		dplyr::group_by(property) |>
		dplyr::filter(
			primary_period == min(primary_period) |
				primary_period == max(primary_period)
		) |>
		dplyr::ungroup()

	properties <- unique(pp_min_max$property)

	all_pp <- tibble::tibble()
	message("\nInclude all primary periods")

	if (!getOption("pbStyle", default = 3)) {
		# will be used if pbStyle is not set by the user (i.e. the default)
		pb <- utils::txtProgressBar(max = length(properties), style = 3)
	} else {
		# if the user calls: options(list(pbStyle = #))
		# the function will use style pbStyle.
		pb <- utils::txtProgressBar(max = length(properties), style = pbStyle)
	}

	for (i in seq_along(properties)) {
		pid <- pp_min_max |> dplyr::filter(property == properties[i])
		p_min <- min(pid$primary_period)
		p_max <- max(pid$primary_period)
		pp <- tibble::tibble(
			property = properties[i],
			primary_period = p_min:p_max,
			timestep = seq_along(p_min:p_max)
		)
		all_pp <- dplyr::bind_rows(all_pp, pp)
		utils::setTxtProgressBar(pb, i)
	}
	close(pb)
	all_pp |> dplyr::mutate(n_id = seq_len(n()))
}

# need to know the total number of timesteps in each property (sampled or not) for indexing
n_timesteps <- function(df) {
	df |>
		dplyr::group_by(property) |>
		dplyr::filter(timestep == max(timestep)) |>
		dplyr::pull(timestep)
}

# index (as a matrix) for tracking abundance in long format (converting wide to long)
# used in the process model
N_lookup_table <- function(df) {
	df |>
		dplyr::mutate(n_id = seq_len(n())) |>
		dplyr::select(-primary_period) |>
		tidyr::pivot_wider(names_from = timestep, values_from = n_id) |>
		dplyr::select(-property) |>
		as.matrix()
}

# calculate the cumulative number of pigs taken as a primary period progresses
removed_in_pp_cumsum <- function(df) {
	df |>
		dplyr::group_by(property, primary_period) |>
		dplyr::mutate(ysum = cumsum(take) - take) |>
		dplyr::ungroup() |>
		dplyr::pull(ysum)
}

# the total number of pigs taken in a primary period across all methods
# including periods without removals effort (equal to 0)
# wide format
total_take <- function(df_take, df_pp) {
	sum_take <- df_take |>
		dplyr::group_by(property, primary_period) |>
		dplyr::summarise(sum_take = sum(take)) |>
		dplyr::ungroup()

	dplyr::left_join(df_pp, sum_take, by = join_by(property, primary_period)) |>
		dplyr::mutate(sum_take = if_else(is.na(sum_take), 0, sum_take)) |>
		dplyr::select(-primary_period, -n_id) |>
		tidyr::pivot_wider(names_from = timestep, values_from = sum_take) |>
		dplyr::select(-property) |>
		as.matrix()
}

# index (long format) for which primary periods have removal effort
# and therefore included in the data model
N_lookup_data <- function(df_take, df_pp) {
	tH <- df_take |>
		dplyr::select(property, primary_period)

	df_pp |>
		dplyr::select(property, primary_period, n_id) |>
		dplyr::right_join(tH, by = join_by(property, primary_period)) |>
		dplyr::pull(n_id)
}

# need start and end indicies for data model
# used for estimiating p
create_start_end <- function(df_take, df_pp) {
	start <- end <- numeric(nrow(df_take))

	df <- dplyr::left_join(df_take, df_pp, by = join_by(primary_period, property))

	message("Creating start/end indicies")

	if (!getOption("pbStyle", default = 3)) {
		# will be used if pbStyle is not set by the user (i.e. the default)
		pb <- utils::txtProgressBar(max = nrow(df), style = 3)
	} else {
		# if the user calls: options(list(pbStyle = #))
		# the function will use style pbStyle.
		pb <- utils::txtProgressBar(max = nrow(df), style = pbStyle)
	}

	for (i in seq_len(nrow(df))) {
		if (df$order[i] > 1) {
			idx <- which(
				df$county == df$county[i] &
					df$property == df$property[i] &
					df$timestep == df$timestep[i] &
					df$order < df$order[i]
			)
			start[i] <- idx[1]
			end[i] <- idx[length(idx)]
			testthat::expect_equal(idx, start[i]:end[i])
		}
		utils::setTxtProgressBar(pb, i)
	}
	close(pb)
	tibble::tibble(start = start, end = end)
}

# informed hyper parameters for beta distribution on global pig survival
create_surv_prior <- function(interval) {
	data_usa <- vital_rate_data |>
		dplyr::filter(
			country == "USA",
			time.period.end != "null",
			time.period.start != "null",
			!paper.ID %in% c(128, 1007, 130, 136)
		) |> # these papers don't have specified date ranges or are meta-analysis
		mutate(
			time.period.end = mdy(time.period.end),
			time.period.start = mdy(time.period.start)
		)

	surv_data <- data_usa |>
		dplyr::filter(!is.na(survival.prop)) |>
		dplyr::select(
			unique.ID,
			paper.ID,
			N.hogs.in.study,
			contains("survival"),
			contains("hunting"),
			state,
			contains("time"),
			method.for.data
		)

	surv_mu <- surv_data |>
		dplyr::mutate(
			weeks = as.numeric(time.period.end - time.period.start) / 7,
			weeks4 = weeks / interval,
			survival.per.4week = survival.prop^(1 / weeks4),
			logit.survival.per.4week = boot::logit(survival.per.4week)
		) |>
		dplyr::filter(survival.per.4week > 0) |>
		dplyr::mutate(scale_factor = survival.per.4week / survival.prop)

	surv_mu_summary <- surv_mu |>
		dplyr::summarise(
			mu = mean(survival.per.4week),
			mu.logit = mean(logit.survival.per.4week)
		)

	surv_var <- surv_data |>
		dplyr::filter(survival.var.type %in% c("SD", "95% CI"))

	surv_sd <- surv_var |>
		dplyr::filter(survival.var.type == "SD") |>
		dplyr::mutate(sd = as.numeric(survival.var))

	surv_sd_calc <- surv_var |>
		dplyr::filter(survival.var.type == "95% CI") |>
		dplyr::mutate(
			low.CI = as.numeric(stringr::str_extract(
				survival.var,
				"[[:graph:]]*(?=\\-)"
			)),
			high.CI = as.numeric(stringr::str_extract(
				survival.var,
				"(?<=\\-)[[:graph:]]*"
			)),
			sd_low = (low.CI - survival.prop) / -1.96,
			sd_high = (high.CI - survival.prop) / 1.96
		) |>
		dplyr::group_by(unique.ID) |>
		dplyr::summarise(sd = max(sd_high, sd_low))

	surv_var_join <- dplyr::left_join(
		surv_var,
		surv_sd_calc,
		by = join_by(unique.ID)
	) |>
		dplyr::filter(survival.var.type != "SD")

	scale_ids <- surv_mu |>
		dplyr::select(unique.ID, scale_factor)

	surv_variance <- dplyr::bind_rows(surv_var_join, surv_sd) |>
		dplyr::left_join(scale_ids, by = join_by(unique.ID)) |>
		dplyr::mutate(
			variance = sd^2,
			variance.4week = variance * scale_factor^2,
			sd.4week = sqrt(variance.4week)
		)

	surv_sd_summary <- surv_variance |>
		dplyr::pull(sd.4week) |>
		mean()

	mu <- surv_mu_summary$mu
	psi <- 1 / mean(surv_variance$variance.4week)
	alpha <- mu * psi
	beta <- (1 - mu) * psi

	list(
		alpha = alpha,
		beta = beta
	)
}

# matrix of landscape covariates for data model
create_X <- function(df, cols = c("c_road_den", "c_rugged", "c_canopy")) {
	df |>
		dplyr::select(all_of(cols)) |>
		as.matrix()
}

get_prior_hyperparams <- function(post_round, posterior_path = NULL) {
	if (post_round == "first") {
		post <- post_first
	} else if (post_round == "last") {
		post <- post_last
	} else if (post_round == "create_new") {
		post <- read_rds(posterior_path)

		get_vec <- function(df, node) {
			df |>
				dplyr::select(contains(node)) |>
				tidyr::pivot_longer(cols = everything(), names_to = "node") |>
				dplyr::group_by(node) |>
				dplyr::summarise(mu = mean(value), tau = 1 / var(value))
		}

		log_rho <- get_vec(post, "log_rho")
		p_mu <- get_vec(post, "p_mu")
		log_gamma <- get_vec(post, "log_gamma")
		beta1 <- get_vec(post, "beta1")
		beta_p <- get_vec(post, "beta_p")
		log_nu <- get_vec(post, "log_nu")

		phi_mu <- get_vec(post, "phi_mu")
		mu <- phi_mu$mu
		v <- 1 / phi_mu$tau
		w <- ((mu * (1 - mu)) / v) - 1
		alpha <- mu * w
		beta <- (1 - mu) * w

		psi_phi <- get_vec(post, "psi_phi")
		mu <- psi_phi$mu
		v <- 1 / psi_phi$tau
		shape <- (mu^2) / v
		rate <- mu / v

		list(
			log_rho_mu = log_rho$mu,
			log_rho_tau = log_rho$tau,
			p_mu_mu = p_mu$mu,
			p_mu_tau = p_mu$tau,
			log_gamma_mu = log_gamma$mu,
			log_gamma_tau = log_gamma$tau,
			beta1_mu = beta1$mu,
			beta1_tau = beta1$tau,
			beta_p_mu = beta_p$mu,
			beta_p_tau = beta_p$tau,
			phi_mu_a = alpha,
			phi_mu_b = beta,
			psi_shape = shape,
			psi_rate = rate,
			log_nu_mu = log_nu$mu,
			log_nu_tau = log_nu$tau
		)
	}
}


# get max/min for prior for first latent state
# min = total take in the first primary period (abundance cannot be lower)
# max = total take assuming observation probability = 0.05
# if max is greater than a density of 100 pigs/km2 then set to
# max density of 100
# prior is on the log scale
get_n1_prior <- function(df) {
	df |>
		dplyr::filter(observed_timestep == 1) |>
		dplyr::group_by(property, property_area_km2) |>
		dplyr::reframe(total_take = sum(take)) |>
		dplyr::mutate(
			min_n = total_take + 5, # add a small buffer to allow for some unobserved pigs
			abundance_min = min_n,
			abundance_max = round(property_area_km2 * 100) + min_n
		)
}

create_ids <- function(df) {
	df |>
		dplyr::mutate(
			property = as.numeric(as.factor(propertyID)),
			county = as.numeric(as.factor(county_code)),
			ppID = paste0(propertyID, "-", primary_period),
			method_fac = as.numeric(as.factor(method)),
			eventID = paste0(ppID, "-", order, "-", method_fac)
		)
}


# figure out if we need to keep sampling
continue_mcmc <- function(mcmc, effective_size, max_psrf, verbose) {
	message("Checking convergence and sample size")
	psrf <- coda::gelman.diag(mcmc, multivariate = FALSE)$psrf
	effective_samples <- coda::effectiveSize(mcmc)

	converged <- all(psrf[, 2] < 1.1)
	enough_samples <- all(effective_samples >= effective_size)
	funky <- any(is.nan(psrf)) | max(psrf) > max_psrf

	message("Convergence [", converged, "]")
	message("Enough effective samples [", enough_samples, "]")
	message("Mixing [", !funky, "]")

	done <- converged & enough_samples
	if (done) {
		message("MCMC complete!")
	}

	# TODO determine burnin and effective sample size POST burnin

	if (funky) {
		message("\n*** Something is wrong with the mcmc! ***")
		done <- FALSE
	}

	if (verbose) {
		print(psrf)
		print(effective_samples)
	}

	list(
		done = done,
		psrf = psrf
	)
}

# for missing values, impute with state mean
mean_impute <- function(x) {
	ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}

# generate centered and scaled versions of these numeric variables
center_scale <- function(x) {
	(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

create_primary_periods <- function(df, interval, create_new = FALSE) {
	find_timestep <- function(dfpp, start, end) {
		after_start <- which(dfpp$start_dates <= start) |> max()
		before_end <- which(dfpp$end_dates >= end) |> min()
		if (after_start == before_end) {
			# then the start and end date is contained within a primary period
			timestep <- dfpp$timestep[before_end]
		} else {
			timestep <- NA
		} # otherwise, timestep[i] will be left as NA and filtered out later
	}

	message("Creating primary periods from start and end dates...")
	if (interval == 28 && !create_new) {
		timesteps <- purrr::map(
			seq_len(nrow(df)),
			~ find_timestep(primary_periods, df$start.date[.x], df$end.date[.x]),
			.progress = TRUE
		) |>
			unlist() |>
			as.integer()
	} else if (create_new) {
		end_dates <- unique(sort(df$end.date))
		min_date <- min(start_dates)
		max_date <- max(end_dates)

		start_dates <- seq(min_date, max_date, by = paste(interval, "day"))
		end_dates <- c(start_dates[-1] - 1, max_date)

		primary_periods <- tibble::tibble(start_dates, end_dates) |>
			mutate(timestep = seq_len(n()))
		primary_periods$month <- month(timestep_df$end_dates)
		primary_periods$year <- year(timestep_df$end_dates)
	} else {
		stop("Interval not supported")
	}

	timesteps <- purrr::map(
		seq_len(nrow(df)),
		~ find_timestep(primary_periods, df$start.date[.x], df$end.date[.x]),
		.progress = TRUE
	) |>
		unlist() |>
		as.integer()

	tmp <- df |>
		dplyr::mutate(timestep = timesteps)

	df_time <- dplyr::left_join(tmp, primary_periods)
	df_test <- df_time |> dplyr::filter(!is.na(timestep))

	testthat::expect_all_true(df_test$start.date >= df_test$start_dates)
	testthat::expect_all_true(df_test$end.date <= df_test$end_dates)

	df_test |>
		dplyr::arrange(propertyID, timestep)
}

# resolve duplicate property values - when there are multiple values, take max
resolve_duplicate <- function(insitu_data) {
	property_areas <- insitu_data |>
		dplyr::distinct(propertyID, property_area_km2) |>
		dplyr::group_by(propertyID) |>
		dplyr::summarize(
			n_areas = length(unique(property_area_km2)),
			property_area_max = max(property_area_km2, na.rm = TRUE),
			# properties with all NA areas get -Inf
			# but the following line changes -Inf to NA
			property_area_km2 = ifelse(
				is.infinite(property_area_max),
				NA,
				property_area_km2
			)
		) |>
		dplyr::ungroup()

	insitu_data |>
		dplyr::left_join(
			property_areas,
			by = join_by(propertyID, property_area_km2)
		) |>
		dplyr::filter(
			!is.na(property_area_km2),
			property_area_km2 >= 1.8,
			effort > 0
		)
}

# need to filter events so that there are at least two events in a timestep
# and there are at least 2 timesteps for each property
dynamic_filter <- function(df) {
	good_events <- df |>
		dplyr::select(propertyID, timestep) |>
		dplyr::group_by(propertyID, timestep) |>
		dplyr::mutate(two_plus_takes = n() >= 2) |>
		dplyr::filter(two_plus_takes) |>
		dplyr::group_by(propertyID) |>
		dplyr::arrange(propertyID, timestep) |>
		dplyr::distinct() |>
		dplyr::mutate(n_timesteps = length(unique(timestep))) |>
		dplyr::filter(n_timesteps >= 2) |>
		dplyr::ungroup() |>
		dplyr::mutate(event_id = paste0(propertyID, "-", timestep)) |>
		dplyr::pull(event_id)

	df |>
		dplyr::mutate(event_id = paste0(propertyID, "-", timestep)) |>
		dplyr::filter(event_id %in% good_events) |>
		dplyr::arrange(propertyID, timestep)
}

take_filter <- function(df) {
	zero_take_prp <- df |>
		dplyr::group_by(propertyID) |>
		dplyr::summarise(sum_take = sum(take)) |>
		dplyr::ungroup() |>
		dplyr::filter(sum_take == 0) |>
		dplyr::pull(propertyID)

	df |> dplyr::filter(!propertyID %in% zero_take_prp)
}

# compute ordering based on time interval midpoints
order_interval <- function(df) {
	df |>
		dplyr::distinct() |>
		dplyr::mutate(midpoint = as.numeric(end.date)) |>
		dplyr::ungroup()
}

# impose stochastic ordering of events by adding jitter
# we are assuming the following order of events when the events have the same midpoint
# e.g., are on the same day:
# 1. (trap or snare), with order random
# 2. (heli or plane), with order random
# 3. hunting
order_stochastic <- function(order.df) {
	n_methods_mid <- order.df |>
		dplyr::select(propertyID, midpoint, method) |>
		dplyr::group_by(propertyID, midpoint) |>
		dplyr::count() |>
		dplyr::arrange(propertyID, midpoint)

	set.seed(8)

	order_df <- order.df
	order_df$jittered_midpoint <- NA
	message("Stochastic ordering...")

	if (!getOption("pbStyle", default = 3)) {
		# will be used if pbStyle is not set by the user (i.e. the default)
		pb <- utils::txtProgressBar(max = nrow(order_df), style = 3)
	} else {
		# if the user calls: options(list(pbStyle = #))
		# the function will use style pbStyle.
		pb <- utils::txtProgressBar(max = nrow(order_df), style = pbStyle)
	}

	for (i in seq_len(nrow(order_df))) {
		if (order_df$method[i] %in% c('Trap', 'Snare')) {
			order_df$jittered_midpoint[i] <- order_df$midpoint[i] +
				runif(1, min = 0, max = .01)
		} else if (order_df$method[i] %in% c('Fixed Wing', 'Helicopter')) {
			order_df$jittered_midpoint[i] <- order_df$midpoint[i] +
				runif(1, min = .02, max = .03)
		} else {
			order_df$jittered_midpoint[i] <- order_df$midpoint[i] +
				runif(1, min = .04, max = .05)
		}
		utils::setTxtProgressBar(pb, i)
	}
	order_df
}

# now compute orders of events based on jittered midpoints
order_of_events <- function(order_df) {
	df <- order_df |>
		dplyr::ungroup() |>
		dplyr::group_by(propertyID, timestep) |>
		dplyr::mutate(
			order = order(jittered_midpoint),
			has_multi = any(order > 1),
			any_ties = any(duplicated(jittered_midpoint)),
			n_survey = n()
		) |>
		dplyr::ungroup() |>
		dplyr::arrange(propertyID, timestep, order) |>
		dplyr::mutate(p = seq_len(n()))

	testthat::expect_all_true(df$has_multi)
	testthat::expect_all_true(!df$any_ties)
	testthat::expect_all_true(df$n_survey > 1)

	df
}

county_codes <- function(df) {
	df |>
		dplyr::rename(statefp = st_gsa_state_cd, countyfp = cnty_gsa_cnty_cd) |>
		dplyr::mutate(
			countyfp = sprintf("%03d", countyfp),
			countyfp = ifelse(cnty_name == "HUMBOLDT (E)", "013", countyfp),
			county_code = as.numeric(paste0(statefp, countyfp)),
			county_code = sprintf("%05d", county_code)
		) |> # need this for joining downstream
		dplyr::select(
			propertyID,
			agrp_prp_id,
			alws_agrprop_id,
			start_dates,
			end_dates,
			st_name,
			cnty_name,
			county_code,
			method,
			trap_count,
			take,
			property_area_km2,
			effort,
			effort_per,
			timestep,
			order,
			n_survey,
			p
		) |>
		dplyr::rename(primary_period = timestep)
}

get_fips <- function(file = "data/fips/national_county.txt") {
	fips <- readr::read_csv(
		file,
		show_col_types = FALSE,
		col_names = c("state", "statefp", "countyfp", "countyname", "classfp"),
		comment = '#'
	)
}

get_ts_length <- function(df) {
	df |>
		dplyr::filter(
			!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS")
		) |>
		dplyr::select(propertyID, primary_period) |>
		dplyr::distinct() |>
		dplyr::group_by(propertyID) |>
		dplyr::filter(
			primary_period == min(primary_period) |
				primary_period == max(primary_period)
		) |>
		dplyr::mutate(delta = c(0, diff(primary_period) + 1)) |>
		dplyr::ungroup() |>
		dplyr::filter(delta != 0)
}

get_n_observations <- function(df, good_props) {
	df |>
		dplyr::filter(propertyID %in% good_props) |>
		dplyr::group_by(propertyID, primary_period) |>
		dplyr::summarise(take = sum(take)) |>
		dplyr::ungroup() |>
		dplyr::group_by(propertyID) |>
		dplyr::count()
}

condition_first_capture <- function(df) {
	good_events <- df |>
		dplyr::group_by(propertyID, timestep) |>
		dplyr::summarise(take = sum(take)) |>
		dplyr::ungroup() |>
		dplyr::group_by(propertyID) |>
		dplyr::mutate(cumulative_take = cumsum(take)) |>
		dplyr::ungroup() |>
		dplyr::filter(cumulative_take > 0) |>
		dplyr::mutate(event_id = paste0(propertyID, "-", timestep)) |>
		dplyr::pull(event_id)

	df |>
		dplyr::mutate(event_id = paste0(propertyID, "-", timestep)) |>
		dplyr::filter(event_id %in% good_events) |>
		dplyr::arrange(propertyID, timestep)
}

create_timestep_df <- function(df) {
	df |>
		dplyr::select(propertyID, primary_period) |>
		unique() |>
		dplyr::arrange(propertyID, primary_period) |>
		dplyr::group_by(propertyID) |>
		dplyr::mutate(observed_timestep = seq_len(n())) |> # timestep is the sequence of primary periods within a property
		dplyr::ungroup() |>
		dplyr::mutate(primary_period = primary_period - min(primary_period) + 1)
}


collate_mcmc_chunks <- function(dest, start = 1) {
	mcmc_dirs <- list.files(dest)
	mcmc_dirs <- setdiff(mcmc_dirs, "modelData.rds")
	mcmc_dirs <- mcmc_dirs[start:length(mcmc_dirs)]
	param_file_name <- "paramSamples.rds"
	state_file_name <- "observedAbundanceSamples.rds"

	# use the first mcmc chunk to initialize storage for each chain
	mcmc_rds <- file.path(dest, mcmc_dirs[1], param_file_name)
	rds <- readr::read_rds(mcmc_rds)
	mcmc <- rds$params
	n_chains <- length(mcmc)
	store_mcmc <- list()
	store_mcmc_state <- list()
	state_count <- 0

	# store each chain, will append below
	for (j in seq_len(n_chains)) {
		store_mcmc[[j]] <- as.matrix(mcmc[[j]])
	}

	state_rds <- file.path(dest, mcmc_dirs[1], state_file_name)
	if (file.exists(state_rds)) {
		state_count <- 1
		mcmc <- readr::read_rds(state_rds)
		for (j in seq_len(n_chains)) {
			store_mcmc_state[[j]] <- as.matrix(mcmc[[j]])
		}
	}

	if (length(mcmc_dirs) >= 2) {
		use_pb <- if_else(length(mcmc_dirs) == 2, FALSE, TRUE)
		if (use_pb) {
			if (!getOption("pbStyle", default = 3)) {
				# will be used if pbStyle is not set by the user (i.e. the default)
				pb <- utils::txtProgressBar(min = 2, max = length(mcmc_dirs), style = 3)
			} else {
				# if the user calls: options(list(pbStyle = #))
				# the function will use style pbStyle.
				pb <- utils::txtProgressBar(
					min = 2,
					max = length(mcmc_dirs),
					style = pbStyle
				)
			}
		}

		# read each mcmc chunk, store each chain from the chunk as a matrix
		for (i in 2:length(mcmc_dirs)) {
			mcmc_rds <- file.path(dest, mcmc_dirs[i], param_file_name)
			rds <- readr::read_rds(mcmc_rds)
			mcmc <- rds$params
			for (j in seq_len(n_chains)) {
				store_mcmc[[j]] <- rbind(store_mcmc[[j]], as.matrix(mcmc[[j]]))
			}

			state_rds <- file.path(dest, mcmc_dirs[i], state_file_name)
			if (file.exists(state_rds)) {
				mcmc <- readr::read_rds(state_rds)

				if (state_count == 0) {
					for (j in seq_len(n_chains)) {
						store_mcmc_state[[j]] <- as.matrix(mcmc[[j]])
					}
				} else {
					for (j in seq_len(n_chains)) {
						store_mcmc_state[[j]] <- rbind(
							store_mcmc_state[[j]],
							as.matrix(mcmc[[j]])
						)
					}
				}
				state_count <- 1
			}
			if (use_pb) utils::setTxtProgressBar(pb, i)
		}
		if (use_pb) close(pb)
	}

	# need to create an mcmc list object to check for convergence
	params <- as.mcmc.list(lapply(store_mcmc, as.mcmc))

	if (state_count == 1) {
		states <- as.mcmc.list(lapply(purrr::compact(store_mcmc_state), as.mcmc))
	} else {
		states <- list()
	}

	return(list(params = params, states = states))
}
