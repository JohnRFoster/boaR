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
		as_tibble()
	if (ncol(s) == 1) {
		colnames(s) <- node
	}
	s
}

# get parameter nodes using subset_mcmc
subset_params <- function() {
	require(dplyr)
	require(purrr)
	map_dfc(lapply(params_check, subset_mcmc), as_tibble) |> as.matrix()
}

# get observed abundance nodes using subset_mcmc
subset_N_observed <- function() {
	nodes <- paste0("N[", model_constants$N_full_unique, "]")
	map_dfc(lapply(nodes, subset_mcmc), as_tibble) |> as.matrix()
}

# get unobserved abundance nodes using subset_mcmc
subset_N_unobserved <- function() {
	nodes <- paste0("N[", model_constants$N_quant_unique, "]")
	N_unobserved <- map_dfc(lapply(nodes, subset_mcmc), as_tibble) |>
		as.matrix()
	t(apply(N_unobserved, 2, quantile, c(0.025, 0.5, 0.975)))
}

write_out_p <- function(p, d, dest) {
	f <- "paramSamples.rds"
	write_rds(list(params = p, diagnostic = d), file.path(dest, f))
}

write_abundance <- function(no, nu, dest) {
	f <- "observedAbundanceSamples.rds"
	write_rds(no, file.path(dest, f))

	f <- "unobservedAbundanceQuantiles.rds"
	write_rds(nu, file.path(dest, f))
}

# miscellaneous helper functions

# need all timesteps whether there are observations or not
create_all_primary_periods <- function(df) {
	pp_min_max <- df |>
		select(property, primary_period) |>
		distinct() |>
		group_by(property) |>
		dplyr::filter(
			primary_period == min(primary_period) |
				primary_period == max(primary_period)
		) |>
		ungroup()

	properties <- unique(pp_min_max$property)

	all_pp <- tibble()
	message("\nInclude all primary periods")
	pb <- txtProgressBar(max = length(properties), style = 1)
	for (i in seq_along(properties)) {
		pid <- pp_min_max |> dplyr::filter(property == properties[i])
		p_min <- min(pid$primary_period)
		p_max <- max(pid$primary_period)
		pp <- tibble(
			property = properties[i],
			primary_period = p_min:p_max,
			timestep = 1:length(p_min:p_max)
		)
		all_pp <- bind_rows(all_pp, pp)
		setTxtProgressBar(pb, i)
	}
	close(pb)
	all_pp |> mutate(n_id = 1:n())
}

# need to know the total number of timesteps in each property (sampled or not) for indexing
n_timesteps <- function(df) {
	df |>
		group_by(property) |>
		dplyr::filter(timestep == max(timestep)) |>
		pull(timestep)
}

# index (as a matrix) for tracking abundance in long format (converting wide to long)
# used in the process model
N_lookup_table <- function(df) {
	df |>
		mutate(n_id = 1:n()) |>
		select(-primary_period) |>
		pivot_wider(names_from = timestep, values_from = n_id) |>
		select(-property) |>
		as.matrix()
}

# calculate the cumulative number of pigs taken as a primary period progresses
removed_in_pp_cumsum <- function(df) {
	df |>
		group_by(property, primary_period) |>
		mutate(ysum = cumsum(take) - take) |>
		ungroup() |>
		pull(ysum)
}

# the total number of pigs taken in a primary period across all methods
# including periods without removals effort (equal to 0)
# wide format
total_take <- function(df_take, df_pp) {
	sum_take <- df_take |>
		group_by(property, primary_period) |>
		summarise(sum_take = sum(take)) |>
		ungroup()

	left_join(df_pp, sum_take, by = join_by(property, primary_period)) |>
		mutate(sum_take = if_else(is.na(sum_take), 0, sum_take)) |>
		select(-primary_period, -n_id) |>
		pivot_wider(names_from = timestep, values_from = sum_take) |>
		select(-property) |>
		as.matrix()
}

# index (long format) for which primary periods have removal effort
# and therefore included in the data model
N_lookup_data <- function(df_take, df_pp) {
	tH <- df_take |>
		select(property, primary_period)

	df_pp |>
		select(property, primary_period, n_id) |>
		right_join(tH, by = join_by(property, primary_period)) |>
		pull(n_id)
}

# need start and end indicies for data model
# used for estimiating p
create_start_end <- function(df_take, df_pp) {
	start <- end <- numeric(nrow(df_take))

	df <- left_join(df_take, df_pp, by = join_by(primary_period, property))

	message("Creating start/end indicies")
	pb <- txtProgressBar(max = nrow(df), style = 1)
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
			assertthat::are_equal(idx, start[i]:end[i])
		}
		setTxtProgressBar(pb, i)
	}
	close(pb)
	tibble(start = start, end = end)
}

# informed hyper parameters for beta distribution on global pig survival
create_surv_prior <- function(interval, data_repo) {
	data <- read_csv(
		file.path(data_repo, "insitu/Vital_Rate_Data.csv"),
		show_col_types = FALSE
	)

	data_usa <- data |>
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
		select(
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
		mutate(
			weeks = as.numeric(time.period.end - time.period.start) / 7,
			weeks4 = weeks / interval,
			survival.per.4week = survival.prop^(1 / weeks4),
			logit.survival.per.4week = boot::logit(survival.per.4week)
		) |>
		dplyr::filter(survival.per.4week > 0) |>
		mutate(scale_factor = survival.per.4week / survival.prop)

	surv_mu_summary <- surv_mu |>
		dplyr::summarise(
			mu = mean(survival.per.4week),
			mu.logit = mean(logit.survival.per.4week)
		)

	surv_var <- surv_data |>
		dplyr::filter(survival.var.type %in% c("SD", "95% CI"))

	surv_sd <- surv_var |>
		dplyr::filter(survival.var.type == "SD") |>
		mutate(sd = as.numeric(survival.var))

	surv_sd_calc <- surv_var |>
		dplyr::filter(survival.var.type == "95% CI") |>
		mutate(
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
		group_by(unique.ID) |>
		summarise(sd = max(sd_high, sd_low))

	surv_var_join <- left_join(surv_var, surv_sd_calc, by = join_by(unique.ID)) |>
		dplyr::filter(survival.var.type != "SD")

	scale_ids <- surv_mu |>
		select(unique.ID, scale_factor)

	surv_variance <- bind_rows(surv_var_join, surv_sd) |>
		left_join(scale_ids, by = join_by(unique.ID)) |>
		mutate(
			variance = sd^2,
			variance.4week = variance * scale_factor^2,
			sd.4week = sqrt(variance.4week)
		)

	surv_sd_summary <- surv_variance |>
		pull(sd.4week) |>
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
		select(all_of(cols)) |>
		as.matrix()
}

get_prior_hyperparams <- function(post_round, posterior_path = NULL) {
	if (post_round == "first") {
		read_rds("data/originalFirstRoundHyperparameters.rds")
	} else if (post_round == "last") {
		read_rds("data/originalLastRoundHyperparameters.rds")
	} else if (post_round == "create_new") {
		post <- read_rds(posterior_path)

		get_vec <- function(df, node) {
			df |>
				select(contains(node)) |>
				pivot_longer(cols = everything(), names_to = "node") |>
				group_by(node) |>
				summarise(mu = mean(value), tau = 1 / var(value))
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
		group_by(property, property_area_km2) |>
		reframe(total_take = sum(take)) |>
		mutate(
			min_n = total_take + 5, # add a small buffer to allow for some unobserved pigs
			abundance_min = min_n,
			abundance_max = round(property_area_km2 * 100) + min_n
		)
}

create_ids <- function(df) {
	df |>
		mutate(
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
	psrf <- gelman.diag(mcmc, multivariate = FALSE)$psrf
	effective_samples <- effectiveSize(mcmc)

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
	if (interval == 28 & !create_new) {
		primary_periods <- read_rds("data/primaryPeriods28Days.rds")

		timesteps <- map(
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

		primary_periods <- tibble(start_dates, end_dates) |>
			mutate(timestep = 1:n())
		primary_periods$month <- month(timestep_df$end_dates)
		primary_periods$year <- year(timestep_df$end_dates)
	} else {
		stop("Interval not supported")
	}

	timesteps <- map(
		seq_len(nrow(df)),
		~ find_timestep(primary_periods, df$start.date[.x], df$end.date[.x]),
		.progress = TRUE
	) |>
		unlist() |>
		as.integer()

	tmp <- df |>
		mutate(timestep = timesteps)

	df_time <- left_join(tmp, primary_periods)
	df_test <- df_time |> dplyr::filter(!is.na(timestep))

	testthat::test_that("Start and end dates align with known primary periods", {
		expect_true(all(df_test$start.date >= df_test$start_dates))
		expect_true(all(df_test$end.date <= df_test$end_dates))
	})

	df_test |>
		arrange(propertyID, timestep)
}

# resolve duplicate property values - when there are multiple values, take max
resolve_duplicate <- function(insitu_data) {
	property_areas <- insitu_data |>
		distinct(propertyID, property_area_km2) |>
		group_by(propertyID) |>
		summarize(
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
		ungroup()

	insitu_data |>
		left_join(property_areas, by = join_by(propertyID, property_area_km2)) |>
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
		select(propertyID, timestep) |>
		group_by(propertyID, timestep) |>
		mutate(two_plus_takes = n() >= 2) |>
		dplyr::filter(two_plus_takes) |>
		group_by(propertyID) |>
		arrange(propertyID, timestep) |>
		distinct() |>
		mutate(n_timesteps = length(unique(timestep))) |>
		dplyr::filter(n_timesteps >= 2) |>
		ungroup() |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		pull(event_id)

	df |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		dplyr::filter(event_id %in% good_events) |>
		arrange(propertyID, timestep)
}

take_filter <- function(df) {
	zero_take_prp <- df |>
		group_by(propertyID) |>
		summarise(sum_take = sum(take)) |>
		ungroup() |>
		dplyr::filter(sum_take == 0) |>
		pull(propertyID)

	df |> dplyr::filter(!propertyID %in% zero_take_prp)
}

# compute ordering based on time interval midpoints
order_interval <- function(df) {
	df |>
		distinct() |>
		mutate(midpoint = as.numeric(end.date)) |>
		ungroup()
}

# impose stochastic ordering of events by adding jitter
# we are assuming the following order of events when the events have the same midpoint
# e.g., are on the same day:
# 1. (trap or snare), with order random
# 2. (heli or plane), with order random
# 3. hunting
order_stochastic <- function(order.df) {
	n_methods_mid <- order.df |>
		select(propertyID, midpoint, method) |>
		group_by(propertyID, midpoint) |>
		count() |>
		arrange(propertyID, midpoint)

	set.seed(8)

	order_df <- order.df
	order_df$jittered_midpoint <- NA
	message("Stochastic ordering...")
	pb <- txtProgressBar(max = nrow(order_df), style = 1)
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
		setTxtProgressBar(pb, i)
	}
	order_df
}

# now compute orders of events based on jittered midpoints
order_of_events <- function(order_df) {
	df <- order_df |>
		ungroup() |>
		group_by(propertyID, timestep) |>
		mutate(
			order = order(jittered_midpoint),
			has_multi = any(order > 1),
			any_ties = any(duplicated(jittered_midpoint)),
			n_survey = n()
		) |>
		ungroup() |>
		arrange(propertyID, timestep, order) |>
		mutate(p = seq_len(n()))

	testthat::expect_all_true(df$has_multi)
	testthat::expect_all_true(!df$any_ties)
	testthat::expect_all_true(df$n_survey > 1)

	df
}

county_codes <- function(df) {
	df |>
		rename(statefp = st_gsa_state_cd, countyfp = cnty_gsa_cnty_cd) |>
		mutate(
			countyfp = sprintf("%03d", countyfp),
			countyfp = ifelse(cnty_name == "HUMBOLDT (E)", "013", countyfp),
			county_code = as.numeric(paste0(statefp, countyfp)),
			county_code = sprintf("%05d", county_code)
		) |> # need this for joining downstream
		select(
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
		rename(primary_period = timestep)
}

get_fips <- function(file = "data/fips/national_county.txt") {
	fips <- read_csv(
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
		select(propertyID, primary_period) |>
		distinct() |>
		group_by(propertyID) |>
		dplyr::filter(
			primary_period == min(primary_period) |
				primary_period == max(primary_period)
		) |>
		mutate(delta = c(0, diff(primary_period) + 1)) |>
		ungroup() |>
		dplyr::filter(delta != 0)
}

get_n_observations <- function(df, good_props) {
	df |>
		dplyr::filter(propertyID %in% good_props) |>
		group_by(propertyID, primary_period) |>
		summarise(take = sum(take)) |>
		ungroup() |>
		group_by(propertyID) |>
		count()
}

condition_first_capture <- function(df) {
	good_events <- df |>
		group_by(propertyID, timestep) |>
		summarise(take = sum(take)) |>
		ungroup() |>
		group_by(propertyID) |>
		mutate(cumulative_take = cumsum(take)) |>
		ungroup() |>
		dplyr::filter(cumulative_take > 0) |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		pull(event_id)

	df |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		dplyr::filter(event_id %in% good_events) |>
		arrange(propertyID, timestep)
}

create_timestep_df <- function(df) {
	df |>
		select(propertyID, primary_period) |>
		unique() |>
		arrange(propertyID, primary_period) |>
		group_by(propertyID) |>
		mutate(observed_timestep = 1:n()) |> # timestep is the sequence of primary periods within a property
		ungroup() |>
		mutate(primary_period = primary_period - min(primary_period) + 1)
}
