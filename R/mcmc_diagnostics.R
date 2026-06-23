#' Calculate MCMC diagnostics for a nimble model
#' Currently runs by collating chunks of MCMC
#' samples and then calculating diagnostics on the collated samples.
#'
#'@description Calculate MCMC diagnostics for a nimble model
#' @param mcmc_dir directory where MCMC chunks are stored
#' @param dest directory where diagnostics and posterior samples will be saved
#' @param data input data frame used for nimble fit
#' @param params_check vector of parameters (nodes) to assess convergence
#' @param n_mcmc number random posterior samples to save after burnin
#' @param effective_size minimum effective sample size for each parameter
#' @param max_psrf maximum potential scale reduction factor (PSRF) to assess if mcmc is behaving well.
#' @param verbose print diagnostic messages to the console
#' @param make_traceplot make traceplots for the parameters in params_check
#' @export

mcmc_diagnostics <- function(
	mcmc_dir,
	dest,
	data,
	params_check,
	n_mcmc = 5000,
	effective_size = 1000,
	max_psrf = 15,
	verbose = TRUE,
	make_traceplot = TRUE
) {
	# get samples
	mcmc_list <- collate_mcmc_chunks(mcmc_dir, start = 1)
	params_mcmc_list <- mcmc_list$params

	message("\n==== MCMC Information ====")
	total_iter <- nrow(params_mcmc_list[[1]])
	message("Total iterations: ", total_iter)

	n_chains <- length(params_mcmc_list)
	message("Number of chains: ", n_chains)

	# calculate psrf (convergence stat) and effective sample size
	message("\n==== Pre-burnin diagnostics ====")
	diagnostic <- continue_mcmc(
		params_mcmc_list,
		effective_size = effective_size,
		max_psrf = max_psrf,
		verbose = verbose
	)

	GBR <- coda::gelman.plot(params_mcmc_list)
	burnin <- GBR$last.iter[
		tail(which(apply(GBR$shrink[,, 2] > 1.1, 1, any)), 1) + 1
	]
	message("Burnin: ", burnin)
	if (is.na(burnin)) {
		burnin <- round(total_iter / 2)
		message(
			"Burnin has not occured! Setting to half of total iterations (",
			burnin,
			")"
		)
	}

	if (burnin >= total_iter * 0.9) {
		message(
			"Warning: Burnin is at the end of the MCMC chains! Consider rerunning for more samples."
		)
		burnin <- round(total_iter / 2)
	}

	params_burnin <- window(params_mcmc_list, start = burnin)

	# calculate psrf (convergence stat) and effective sample size
	message("\n==== Post-burnin diagnostics ====")
	diagnostic <- continue_mcmc(
		params_burnin,
		effective_size = effective_size,
		max_psrf = max_psrf,
		verbose = verbose
	)

	# get a random sample of the posterior and save
	posterior_burnin <- params_burnin |>
		as.matrix()

	draws <- sample.int(nrow(posterior_burnin), n_mcmc, replace = TRUE)
	posterior_samples <- posterior_burnin |>
		tibble::as_tibble() |>
		dplyr::slice(draws)

	readr::write_rds(
		posterior_samples,
		file.path(dest, "posteriorSamples.rds")
	)

	states_mcmc_list <- mcmc_list$states

	states_burnin <- window(states_mcmc_list, start = burnin)

	posterior_burnin2 <- states_burnin |>
		as.matrix()

	state_samples <- posterior_burnin2 |>
		tibble::as_tibble() |>
		dplyr::slice(draws)

	readr::write_rds(state_samples, file.path(dest, "stateSamples.rds"))

	if (make_traceplot) {
		message("\n==== Making traceplots ====")
		trace_plot(
			params_mcmc_list = params_mcmc_list,
			nodes_2_plot = params_check,
			n_chains = n_chains,
			dest = dest
		)
	}

	density <- density_stats(state_samples, data)
	write_rds(density, file.path(dest, "densitySummaries.rds"))
}
