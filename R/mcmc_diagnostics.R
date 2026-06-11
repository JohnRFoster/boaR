#' Calculate MCMC diagnostics for a nimble model
#' Currently runs by collating chunks of MCMC
#' samples and then calculating diagnostics on the collated samples.
#'
#'@description Calculate MCMC diagnostics for a nimble model

mcmc_diagnostics <- function(
	mcmc_dir,
	dest,
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
		dplyr::slice(draws) |>
		dplyr::mutate(np = config$np)

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
		dplyr::slice(draws) |>
		dplyr::mutate(np = config$np)

	readr::write_rds(state_samples, file.path(dest, "stateSamples.rds"))

	if (make_traceplot) {
		message("\n==== Making traceplots ====")
		traceplot(
			params_mcmc_list = params_mcmc_list,
			nodes_2_plot = params_check,
			n_chains = n_chains,
			dest = dest
		)
	}
}
