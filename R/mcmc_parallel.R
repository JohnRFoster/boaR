#' Function for running nimble in parallel
#'
#'@description Run nimble on specified cluster, check for convergence and save samples periodically
#'@param cl Cluster made using the parallel package
#'@param model_code the nimble model code
#'@param model_constants model constants, as a list
#'@param model_data model data, as a list
#'@param model_inits initial model values, as a list, one list element for each chain
#'@param params_check vector of parameters (nodes) to assess convergence
#'@param n_iters the number of iterations to run in each chunk after the model has compiled
#'@param dest output filepath. within this directory each 'chunk' will be saved into it's own directory
#'@param monitors_add vector of nodes to monitor in addition to default nodes
#'@param custom_samplers data frame of mcmc config. Nodes in one column, sampler types in another
#'
#'
#'@export

mcmc_parallel <- function(
	cl,
	model_code,
	model_constants,
	model_data,
	model_inits,
	params_check,
	n_iters,
	dest,
	monitors_add = NULL,
	custom_samplers = NULL
) {
	export <- c(
		"single_mcmc_chain",
		"continue_sampling",
		"subset_mcmc",
		"subset_params",
		"subset_N_observed",
		"subset_N_unobserved",
		"model_code",
		"model_data",
		"model_constants",
		"nimble_inits",
		"n_iters",
		"custom_samplers",
		"monitors_add",
		"params_check",
		"calc_log_potential_area"
	)

	parallel::clusterExport(cl, export, envir = environment())

	# for (i in seq_along(cl)) {
	# 	init <- model_inits[[i]]
	# 	parallel::clusterExport(cl[i], "init", envir = environment())
	# }

	# initialize model and first samples
	c <- 1
	start <- Sys.time()
	out <- parallel::clusterEvalQ(
		cl,
		single_mcmc_chain(
			model_constants = model_constants,
			model_data = model_data,
			model_code = model_code,
			init = nimble_inits(model_constants, model_data, buffer = 200),
			n_iter = n_iters,
			custom_samplers = NULL,
			monitors_add = monitors_add
		)
	)
	message("Model compile and initial 1000 iterations completed in:")
	print(round(Sys.time() - start, 2))

	start2 <- Sys.time()
	out2 <- parallel::clusterEvalQ(cl, continue_sampling())
	message("Additional ", n_iters, " iterations completed in:")
	print(round(Sys.time() - start2, 2))

	message("\n", n_iters * c, " total iterations completed in:")
	print(round(Sys.time() - start, 2))

	# use mcmc on clusters to subset parameters, observed states, and unobserved states
	params <- parallel::clusterEvalQ(cl, subset_params())
	params <- coda::as.mcmc.list(lapply(params, as.mcmc))

	N_observed <- parallel::clusterEvalQ(cl, subset_N_observed())
	N_observed <- coda::as.mcmc.list(lapply(N_observed, as.mcmc))

	N_unobserved <- parallel::clusterEvalQ(cl, subset_N_unobserved()) |>
		as.matrix()

	params_mcmc_list <- coda::as.mcmc.list(lapply(params, as.mcmc))
	efsize <- 1000
	diagnostic <- continue_mcmc(
		mcmc = params_mcmc_list,
		effective_size = efsize,
		max_psrf = 15,
		verbose = TRUE
	)

	c_dir <- sprintf("%04d", c)
	path <- file.path(dest, c_dir)
	if (!dir.exists(path)) {
		dir.create(path, showWarnings = FALSE, recursive = TRUE)
	}

	write_out_p(params, diagnostic, path)

	converged <- all(diagnostic$psrf[, 2] < 1.1)
	if (converged) {
		write_abundance(N_observed, N_unobserved, path)
	}

	continue <- !diagnostic$done
	write_N <- FALSE
	while (continue) {
		c <- c + 1
		resetMV <- TRUE
		parallel::clusterExport(cl, "resetMV", envir = environment())

		start2 <- Sys.time()
		out2 <- parallel::clusterEvalQ(cl, continue_sampling())
		message("Additional ", n_iters, " iterations completed in:")
		print(round(Sys.time() - start2, 2))

		total_iters <- n_iters * c
		message("\n", total_iters, " total iterations completed in:")
		print(round(Sys.time() - start, 2))

		# use mcmc on clusters to subset parameters, observed states, and unobserved states
		params <- parallel::clusterEvalQ(cl, subset_params())
		params <- coda::as.mcmc.list(lapply(params, as.mcmc))

		c_dir <- sprintf("%04d", c)
		path <- file.path(dest, c_dir)
		if (!dir.exists(path)) {
			dir.create(path, showWarnings = FALSE, recursive = TRUE)
		}

		diagnostic <- NULL
		write_out_p(params, diagnostic, path)

		N_observed <- parallel::clusterEvalQ(cl, subset_N_observed())
		N_observed <- coda::as.mcmc.list(lapply(N_observed, as.mcmc))
		N_unobserved <- parallel::clusterEvalQ(cl, subset_N_unobserved()) |>
			as.matrix()

		params_mcmc_list <- collate_mcmc_chunks(dest)$params
		diagnostic <- continue_mcmc(
			mcmc = params_mcmc_list,
			effective_size = efsize,
			max_psrf = 15,
			verbose = TRUE
		)
		converged <- all(diagnostic$psrf[, 2] < 1.1)
		write_N <- dplyr::if_else(converged, TRUE, FALSE)

		if (converged || write_N) {
			write_abundance(N_observed, N_unobserved, path)
		}

		continue <- !diagnostic$done
		message("=================================================")

		if (c == 100) continue <- FALSE
	}

	if_else(diagnostic$done, TRUE, FALSE)
}
