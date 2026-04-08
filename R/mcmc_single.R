#' Run a single MCMC chain for the model
#'
#' @param model_constants A list of constants for the nimble model.
#' @param model_data A list of data for the nimble model.
#' @param model_code The nimble model code.
#' @param init A list of initial values for the nimble model.
#' @param n_iter The number of MCMC iterations to run.
#' @param custom_samplers A data frame specifying custom samplers to use. Should have columns "node" and "type".
#' @param monitors_add A character vector of additional nodes to monitor.
#'
#' @import nimble
#' @return A matrix of MCMC samples.
#' @export

single_mcmc_chain <- function(
  model_constants,
  model_data,
  model_code,
  init,
  n_iter,
  custom_samplers = NULL,
  monitors_add = NULL
) {
  require(nimble)

  Rmodel <- nimbleModel(
    code = model_code, # comes from sourcing the model above
    constants = model_constants,
    data = model_data,
    inits = init,
    calculate = TRUE
  )

  Rmodel$initializeInfo()

  N <- Rmodel$N
  nH_p <- model_constants$nH_p
  n_survey <- model_constants$n_survey
  y_sum <- model_data$y_sum

  for (i in 1:n_survey) {
    N_model <- N[nH_p[i]]
    n <- round(N_model - y_sum[i])
    if (n <= 0) {
      print(i)
      Rmodel$N[nH_p[i]] <- N_model + (abs(n)^2)
    }
  }

  # calc <- Rmodel$calculate()
  # if(is.infinite(calc) | is.nan(calc) | is.na(calc)){
  #   stop(paste0("Model log probability is ", calc))
  # }

  # default MCMC configuration
  mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

  control_rw <- list(
    adaptInterval = 100,
    adaptFactorExponent = 0.6
  )

  mcmcConf$removeSamplers("beta1")
  mcmcConf$addSampler(
    target = "beta1",
    type = "RW_block",
    control = control_rw
  )

  mcmcConf$removeSamplers("log_nu")
  mcmcConf$addSampler(
    target = "log_nu",
    type = "slice"
  )

  # if specified, change nodes to specified parameters
  if (!is.null(custom_samplers)) {
    for (i in seq_len(nrow(custom_samplers))) {
      node <- custom_samplers$node[i]
      type <- custom_samplers$type[i]
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, type)
    }
  }

  if (!is.null(monitors_add)) {
    mcmcConf$addMonitors(monitors_add)
  }

  Rmcmc <- buildMCMC(mcmcConf)
  Cmodel <- compileNimble(Rmodel)

  # need to define this in the global environment to continue mcmc when running in parallel
  Cmcmc <<- compileNimble(Rmcmc)

  Cmcmc$run(niter = n_iter, nburnin = n_iter / 2, thin = 10)
  samples <- as.matrix(Cmcmc$mvSamples)
  return(samples)
}
