message("\n==== check_mcmc.R ====")

if (!exists("params_mcmc")) {
  stop("check_mcmc.R requires params_mcmc in the source environment.")
}

if (exists("mcmc_dir")) {
  message("MCMC directory: ", mcmc_dir)
}

if (exists("stuck_params") && nrow(stuck_params) > 0) {
  message("Stuck monitored parameters:")
  print(stuck_params)
}

chain_ids <- seq_along(params_mcmc)
print_first_mcmc_iteration(params_mcmc, params_check, chain_ids = chain_ids)

message("\n==== Unique values by monitored node ====")
for (i in seq_along(params_mcmc)) {
  chain_samples <- as.matrix(params_mcmc[[i]])
  message("Chain ", chain_ids[i])

  for (node in params_check) {
    matched_nodes <- grep(node, colnames(chain_samples), value = TRUE, fixed = TRUE)

    if (length(matched_nodes) == 0) {
      message("  ", node, ": <not found>")
      next
    }

    unique_counts <- vapply(
      matched_nodes,
      function(node_name) {
        node_values <- chain_samples[, node_name]
        length(unique(node_values[!is.na(node_values)]))
      },
      integer(1)
    )

    values <- paste0(matched_nodes, " unique=", unique_counts)
    message("  ", node, ": ", paste(values, collapse = ", "))
  }
}
