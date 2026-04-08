#' Create A List For Nimble
#'
#' @param df A data frame containing the survey data.
#' @param interval The length of the primary period in days.
#' @param data_repo The path to the data repository.
#' @param post_round "first" or "last"; determines which set of prior hyperparameters to use.
#'
#' @return A list of constants for the nimble model.
#' @export

nimble_constants <- function(
  df,
  interval,
  data_repo,
  post_round
) {
  all_primary_periods <- create_all_primary_periods(df)
  n_time_prop <- n_timesteps(all_primary_periods)
  n_method <- length(unique(df$method))
  nH <- N_lookup_table(all_primary_periods)
  nH_p <- N_lookup_data(df, all_primary_periods)
  N_full_unique <- nH_p |> unique()
  N_quant_unique <- setdiff(seq(1, max(N_full_unique)), N_full_unique)
  rem <- total_take(df, all_primary_periods)
  targets::tar_assert_true(all(rem[, 1] != 0))
  X <- create_X(df)
  start_end <- create_start_end(df, all_primary_periods)
  n1_priors <- get_n1_prior(df)

  constants <- list(
    n_survey = nrow(df),
    n_ls = length(data_ls), # mean litter size year from VerCauteren et al. 2019 pg 63
    n_property = length(unique(df$property)),
    n_first_survey = length(which(df$order == 1)),
    n_not_first_survey = length(which(df$order != 1)),
    n_method = n_method,
    n_betaP = n_method * ncol(X),
    beta_p_row = rep(seq_len(n_method), each = ncol(X)),
    beta_p_col = rep(seq_len(ncol(X)), n_method),
    n_time_prop = n_time_prop,
    nH = nH,
    nH_p = nH_p,
    N_full_unique = N_full_unique,
    N_quant_unique = N_quant_unique,
    rem = rem,
    log_pi = log(pi),
    first_survey = which(df$order == 1),
    not_first_survey = which(df$order != 1),
    m_p = ncol(X),
    start = start_end$start,
    end = start_end$end,
    method = as.numeric(as.factor(df$method)),
    pp_len = interval,
    n1_min = n1_priors$abundance_min,
    n1_max = n1_priors$abundance_max
  )

  prior_hyperparams <- get_prior_hyperparams(post_round)

  append(constants, prior_hyperparams)
}
