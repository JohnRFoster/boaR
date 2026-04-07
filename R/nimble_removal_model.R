#' The Nimble code for the removal model
# TODO update this when published
#' @description This model is based on the removal model described in Foster et al. 2025 (doi.org/10.1101/2025.10.16.680799)
#' @return model code
#' @export

nimble_removal_model <- function() {
  nimble::nimbleCode({
    # priors
    for (i in 1:n_method) {
      log_rho[i] ~ dnorm(log_rho_mu[i], tau = log_rho_tau[i])
    }

    for (i in 1:2) {
      p_mu[i] ~ dnorm(p_mu_mu[i], tau = p_mu_tau[i])
      logit(p_unique[i]) <- p_mu[i]

      log_gamma[i] ~ dnorm(log_gamma_mu[i], tau = log_gamma_tau[i])
    }

    # non time varying coefficients - observation model
    for (i in 1:n_method) {
      beta1[i] ~ dnorm(beta1_mu[i], tau = beta1_tau[i])
    }

    for (i in 1:n_betaP) {
      beta_p[beta_p_row[i], beta_p_col[i]] ~ dnorm(
        beta_p_mu[i],
        tau = beta_p_tau[i]
      )
    }

    # estimate apparent survival
    phi_mu ~ dbeta(phi_mu_a, phi_mu_b)
    psi_phi ~ dgamma(psi_shape, psi_rate)
    a_phi <- phi_mu * psi_phi
    b_phi <- (1 - phi_mu) * psi_phi

    log_nu ~ dnorm(log_nu_mu, tau = log_nu_tau) # mean litter size
    log(nu) <- log_nu

    ## convert to expected number of pigs per primary period
    log_zeta <- log(pp_len) + log_nu - log(365)
    log(zeta) <- log_zeta
    for (i in 1:n_ls) {
      J[i] ~ dpois(nu)
    }

    for (i in 1:n_survey) {
      log_potential_area[i] <- calc_log_potential_area(
        log_rho = log_rho[1:n_method],
        log_gamma = log_gamma[1:2],
        p_unique = p_unique[1:2],
        log_effort_per = log_effort_per[i],
        effort_per = effort_per[i],
        n_trap_m1 = n_trap_m1[i],
        log_pi = log_pi,
        method = method[i]
      )

      # probability of capture, given that an individual is in the surveyed area
      log_theta[i] <- log(
        ilogit(
          beta1[method[i]] + inprod(X_p[i, 1:m_p], beta_p[method[i], 1:m_p])
        )
      ) +
        min(0, log_potential_area[i] - log_survey_area_km2[i])

      # likelihood
      y[i] ~ dpois(p[i] * (N[nH_p[i]] - y_sum[i]))
    }

    # the probability an individual is captured on the first survey
    for (i in 1:n_first_survey) {
      log(p[first_survey[i]]) <- log_theta[first_survey[i]]
    }

    # the probability an individual is captured after the first survey
    for (i in 1:n_not_first_survey) {
      log(p[not_first_survey[i]]) <- log_theta[start[not_first_survey[i]]] +
        sum(log(
          1 -
            exp(log_theta[start[not_first_survey[i]]:end[not_first_survey[i]]])
        ))
    }

    for (i in 1:n_property) {
      lambda_1[i] ~ dunif(n1_min[i], n1_max[i])
      N[nH[i, 1]] ~ dpois(round(lambda_1[i]))

      # population growth across time steps
      for (j in 2:n_time_prop[i]) {
        # loop through every PP, including missing ones

        lambda[nH[i, j - 1]] <- (N[nH[i, j - 1]] - rem[i, j - 1]) *
          zeta /
          2 +
          (N[nH[i, j - 1]] - rem[i, j - 1]) * phi[nH[i, j - 1]]

        N[nH[i, j]] ~ dpois(lambda[nH[i, j - 1]])
        phi[nH[i, j - 1]] ~ dbeta(a_phi, b_phi)
      }
    }
  })
}


#' Calculate the log of the potential area surveyed
#' @param log_rho vector of log rho coefficients (scaling effort to area) for each method
#' @param log_gamma vector of log gamma coefficients (saturating effect) for traps and snares
#' @param p_unique vector of probabilities that a trap captures a unique individual for traps and snares
#' @param log_effort_per scalar of log effort per survey
#' @param effort_per scalar of effort per survey
#' @param n_trap_m1 scalar of number of traps minus 1 for each survey
#' @param log_pi log of pi
#' @param method method index for each survey
#' @return log of the potential area surveyed
#' @export
calc_log_potential_area <- nimble::nimbleFunction(
  run = function(
    log_rho = double(1),
    log_gamma = double(1),
    p_unique = double(1),
    log_effort_per = double(0),
    effort_per = double(0),
    n_trap_m1 = double(0),
    log_pi = double(0),
    method = double(0)
  ) {
    m <- method

    if (m <= 3) {
      # firearms, fixed wing, and helicopter
      log_potential_area <- log_rho[m] + log_effort_per
    } else {
      # traps and snares
      log_potential_area <- log_pi +
        (2 *
          (log_rho[m] +
            log_effort_per -
            log(exp(log_gamma[m - 3]) + effort_per))) +
        log(1 + (p_unique[m - 3] * n_trap_m1))
    }
    return(log_potential_area)
    returnType(double(0))
  },
  buildDerivs = TRUE
)
