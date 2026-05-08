#' The Nimble code for the removal model
# TODO update this when published
#' @description This model is based on the removal model described in Foster et al. 2025 (doi.org/10.1101/2025.10.16.680799)
#' @return model code
#' @import nimble
#' @export

nimble_removal_model <- function() {
  nimble::nimbleCode({
    # priors
    # no i index when using a single method
    # single method used in the data, scenarios are:
    # shooting only
    # traps only
    # snares only
    if (single_method) {
      if (use_traps_or_snares) {
        p_mu ~ dnorm(p_mu_mu, tau = p_mu_tau)
        logit(p_unique) <- p_mu

        log_gamma ~ dnorm(log_gamma_mu, tau = log_gamma_tau)
      }

      # log_rho is needed for every method
      log_rho ~ dnorm(log_rho_mu, tau = log_rho_tau)

      # observation model is used for every method
      beta1 ~ dnorm(beta1_mu, tau = beta1_tau)

      for (i in 1:n_betaP) {
        beta_p[1, i] ~ dnorm(beta_p_mu[i], tau = beta_p_tau[i])
      }
    } else {
      # multiple methods used in the data, scenarios are:
      # shooting method(s) with traps or snares
      # shooting method(s) with traps and snares
      # shooting only
      # traps and snares only

      # need both use_traps_or_snares and use_traps_and_snares
      # so if traps aren't used p and gamma will be skipped
      if (use_traps_or_snares) {
        p_mu ~ dnorm(p_mu_mu, tau = p_mu_tau)
        logit(p_unique) <- p_mu

        log_gamma ~ dnorm(log_gamma_mu, tau = log_gamma_tau)
      }

      if (use_traps_and_snares) {
        for (i in 1:2) {
          p_mu[i] ~ dnorm(p_mu_mu[i], tau = p_mu_tau[i])
          logit(p_unique[i]) <- p_mu[i]

          log_gamma[i] ~ dnorm(log_gamma_mu[i], tau = log_gamma_tau[i])
        }
      }

      # regardless of the combination of methods used
      # log_rho is needed for each method
      for (i in 1:n_method) {
        log_rho[i] ~ dnorm(log_rho_mu[i], tau = log_rho_tau[i])
      }

      # observation model is used for every method

      for (i in 1:n_method) {
        beta1[i] ~ dnorm(beta1_mu[i], tau = beta1_tau[i])
      }

      # row indexes methods
      # column indexes land cover covariates
      # if all methods were used this would be [5, 3]
      # row 1 = FIREARMS, [c_road_den, c_rugged, c_canopy]
      # row 2 = FIXED WING, [c_road_den, c_rugged, c_canopy]
      # row 3 = HELICOPTER, [c_road_den, c_rugged, c_canopy]
      # row 4 = SNARE, [c_road_den, c_rugged, c_canopy]
      # row 5 = TRAPS, [c_road_den, c_rugged, c_canopy]
      # this loop builds by row
      for (i in 1:n_betaP) {
        beta_p[beta_p_row[i], beta_p_col[i]] ~ dnorm(
          beta_p_mu[i],
          tau = beta_p_tau[i]
        )
      }
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
      if (single_method) {
        if (use_shooting) {
          log_potential_area[i] <- log_rho + log_effort_per[i]
        } else {
          log_potential_area[i] <- log_pi +
            (2 *
              (log_rho +
                log_effort_per[i] -
                log(exp(log_gamma) + effort_per[i]))) +
            log(1 + (p_unique * n_trap_m1[i]))
        }
      } else {
        # multiple methods used in the data, scenarios are:
        # shooting method(s) with traps or snares
        # shooting method(s) with traps and snares
        # shooting only
        # traps and snares only

        if (use_traps_or_snares) {
          # either traps or snares so 1 gamma and 1 p
          log_potential_area[i] <- calc_log_potential_area(
            log_rho = log_rho[1:n_method],
            log_gamma = log_gamma,
            p_unique = p_unique,
            log_effort_per = log_effort_per[i],
            effort_per = effort_per[i],
            n_trap_m1 = n_trap_m1[i],
            log_pi = log_pi,
            shooting = shooting[i],
            method = method[i],
            ts_id = ts_id[i]
          )
        }
        if (use_traps_and_snares) {
          # both traps and snares so 2 gamma and 2 p
          log_potential_area[i] <- calc_log_potential_area(
            log_rho = log_rho[1:n_method],
            log_gamma = log_gamma[1:2],
            p_unique = p_unique[1:2],
            log_effort_per = log_effort_per[i],
            effort_per = effort_per[i],
            n_trap_m1 = n_trap_m1[i],
            log_pi = log_pi,
            shooting = shooting[i],
            method = method[i],
            ts_id = ts_id[i]
          )
        }

        if (!use_traps_or_snares & !use_traps_and_snares) {
          # no gamma no p
          log_potential_area[i] <- calc_log_potential_area(
            log_rho = log_rho[1:n_method],
            log_gamma = 0,
            p_unique = 0,
            log_effort_per = log_effort_per[i],
            effort_per = effort_per[i],
            n_trap_m1 = n_trap_m1[i],
            log_pi = log_pi,
            shooting = shooting[i],
            method = method[i],
            ts_id = ts_id[i]
          )
        }
      }
      # probability of capture, given that an individual is in the surveyed area
      if (single_method) {
        log_theta[i] <- log(
          ilogit(
            beta1 + inprod(X_p[i, 1:m_p], beta_p[1, 1:m_p])
          )
        ) +
          min(0, log_potential_area[i] - log_survey_area_km2[i])
      } else {
        log_theta[i] <- log(
          ilogit(
            beta1[method[i]] + inprod(X_p[i, 1:m_p], beta_p[method[i], 1:m_p])
          )
        ) +
          min(0, log_potential_area[i] - log_survey_area_km2[i])
      }

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

    if (single_property) {
      N[nH[1, 1]] ~ dunif(n1_min, n1_max)
      # N[nH[1, 1]] ~ dpois(round(lambda_1))

      # population growth across time steps
      for (j in 2:n_time_prop) {
        # loop through every PP, including missing ones

        lambda[nH[1, j - 1]] <- (N[nH[1, j - 1]] - rem[1, j - 1]) *
          zeta /
          2 +
          (N[nH[1, j - 1]] - rem[1, j - 1]) * phi[nH[1, j - 1]]

        N[nH[1, j]] ~ dpois(lambda[nH[1, j - 1]])
        phi[nH[1, j - 1]] ~ dbeta(a_phi, b_phi)
      }
    } else {
      for (i in 1:n_property) {
        N[nH[i, 1]] ~ dunif(n1_min[i], n1_max[i])
        # N[nH[i, 1]] ~ dpois(round(lambda_1[i]))

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
    shooting = double(0),
    method = double(0),
    ts_id = double(0)
  ) {
    if (shooting == 1) {
      # firearms, fixed wing, and helicopter
      log_potential_area <- log_rho[method] + log_effort_per
    } else {
      # traps and snares
      log_potential_area <- log_pi +
        (2 *
          (log_rho[ts_id] +
            log_effort_per -
            log(exp(log_gamma[ts_id]) + effort_per))) +
        log(1 + (p_unique[ts_id] * n_trap_m1))
    }
    return(log_potential_area)
    returnType(double(0))
  },
  buildDerivs = TRUE
)
