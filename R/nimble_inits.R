#' Prepare initial values for nimble model
#' @param constants_nimble A list of constants for the nimble model.
#' @param data_nimble A list of data for the nimble model.
#' @param buffer A numeric value to add to the initial values of N and lambda_1 to ensure they comply with take values.
#' @return A list of initial values for the nimble model. Works for one chain.
#' @export

nimble_inits <- function(constants_nimble, data_nimble, buffer = 1000) {
	with(append(constants_nimble, data_nimble), {
		beta1 <- rnorm(n_method, beta1_mu, sqrt(1 / beta1_tau))
		beta_p <- matrix(
			rnorm(n_betaP, beta_p_mu, sqrt(1 / beta_p_tau)),
			n_method,
			m_p,
			byrow = TRUE
		)
		p_mu <- rnorm(2, p_mu_mu, sqrt(1 / p_mu_tau))
		log_gamma <- rnorm(2, log_gamma_mu, sqrt(1 / log_gamma_tau))
		log_rho <- rnorm(n_method, log_rho_mu, sqrt(1 / log_rho_tau))
		psi_phi <- rgamma(1, shape = psi_shape, rate = psi_rate)
		phi_mu <- rbeta(1, phi_mu_a, phi_mu_b)
		log_nu <- rnorm(1, log_nu_mu, sqrt(1 / log_nu_tau))
		nu <- exp(log_nu)

		a <- phi_mu * psi_phi
		b <- (1 - phi_mu) * psi_phi
		zeta <- pp_len * nu / 365
		N <- phi <- lambda <- rep(NA, max(nH, na.rm = TRUE))
		n_init <- rep(NA, n_property)
		for (i in 1:n_property) {
			n_init[i] <- round(runif(1, n1_max[i] * 0.25, n1_max[i] * 0.75)) +
				sum(rem[i, ], na.rm = TRUE)
			N[nH[i, 1]] <- n_init[i]
			for (j in 2:n_time_prop[i]) {
				phi[nH[i, j - 1]] <- max(0.05, min(rbeta(1, a, b), 0.95))
				z <- N[nH[i, j - 1]] - rem[i, j - 1]
				z <- max(1, z)
				lambda[nH[i, j - 1]] <- z * zeta / 2 + z * phi[nH[i, j - 1]]

				N[nH[i, j]] <- rpois(1, lambda[nH[i, j - 1]])
			}
		}

		list(
			lambda_1 = pmin(n_init + buffer, n1_max),
			beta_p = beta_p,
			beta1 = beta1,
			p_mu = p_mu,
			p_unique = boot::inv.logit(p_mu),
			phi_mu = phi_mu,
			psi_phi = psi_phi,
			a_phi = a,
			b_phi = b,
			N = N + buffer,
			log_nu = log_nu,
			log_gamma = log_gamma,
			log_rho = log_rho,
			phi = phi,
			zeta = zeta,
			log_zeta = log(zeta)
		)
	})
}
