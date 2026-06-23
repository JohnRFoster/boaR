#' Prepare initial values for nimble model
#' @param constants_nimble A list of constants for the nimble model.
#' @param data_nimble A list of data for the nimble model.
#' @param buffer A numeric value to add to the initial values of N and lambda_1 to ensure they comply with take values.
#' @param beta1 vector for initial values of beta1. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @param betap vector for initial values of beta_p. This matrix is specified by row. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @param p_mu vector for initial values of p_mu. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @param log_gamma vector for initial values of log_gamma. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @param log_rho vector for initial values of log_rho. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @param psi_phi range vector (min, max) for the bounds of the initial value for psi_phi. A random value will be drawn between this range. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @param phi_mu range vector (min, max) for the bounds of the initial value for phi_mu. A random value will be drawn between this range. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @param log_nu range vector (min, max) for the bounds of the initial value for log_nu. A random value will be drawn between this range. If NULL (default) uses prior distribution parameters specfied in constants_nimble.
#' @return A list of initial values for the nimble model. Works for one chain.
#' @export

nimble_inits <- function(
	constants_nimble,
	data_nimble,
	buffer = 1000,
	beta1 = NULL,
	beta_p = NULL,
	p_mu = NULL,
	log_gamma = NULL,
	log_rho = NULL,
	psi_phi = NULL,
	phi_mu = NULL,
	log_nu = NULL
) {
	with(append(constants_nimble, data_nimble), {
		if (is.null(beta1)) {
			beta1 <- rnorm(n_method, beta1_mu, sqrt(1 / beta1_tau))
		} else {
			beta1 <- rnorm(length(beta1), beta1, 0.25)
		}

		if (is.null(beta_p)) {
			beta_p <- matrix(
				rnorm(n_betaP, beta_p_mu, sqrt(1 / beta_p_tau)),
				n_method,
				m_p,
				byrow = TRUE
			)
		} else {
			beta_p <- matrix(
				rnorm(length(beta_p), beta_p, 0.25),
				n_method,
				m_p,
				byrow = TRUE
			)
		}

		if (is.null(p_mu)) {
			p_mu <- rnorm(length(p_mu_mu), p_mu_mu, sqrt(1 / p_mu_tau))
		} else {
			p_mu <- rnorm(length(p_mu), p_mu, 0.1)
		}

		if (is.null(log_gamma)) {
			log_gamma <- rnorm(
				length(log_gamma_mu),
				log_gamma_mu,
				sqrt(1 / log_gamma_tau)
			)
		} else {
			log_gamma <- rnorm(
				length(log_gamma),
				log_gamma,
				0.1
			)
		}

		if (is.null(log_rho)) {
			log_rho <- rnorm(n_method, log_rho_mu, sqrt(1 / log_rho_tau))
		} else {
			log_rho <- rnorm(n_method, log_rho, 0.1)
		}

		if (is.null(psi_phi)) {
			psi_phi <- rgamma(1, shape = psi_shape, rate = psi_rate)
		} else {
			psi_phi <- runif(1, psi_phi[1], psi_phi[2])
		}

		if (is.null(phi_mu)) {
			phi_mu <- rbeta(1, phi_mu_a, phi_mu_b)
		} else {
			phi_mu <- runif(1, phi_mu[1], phi_mu[2])
		}

		if (is.null(log_nu)) {
			log_nu <- rnorm(1, log_nu_mu, sqrt(1 / log_nu_tau))
		} else {
			log_nu <- runif(1, log_nu[1], log_nu[2])
		}

		a <- phi_mu * psi_phi
		b <- (1 - phi_mu) * psi_phi

		N <- phi <- lambda <- rep(NA, max(nH, na.rm = TRUE))
		n_init <- rep(NA, n_property)
		for (i in 1:n_property) {
			n_init[i] <- round(runif(1, n1_max[i] * 0.25, n1_max[i] * 0.75)) +
				sum(rem[i, ], na.rm = TRUE)
			N[nH[i, 1]] <- n_init[i]
			for (j in 2:n_time_prop[i]) {
				phi[nH[i, j - 1]] <- max(0.5, min(rbeta(1, a, b), 0.99))
				z <- N[nH[i, j - 1]] - rem[i, j - 1]
				z <- max(1, z)

				nu <- exp(log_nu)
				zeta <- pp_len * nu / 365

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
			lambda = lambda[1:(n_time_prop - 1)] + buffer,
			N = N + buffer,
			log_nu = log_nu,
			log_gamma = log_gamma,
			log_rho = log_rho,
			phi = phi[1:(n_time_prop - 1)],
			zeta = zeta,
			log_zeta = log(zeta)
		)
	})
}
