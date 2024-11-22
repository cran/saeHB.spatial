#' Small Area Estimation under Spatial Simultaneous Autoregressive (SAR) Model and Normal Distribution using Hierarchical Bayesian Method
#'
#' @description This function gives small area estimator under Spatial SAR Model and is implemented to variable of interest (y) that assumed to be a Normal Distribution. The range of data is \eqn{(-\infty < y < \infty)}.
#'
#' @param formula formula that describe the fitted model.
#' @param vardir sampling variances of direct estimations.
#' @param proxmat \code{D*D} proximity matrix with values in the interval \code{[0,1]} containing the proximities between the row and column domains. The rows add up to \code{1}.
#' @param iter.update number of updates with default \code{3}.
#' @param iter.mcmc number of total iterations per chain with default \code{2000}.
#' @param thin thinning rate, must be a positive integer with default \code{1}.
#' @param burn.in number of iterations to discard at the beginning with default \code{1000}.
#' @param coef optional vector containing the mean of the prior distribution of the regression model coefficients.
#' @param var.coef optional vector containing the variances of the prior distribution of the regression model coefficients.
#' @param data the data frame.
#'
#'
#' @returns This function returns a list of the following objects:
#' \item{Est}{A data frame of Small Area mean Estimates using Hierarchical Bayesian Method}
#' \item{refVar}{Estimated random effect variances}
#' \item{coefficient}{A data frame with estimated model coefficient}
#' \item{plot}{Trace, Density, and Autocorrelation Function Plot of MCMC samples}
#'
#' @examples
#' ## For data without any non-sampled area
#' data(sp.norm)       # Load dataset
#' data(prox.mat)      # Load proximity Matrix
#'
#' \donttest{result <- sar.normal(y ~ x1 + x2, "vardir", prox.mat, data = sp.norm)}
#'
#' \donttest{result$Est}          # Small Area mean Estimates
#' \donttest{result$refVar}       # Estimated random effect variances
#' \donttest{result$coefficient}  # Estimated model coefficient
#'
#'
#' # Load library 'coda' to execute the plot
#' # autocorr.plot(result$plot[[3]])    # Generate ACF Plot
#' # plot(result$plot[[3]])             # Generate Density and Trace plot
#'
#'
#'
#' ## For data with non-sampled area use sp.normNs
#'
#'
#' @import stats
#' @import rjags
#' @import stringr
#' @import grDevices
#' @import graphics
#' @import coda
#'
#' @export
sar.normal <- function(formula, vardir, proxmat, iter.update = 3, iter.mcmc = 2000, thin = 1, burn.in = 1000, coef, var.coef, data) {
  result <- list(Est = NA, refVar = NA, coefficient = NA, plot = NA)
  formuladata <- model.frame(formula, data, na.action = NULL)

  if (any(is.na(formuladata[, -1]))) stop("Auxiliary Variables contains NA values.")

  W <- as.matrix(proxmat)
  if (any(is.na(W))) stop("Proximity Matrix contains NA values.")
  if (nrow(W) != nrow(formuladata) | ncol(W) != nrow(formuladata))
    stop("Proximity Matrix must be D*D matrix")

  auxVar <- as.matrix(formuladata[, -1])
  nvar <- ncol(auxVar) + 1
  formuladata <- data.frame(formuladata, vardir = data[, vardir])

  if (!missing(coef)){
    if (length(coef) != nvar){
      stop("Length of vector coef does not match the number of model coefficients, the length must be ", nvar)
    }
    mu.b = coef
  } else {
    mu.b = rep(0, nvar)
  }

  if (!missing(var.coef)){
    if (length(var.coef) != nvar){
      stop("Length of vector var.coef does not match the number of model coefficients, the length must be ", nvar)
    }
    tau.b = 1/var.coef
  } else {
    tau.b = rep(1, nvar)
  }

  if (iter.update < 3){
    stop("The number of iteration updates at least 3 times")
  }

  # Tersampel
  if (!any(is.na(formuladata[, 1]))) {
    formuladata <- as.matrix(na.omit(formuladata))
    x <- model.matrix(formula, data = as.data.frame(formuladata))
    n <- nrow(formuladata)

    tau.ua = tau.ub = 1
    a.var = rep(1, n)
    I <- diag(1, n)
    O <- rep(0, n)

    for (i in 1:iter.update) {
      dat <- list(n = n, nvar = nvar, y = formuladata[, 1], x = as.matrix(x[, -1]), vardir = formuladata[, (nvar + 1)],
                  mu.b = mu.b, tau.b = tau.b, tau.ua = tau.ua, tau.ub = tau.ub, I = I, O = O, W = W)
      inits <- list(v = rep(0, n), b = mu.b, tau.u = 1, rho = 0)

      cat("model {
          for (i in 1:n) {
		          y[i] ~ dnorm(mu[i], tau[i])
		          mu[i] <- b[1] + sum(b[2:nvar]*x[i,]) + v[i]
		          tau[i] <- 1/vardir[i]
		          a.var[i] <- sig.v[i, i]
          }

          C <- (I - rho * W)
	        tau.v <- tau.u * t(C) %*% C
	        sig.v <- inverse(tau.v)

	        v ~ dmnorm(O, tau.v)

	        for (k in 1:nvar) {
	            b[k] ~ dnorm(mu.b[k],tau.b[k])
	        }

	        tau.u ~ dgamma(tau.ua, tau.ub)

	        rho ~ dunif(-0.9999, 0.9999)

			}", file = "saeHBspatial.txt")

      jags.m <- jags.model(file = "saeHBspatial.txt", data = dat, inits = inits, n.chains = 1, n.adapt = 500)
      file.remove("saeHBspatial.txt")
      params <- c("mu", "a.var", "b", "tau.u", "rho")
      samps <- coda.samples(jags.m, params, n.iter = iter.mcmc, thin = thin)
      samps1 <- window(samps, start = burn.in + 1, end = iter.mcmc)

      result_samps = summary(samps1)
      a.var = result_samps$statistics[1:n]
      beta = result_samps$statistics[(n + 1):(n + nvar), 1:2]

      for (i in 1:nvar) {
        mu.b[i] = beta[i, 1]
        tau.b[i] = 1/(beta[i, 2]^2)
      }
      tau.ua = result_samps$statistics[2*n + nvar + 2, 1]^2/result_samps$statistics[2*n + nvar + 2, 2]^2
      tau.ub = result_samps$statistics[2*n + nvar + 2, 1]/result_samps$statistics[2*n + nvar + 2, 2]^2
    }

    result_samps = summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i - 1)
      b.varnames[i] <- str_replace_all(paste("b[", idx.b.varnames, "]"), pattern = " ", replacement = "")
    }
    b.varnames[(nvar + 1)] <- "rho"

    result_mcmc <- samps1[, c((n + 1):(n + nvar), (2*n + nvar + 1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var = result_samps$statistics[1:n]

    beta = result_samps$statistics[(n + 1):(n + nvar), 1:2]
    rho = result_samps$statistics[2*n + nvar + 1, 1:2]
    coef = rbind(beta, rho)
    rownames(coef) <- b.varnames

    mu = result_samps$statistics[(n + nvar + 1):(2*n + nvar), 1:2]
    Estimation = data.frame(mu)

    Quantiles <- as.data.frame(result_samps$quantiles)
    q_beta <- (Quantiles[(n + 1):(n + nvar), ])
    q_rho <- Quantiles[2*n + nvar + 1, ]
    q_coef <- rbind(q_beta, q_rho)
    rownames(q_coef) <- b.varnames
    coef <- cbind(coef, q_coef)

    q_mu <- Quantiles[(n + nvar + 1):(2*n + nvar), ]
    Estimation <- data.frame(Estimation, q_mu)
    colnames(Estimation) <- c("MEAN", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

  } else { # Tidak Tersampel
    formuladata <- as.data.frame(formuladata)
    x <- as.matrix(formuladata[, 2:nvar])
    n <- nrow(formuladata)

    tau.ua = tau.ub = 1
    a.var = rep(1, n)
    I <- diag(1, n)
    O <- rep(0, n)

    formuladata$idx <- rep(1:n)
    data_sampled <- na.omit(formuladata)
    data_nonsampled <- formuladata[-data_sampled$idx, ]
    idx_sampled = data_sampled$idx
    idx_nonsampled = data_nonsampled$idx

    for (i in 1:iter.update) {
      dat <- list(nvar = nvar, y = formuladata[, 1], x = as.matrix(formuladata[, 2:nvar]), vardir = formuladata$vardir,
                  idx_sampled = idx_sampled, idx_nonsampled = idx_nonsampled,
                  mu.b = mu.b, tau.b = tau.b, tau.ua = tau.ua, tau.ub = tau.ub, I = I, O = O, W = W)
      inits <- list(v = rep(0, n), b = mu.b, tau.u = 1, rho = 0)

      cat("model {
          for (i in idx_sampled) {
              y[i] ~ dnorm(mu[i], tau[i])
              mu[i] <- b[1] + sum(b[2:nvar]*x[i,]) + v[i]
              tau[i] <- 1/vardir[i]
              a.var[i] <- sig.v[i, i]
          }

          for (j in idx_nonsampled) {
              mu[j] <- mu.b[1] + sum(mu.b[2:nvar]*x[j,]) + v[j]
              a.var[j] <- sig.v[j, j]
          }

          C <- (I - rho * W)
	        tau.v <- tau.u * t(C) %*% C
	        sig.v <- inverse(tau.v)

	        v ~ dmnorm(O, tau.v)

          for (k in 1:nvar) {
              b[k] ~ dnorm(mu.b[k], tau.b[k])
          }

          tau.u ~ dgamma(tau.ua, tau.ub)

	        rho ~ dunif(-0.9999, 0.9999)

      }", file = "saeHBspatial.txt")

      jags.m <- jags.model(file = "saeHBspatial.txt", data = dat, inits = inits, n.chains = 1, n.adapt = 500)
      file.remove("saeHBspatial.txt")
      params <- c("mu", "a.var", "b", "tau.u", "rho")
      samps <- coda.samples(jags.m, params, n.iter = iter.mcmc, thin = thin)
      samps1 <- window(samps, start = burn.in + 1, end = iter.mcmc)

      result_samps = summary(samps1)
      a.var = result_samps$statistics[1:n]
      beta = result_samps$statistics[(n + 1):(n + nvar), 1:2]

      for (i in 1:nvar) {
        mu.b[i] = beta[i, 1]
        tau.b[i] = 1/(beta[i, 2]^2)
      }
      tau.ua = result_samps$statistics[2*n + nvar + 2, 1]^2/result_samps$statistics[2*n + nvar + 2, 2]^2
      tau.ub = result_samps$statistics[2*n + nvar + 2, 1]/result_samps$statistics[2*n + nvar + 2, 2]^2
    }

    result_samps = summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i - 1)
      b.varnames[i] <- str_replace_all(paste("b[", idx.b.varnames, "]"), pattern = " ", replacement = "")
    }
    b.varnames[(nvar + 1)] <- "rho"

    result_mcmc <- samps1[, c((n + 1):(n + nvar), (2*n + nvar + 1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var = result_samps$statistics[1:n]

    beta = result_samps$statistics[(n + 1):(n + nvar), 1:2]
    rho = result_samps$statistics[2*n + nvar + 1, 1:2]
    coef = rbind(beta, rho)
    rownames(coef) <- b.varnames

    mu = result_samps$statistics[(n + nvar + 1):(2*n + nvar), 1:2]
    Estimation = data.frame(mu)

    Quantiles <- as.data.frame(result_samps$quantiles)
    q_beta <- (Quantiles[(n + 1):(n + nvar), ])
    q_rho <- Quantiles[2*n + nvar + 1, ]
    q_coef <- rbind(q_beta, q_rho)
    rownames(q_coef) <- b.varnames
    coef <- cbind(coef, q_coef)

    q_mu <- Quantiles[(n + nvar + 1):(2*n + nvar), ]
    Estimation <- data.frame(Estimation, q_mu)
    colnames(Estimation) <- c("MEAN", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

  }

  result$Est = Estimation
  result$refVar = a.var
  result$coefficient = coef
  result$plot = list(graphics.off(), par(mar = c(2, 2, 2, 2)),
                     autocorr.plot(result_mcmc, col = "brown2", lwd = 2),
                     plot(result_mcmc, col = "brown2", lwd = 2))
  on.exit(par(result$plot[[2]]))
  return(result)
}


