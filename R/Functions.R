
#' Layer Mean of the log-affine local Pareto Distribution

#' @description  Calculates the expected loss of a log-affine local Pareto
#'               distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha_0 Numeric. Local Pareto alpha at t.
#' @param t Numeric. Threshold of the distribution.
#' @param gamma Numeric. Parameter \code{gamma}.
#' @param delta Numeric. Parameter \code{delta}. Alternative parameter that can be used instead of \code{gamma}.
#'        \code{delta = alpha_0 * gamma * log(2)} is the increase in the local Pareto alpha when doubling the amount.
#'
#' @return The expected loss of the log-affine local Pareto distribution with parameters \code{t}, \code{alpha_0} and \code{gamma} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' LALocPareto_Layer_Mean(4000, 1000, t = 1000, alpha_0 = 2, gamma = 1)
#' LALocPareto_Layer_Mean(4000, 1000, t = 1000, alpha_0 = 2, delta = 1)
#'
#' @export


LALocPareto_Layer_Mean <- function(Cover, AttachmentPoint, t, alpha_0, gamma = NULL, delta = NULL) {
  if (!is.nonnegative.vector(Cover)) {
    warning("Cover must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.nonnegative.finite.vector(AttachmentPoint)) {
    warning("AttachmentPoint must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(alpha_0)) {
    warning("alpha_0 must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(gamma) && is.positive.finite.vector(delta)) {
    gamma <- delta / alpha_0 / log(2)
  }
  if (!is.positive.finite.vector(gamma)) {
    warning("gamma or delta must be a positive vector.")
    return(NaN)
  }
  delta <- alpha_0 * gamma * log(2)

  vecfun <- Vectorize(LALocPareto_Layer_Mean_s)
  vecfun(Cover, AttachmentPoint, t, alpha_0, gamma, delta)
}



LALocPareto_Layer_Mean_s <- function(Cover, AttachmentPoint, t, alpha_0, gamma, delta) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_0)) {
    warning("alpha must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!valid.parameters.LALocPareto(t, alpha_0, gamma, delta)) {
    warning("Parameter of LALocPareto not valid.")
    return(NaN)
  }

  SF <- function(x) {
    res <- ifelse(x > t, (t/x)^(alpha_0 + 0.5 * alpha_0 * gamma * log(x/t)), 1)
    return(res)
  }
  EP <- integrate(SF, lower = AttachmentPoint, upper = AttachmentPoint + Cover)$value

  temp <- max(AttachmentPoint, t)
  for (i in 1:20) {
    if (temp * 10^i < Cover) {
      EP_temp <- integrate(SF, lower = AttachmentPoint, upper = AttachmentPoint + temp * 10^i)$value
      EP <- max(EP, EP_temp)
    }
  }
  return(EP)
}


#' Second layer moment of the log-affine local Pareto Distribution

#' @description  Calculates the second moment of the loss of a log-affine local Pareto
#'               distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha_0 Numeric. Local Pareto alpha at t.
#' @param t Numeric. Threshold of the distribution.
#' @param gamma Numeric. Parameter \code{gamma}.
#' @param delta Numeric. Parameter \code{delta}. Alternative parameter that can be used instead of \code{gamma}.
#'        \code{delta = alpha_0 * gamma * log(2)} is the increase in the local Pareto alpha when doubling the amount.
#'
#' @return The second moment of the loss of the log-affine local Pareto distribution with parameters \code{t}, \code{alpha_0} and \code{gamma} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' LALocPareto_Layer_SM(4000, 1000, t = 1000, alpha_0 = 2, gamma = 1)
#' LALocPareto_Layer_SM(4000, 1000, t = 1000, alpha_0 = 2, delta = 1)
#'
#' @export


LALocPareto_Layer_SM <- function(Cover, AttachmentPoint, t, alpha_0, gamma = NULL, delta = NULL) {
  if (!is.nonnegative.vector(Cover)) {
    warning("Cover must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.nonnegative.finite.vector(AttachmentPoint)) {
    warning("AttachmentPoint must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(alpha_0)) {
    warning("alpha_0 must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(gamma) && is.positive.finite.vector(delta)) {
    gamma <- delta / alpha_0 / log(2)
  }
  if (!is.positive.finite.vector(gamma)) {
    warning("gamma or delta must be a positive vector.")
    return(NaN)
  }
  delta <- alpha_0 * gamma * log(2)


  vecfun <- Vectorize(LALocPareto_Layer_SM_s)
  vecfun(Cover, AttachmentPoint, t, alpha_0, gamma, delta)
}



LALocPareto_Layer_SM_s <- function(Cover, AttachmentPoint, t, alpha_0, gamma, delta) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.positive.finite.number(alpha_0)) {
    warning("alpha must be a positive number.")
    return(NaN)
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!valid.parameters.LALocPareto(t, alpha_0, gamma, delta)) {
    warning("Parameter of LALocPareto not valid.")
    return(NaN)
  }

  Integrand <- function(x) {
    res <- 2 * (x - AttachmentPoint) * ifelse(x > t, 1 - (t/x)^(alpha_0 + 0.5 * alpha_0 * gamma * log(x/t)), 0)
    return(res)
  }
  Int <- integrate(Integrand, lower = AttachmentPoint, upper = AttachmentPoint + Cover)$value

  temp <- max(AttachmentPoint, t)
  for (i in 1:20) {
    if (temp * 10^i < Cover) {
      Int_temp <- integrate(Integrand, lower = AttachmentPoint, upper = AttachmentPoint + temp * 10^i)$value
      Int <- max(Int, Int_temp)
    }
  }
  CDF <- function(x) {1 - (t/x)^(alpha_0 + 0.5 * alpha_0 * gamma * log(x/t))}

#  EP <- (AttachmentPoint + Cover)^2 * CDF(AttachmentPoint + Cover) - AttachmentPoint^2 * CDF(AttachmentPoint) +
#              Cover^2 * (1 - CDF(AttachmentPoint + Cover)) - Int - 2 * AttachmentPoint * (AttachmentPoint + Cover) * CDF(AttachmentPoint + Cover) +
#              AttachmentPoint^2 * CDF(AttachmentPoint)
  EP <- Cover^2 - Int
  return(EP)
}


#' Layer variance of the log-affine local Pareto Distribution

#' @description  Calculates the variance of the loss of a log-affine local Pareto
#'               distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha_0 Numeric. Local Pareto alpha at t.
#' @param t Numeric. Threshold of the distribution.
#' @param gamma Numeric. Parameter \code{gamma}.
#' @param delta Numeric. Parameter \code{delta}. Alternative parameter that can be used instead of \code{gamma}.
#'        \code{delta = alpha_0 * gamma * log(2)} is the increase in the local Pareto alpha when doubling the amount.
#'
#' @return The variance of the loss of the log-affine local Pareto distribution with parameters \code{t}, \code{alpha_0} and \code{gamma} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' LALocPareto_Layer_Var(4000, 1000, t = 1000, alpha_0 = 2, gamma = 1)
#' LALocPareto_Layer_Var(4000, 1000, t = 1000, alpha_0 = 2, delta = 1)
#'
#' @export


LALocPareto_Layer_Var <- function(Cover, AttachmentPoint, t, alpha_0, gamma = NULL, delta = NULL) {
  if (!is.nonnegative.vector(Cover)) {
    warning("Cover must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.nonnegative.finite.vector(AttachmentPoint)) {
    warning("AttachmentPoint must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(alpha_0)) {
    warning("alpha_0 must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(gamma) && is.positive.finite.vector(delta)) {
    gamma <- delta / alpha_0 / log(2)
  }
  if (!is.positive.finite.vector(gamma)) {
    warning("gamma or delta must be a positive vector.")
    return(NaN)
  }
  delta <- alpha_0 * gamma * log(2)


  LALocPareto_Layer_SM(Cover, AttachmentPoint, t, alpha_0, gamma, delta) - LALocPareto_Layer_Mean(Cover, AttachmentPoint, t, alpha_0, gamma, delta)^2
}



#' Simulation of the log-affine local Pareto Distribution
#'
#' @description Generates random deviates of a log-affine local Pareto distribution
#'
#' @param n Numeric. Number of observations.
#' @param t Numeric vector. Thresholds of the log-affine local Pareto distributions
#' @param alpha_0 Numeric vector. Local Pareto alphas of the log-affine local Pareto distributions at \code{t}.
#' @param gamma NULL or Numeric vector. If \code{gamma} is \code{NULL} then \code{delta} is used.
#' @param delta NULL or Numeric vector. If \code{gamma} is \code{NULL} then \code{delta} is used.
#'
#' @return A vector of \code{n} samples from the log-affine local Pareto distribution
#'
#' @examples
#' rLALocPareto(100, 1000, 2, 1)
#' rLALocPareto(100, 1000, 2, delta = 0.1)
#'
#' @export

rLALocPareto <- function(n, t, alpha_0, gamma = NULL, delta = NULL) {
  if (!is.positive.finite.number(n)) {
    warning("n must be a positive number.")
    return(NaN)
  }
  n <- ceiling(n)
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (n %% length(t) != 0) {
    warning("n is not a multiple of length(t).")
    t <- rep(t, ceiling(n / length(t)))[1:n]
  }
  if (!is.positive.finite.vector(alpha_0)) {
    warning("alpha must be a positive vector.")
    return(NaN)
  }
  if (n %% length(alpha_0) != 0) {
    warning("n is not a multiple of length(alpha_0).")
    alpha_0 <- rep(alpha_0, ceiling(n / length(alpha_0)))[1:n]
  }

  if (is.null(gamma)) {
    gamma = delta / alpha_0 / log(2)
  }


  FinvLALocPareto <- function(x, t, alpha_0, gamma) {
    return(t * exp((-alpha_0 + sqrt(alpha_0^2 - 2 * alpha_0 * gamma * log(1-x))) / (alpha_0 * gamma)))
  }



  sim_losses <- numeric(n)
  sim_losses <- FinvLALocPareto(stats::runif(n, 0, 1), t, alpha_0, gamma)

  return(sim_losses)
}




#' Maximum Likelihood Estimation for the log-affine local Pareto distribution
#'
#' @description Calculates the maximum likelihood estimator for the parameters alpha_0 and gamma of a log-affine local Pareto distribution
#' with a known threshold
#'
#' @param losses Numeric vector. Losses that are used for the ML estimation.
#' @param t Numeric. Threshold of the Pareto distribution.
#' @param gamma Numeric vector or NULL. If NULL then the parameter \code{gamma} is estimated.
#' @param weights Numeric vector or NULL. Weights for the losses. For instance \code{weights[i] = 2} and \code{weights[j] = 1} for \code{j != i} has the same effect as adding another loss of size \code{loss[i]}.
#'
#' @return Maximum likelihood estimator for the parameters \code{alpha_0} and \code{gamma} of a log-affine local Pareto distribution with threshold \code{t} given the observations \code{losses}
#'
#' @examples
#' losses <- rLALocPareto(100, 1000, 2, 2)
#' LALocPareto_ML_Estimator(losses, 1000)
#'
#' @export

LALocPareto_ML_Estimator <- function(losses, t, gamma = NULL, delta = NULL, weights = NULL) {
  if (!is.nonnegative.finite.vector(losses)) {
    warning("losses must be non-negative.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be positive.")
    return(NaN)
  }
  if (length(t) != 1 && length(t) != length(losses)) {
    warning("t must have length 1 or same length as losses.")
    return(NaN)
  }
  if (!is.null(gamma)) {
    if (length(gamma) != 1 && length(gamma) != length(losses)) {
      warning("gamma must be NULL or have length 1 or same length as losses.")
      return(NaN)
    }
  }
  if (!is.null(delta)) {
    if (length(delta) != 1 && length(delta) != length(losses)) {
      warning("delta must be NULL or have length 1 or same length as losses.")
      return(NaN)
    }
  }
  if (!is.null(gamma) & !is.null(delta)) {
    warning("Either gamma or delta must be NULL.")
    return(NaN)
  }


  if (length(t) == 1) {
    t <- rep(t, length(losses))
  }
  if (is.null(weights)) weights <- rep(1, length(losses))
  if (!is.positive.finite.vector(weights)) {
    warning("weights must NULL or positive.")
    return(NaN)
  }
  if (length(weights) != length(losses)) {
    warning("weights must have the same length as losses.")
    return(NaN)
  }

  index <- losses > t
  if (sum(index) == 0) {
    warning("No loss is larger than t and the specific reporting_threshold.")
    return(NaN)
  }
  losses <- losses[index]
  weights <- weights[index]
  t <- t[index]
  if (!is.null(gamma)) {
    if (length(gamma) != 1) {
      gamma <- gamma[index]
    }
  }

  n <- length(losses)

  if (!is.null(gamma)) {
    alpha_0_hat <- n / sum(log(losses / t) + gamma / 2 * log(losses / t)^2)
    MLE <- list(alpha_0_hat = alpha_0_hat, gamma_hat = gamma, delta_hat = gamma * alpha_0_hat * log(2), gamma = gamma, delta = NULL)

    FisherInformation <- matrix(NA, 2, 2)
    FisherInformation[1,1]= n / alpha_0_hat^2
    rownames(FisherInformation) <- c("alpha_0", "gamma")
    colnames(FisherInformation) <- c("alpha_0", "gamma")

    #MLE$ObservedFisherInformation <- FisherInformation

    ApproxCovMatrix <- FisherInformation
    ApproxCovMatrix[1,1] <- alpha_0_hat^2 / n
    ApproxCovMatrix[2,1] <- 0
    ApproxCovMatrix[1,2] <- 0
    ApproxCovMatrix[2,2] <- 0
    MLE$ApproxCovMatrix <- ApproxCovMatrix

    Jacobi <- matrix(0, 2, 2)
    Jacobi[1,1] <- 1
    Jacobi[1,2] <- 0
    Jacobi[2,1] <- gamma * log(2)
    Jacobi[2,2] <- alpha_0_hat * log(2)
    ApproxCovMatrix2 <- Jacobi %*% ApproxCovMatrix %*% t(Jacobi)
    rownames(ApproxCovMatrix2) <- c("alpha_0", "delta")
    colnames(ApproxCovMatrix2) <- c("alpha_0", "delta")

    MLE$ApproxCovMatrix2 <- ApproxCovMatrix2

    return(MLE)
  } else if (!is.null(delta)) {
    a <- sum(log(losses / t))
    b <- sum(log(losses / t)^2)
    nll <- function(alpha_0) {
      alpha_0 * a + delta / 2 * b - sum(log(alpha_0 + delta * log(losses / t) / log(2)))
    }

    optim_result <- NULL
    try(optim_result <- stats::optim(1, nll, method = "Brent", lower = 0, upper = 100))
    if (is.null(optim_result) || optim_result$convergence != 0) {
      warning("optimization did not converge")
      return(NaN)
    }
    alpha_0_hat <- optim_result$par
    MLE <- list(alpha_0_hat = alpha_0_hat, gamma_hat = delta / alpha_0_hat / log(2), delta_hat = delta, gamma = NULL, delta = delta)

    FisherInformation <- matrix(NA, 2, 2)
    FisherInformation[1,1]= sum(1 / (alpha_0_hat + delta * log(losses / t) / log(2))^2)
    rownames(FisherInformation) <- c("alpha_0", "delta")
    colnames(FisherInformation) <- c("alpha_0", "delta")

    #MLE$ObservedFisherInformation <- FisherInformation

    ApproxCovMatrix <- 1 / FisherInformation
    ApproxCovMatrix[1,2] <- 0
    ApproxCovMatrix[2,1] <- 0
    ApproxCovMatrix[2,2] <- 0


    Jacobi <- matrix(0, 2, 2)
    Jacobi[1,1] <- 1
    Jacobi[1,2] <- 0
    Jacobi[2,1] <- -delta / alpha_0_hat^2 / log(2)
    Jacobi[2,2] <- 1 / alpha_0_hat / log(2)
    ApproxCovMatrix2 <- Jacobi %*% ApproxCovMatrix %*% t(Jacobi)
    rownames(ApproxCovMatrix2) <- c("alpha_0", "gamma")
    colnames(ApproxCovMatrix2) <- c("alpha_0", "gamma")

    MLE$ApproxCovMatrix <- ApproxCovMatrix2
    MLE$ApproxCovMatrix2 <- ApproxCovMatrix



    return(MLE)
  } else {
    a <- sum(log(losses / t))
    b <- sum(log(losses / t)^2)
    nlambda <- function(gamma) {
      n * log(a + gamma / 2 * b) - sum(log(1 + gamma * log(losses / t)))
    }

    optim_result <- NULL
    try(optim_result <- stats::optim(1, nlambda, method = "Brent", lower = 0, upper = 100))
    if (is.null(optim_result) || optim_result$convergence != 0) {
      warning("optimization did not converge")
      return(NaN)
    }
    gamma_hat <- optim_result$par
    alpha_0_hat <- n / sum(log(losses / t) + gamma_hat / 2 * log(losses / t)^2)
    MLE <- list(alpha_0_hat = alpha_0_hat, gamma_hat = gamma_hat, delta_hat = gamma_hat * alpha_0_hat * log(2), gamma = NULL, delta = NULL)

    FisherInformation <- matrix(NA, 2, 2)
    FisherInformation[1,1] = n / alpha_0_hat^2
    FisherInformation[1,2] = 0.5 * sum(log(losses / t)^2)
    FisherInformation[2,1] = FisherInformation[1,2]
    FisherInformation[2,2] = sum(log(losses / t)^2 / (1 + gamma_hat * log(losses / t))^2)

    rownames(FisherInformation) <- c("alpha_0", "gamma")
    colnames(FisherInformation) <- c("alpha_0", "gamma")

    MLE$ObservedFisherInformation <- FisherInformation
    ApproxCovMatrix <- solve(FisherInformation)
    MLE$ApproxCovMatrix <- ApproxCovMatrix

    Jacobi <- matrix(0, 2, 2)
    Jacobi[1,1] <- 1
    Jacobi[1,2] <- 0
    Jacobi[2,1] <- gamma_hat * log(2)
    Jacobi[2,2] <- alpha_0_hat * log(2)
    ApproxCovMatrix2 <- Jacobi %*% ApproxCovMatrix %*% t(Jacobi)
    rownames(ApproxCovMatrix2) <- c("alpha_0", "delta")
    colnames(ApproxCovMatrix2) <- c("alpha_0", "delta")

    MLE$ApproxCovMatrix2 <- ApproxCovMatrix2

    return(MLE)
  }

}



#' CDF of the log-affine local Pareto Distribution

#' @description  Calculates the CDF of a log-affine local Pareto distribution
#'
#' @param x Numeric. Value, where the CDF is evaluated.
#' @param alpha_0 Numeric. Local Pareto alpha at t.
#' @param t Numeric. Threshold of the distribution.
#' @param gamma Numeric. Parameter \code{gamma}.
#' @param delta Numeric. Parameter \code{delta}. Alternative parameter that can be used instead of \code{gamma}.
#'        \code{delta = alpha_0 * gamma * log(2)} is the increase in the local Pareto alpha when doubling the amount.
#'
#' @return The CDF of the log-affine local Pareto distribution with parameters \code{t}, \code{alpha_0} and \code{gamma} at \code{x}.
#'
#' @examples
#' pLALocPareto(2000, t = 1000, alpha_0 = 2, gamma = 1)
#' pLALocPareto((1:10)*1000, t = 1000, alpha_0 = 2, delta = 1)
#'
#' @export


pLALocPareto <- function(x, t, alpha_0, gamma = NULL, delta = NULL) {
  if (!is.vector(x)) {
    warning("x must be a vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(alpha_0)) {
    warning("alpha_0 must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(gamma) && is.positive.finite.vector(delta)) {
    gamma <- delta / alpha_0 / log(2)
  }
  if (!is.positive.finite.vector(gamma)) {
    warning("gamma or delta must be a positive vector.")
    return(NaN)
  }
  delta <- alpha_0 * gamma * log(2)

  SF <- ifelse(x > t, (t/x)^(alpha_0 + 0.5 * alpha_0 * gamma * log(x/t)), 1)
  return(1-SF)
}





#' PDF of the log-affine local Pareto Distribution

#' @description  Calculates the PDF of a log-affine local Pareto distribution
#'
#' @param x Numeric. Value, where the CDF is evaluated.
#' @param alpha_0 Numeric. Local Pareto alpha at t.
#' @param t Numeric. Threshold of the distribution.
#' @param gamma Numeric. Parameter \code{gamma}.
#' @param delta Numeric. Parameter \code{delta}. Alternative parameter that can be used instead of \code{gamma}.
#'        \code{delta = alpha_0 * gamma * log(2)} is the increase in the local Pareto alpha when doubling the amount.
#'
#' @return The PDF of the log-affine local Pareto distribution with parameters \code{t}, \code{alpha_0} and \code{gamma} at \code{x}.
#'
#' @examples
#' dLALocPareto(2000, t = 1000, alpha_0 = 2, gamma = 1)
#' dLALocPareto((1:10)*1000, t = 1000, alpha_0 = 2, delta = 1)
#'
#' @export


dLALocPareto <- function(x, t, alpha_0, gamma = NULL, delta = NULL) {
  if (!is.vector(x)) {
    warning("x must be a vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(alpha_0)) {
    warning("alpha_0 must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(gamma) && is.positive.finite.vector(delta)) {
    gamma <- delta / alpha_0 / log(2)
  }
  if (!is.positive.finite.vector(gamma)) {
    warning("gamma or delta must be a positive vector.")
    return(NaN)
  }
  delta <- alpha_0 * gamma * log(2)

  dens <- ifelse(x > t, (alpha_0 + alpha_0 * gamma * log(x/t)) / x * (t/x)^(alpha_0 + 0.5 * alpha_0 * gamma * log(x/t)), 0)
  return(dens)
}








#' Quantile function of the log-affine local Pareto Distribution

#' @description  Calculates the quantile function of a log-affine local Pareto distribution
#'
#' @param p Numeric. The function evaluates the quantile function at \code{p}.
#' @param alpha_0 Numeric. Local Pareto alpha at t.
#' @param t Numeric. Threshold of the distribution.
#' @param gamma Numeric. Parameter \code{gamma}.
#' @param delta Numeric. Parameter \code{delta}. Alternative parameter that can be used instead of \code{gamma}.
#'        \code{delta = alpha_0 * gamma * log(2)} is the increase in the local Pareto alpha when doubling the amount.
#'
#' @return The quantile function of the log-affine local Pareto distribution with parameters \code{t}, \code{alpha_0} and \code{gamma} at \code{p}.
#'
#' @examples
#' qLALocPareto(0.5, 1000, 2, gamma = 1)
#' qLALocPareto(0.5, 1000, 2, delta = 1)
#'
#' @export


qLALocPareto <- function(p, t, alpha_0, gamma = NULL, delta = NULL) {
  if (!is.nonnegative.finite.vector(p)) {
    warning("p must be a non-negative vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(alpha_0)) {
    warning("alpha_0 must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(t)) {
    warning("t must be a positive vector.")
    return(NaN)
  }
  if (!is.positive.finite.vector(gamma) && is.positive.finite.vector(delta)) {
    gamma <- delta / alpha_0 / log(2)
  }
  if (!is.positive.finite.vector(gamma)) {
    warning("gamma or delta must be a positive vector.")
    return(NaN)
  }

  res <- ifelse(p>=1, Inf, ifelse(p<=0, t, t * exp((-alpha_0 + sqrt(alpha_0^2 - 2 * alpha_0 * gamma * log(1-p))) / (alpha_0 * gamma))))
  return(res)
}


































dPanjer <- function(x, mean, dispersion) {
  if (dispersion == 1) {
    # Poisson distribution
    result <- stats::dpois(x, mean)
  } else if (dispersion < 1) {
    # Binomial distribution
    q <- dispersion
    p <- 1 - q
    # Not every dispersion < 1 can be realized with a binomial distribution. Round up number of trys n and recalculate p and q:
    n <-  ceiling(mean / p)
    p <- mean / n
    q <- 1 - p
    if (abs(dispersion - q) > 0.01) {
      warning(paste0("Dispersion has been adjusted from ", round(dispersion, 2)," to ", round(q, 2), " to obtain a matching binomial distribution."))
    }
    result <- stats::dbinom(x, n, p)
  } else {
    # Negative binomial distribution
    p <- 1 / dispersion
    alpha <- mean / (dispersion - 1)
    result <- stats::dnbinom(x, alpha, p)
  }
  return(result)
}

rPanjer <- function(n, mean, dispersion) {
  if (dispersion == 1) {
    # Poisson distribution
    result <- stats::rpois(n, mean)
  } else if (dispersion < 1) {
    # Binomial distribution
    q <- dispersion
    p <- 1 - q
    # Not every dispersion < 1 can be realized with a binomial distribution. Round up number of trys m and recalculate p and q:
    m <-  ceiling(mean / p)
    p <- mean / m
    q <- 1 - p
    if (abs(dispersion - q) > 0.01) {
      warning(paste0("Dispersion has been adjusted from ", round(dispersion, 2)," to ", round(q, 2), " to obtain a matching binomial distribution."))
    }
    result <- stats::rbinom(n, m, p)
  } else {
    # Negative binomial distribution
    p <- 1 / dispersion
    alpha <- mean / (dispersion - 1)
    result <- stats::rnbinom(n, alpha, p)
  }
  return(result)
}


#' Approximate local Pareto distribution by a piecewise Pareto distribution (old version)
#'
#' @description  Calculates the parameters of an approximating piecewise Pareto distribution (old version)
#'
#' @param alpha Function.
#' @param t Numeric. Threshold of the distribution. If \code{t = 0} then \code{t} is derived from \code{alpha}.
#'
#' @return List with parameters and infos.
#'
#' @examples
#' alpha_0 <- 1.5
#' t <- 1000
#' gamma <- 0.5
#' alpha <- function(x) alpha_0 + alpha_0 * gamma * log(x/t)
#' LocalPareto_2_PiecewisePareto_OLD(t, alpha)
#'
#' @export


LocalPareto_2_PiecewisePareto_OLD <- function(alpha, t = NULL) {

  # check if alpha is a vectorized function
  res <- NULL
  try(res <- alpha(c(0,1)), silent = T)
  if (!is.numeric(res)) {
    alpha_s <- alpha
    alpha <- Vectorize(alpha_s)
    res <- NULL
    try(res <- alpha(c(0,1)), silent = T)
    if (!is.numeric(res)) {
      return("alpha is not a vectorized function")
    }
  }
  if (is.positive.finite.number(t)) {
    x_test <- t * 1.9^(0:20)
  } else {
    x_test <- 1.9^((-20):20)
  }
  if (min(alpha(x_test)) < 0) {
    return("alpha is not non-negative")
  }

  if (is.null(t)) {
    t <- find_t(alpha)
  } else{
    if (!is.nonnegative.finite.number(t)) {
      return("t must be a non-negative number")
    }
  }

  tol_alpha <- 0.01
  min_step <- 1.01
  max_step <- 1.1
  target_prob <- 0.95
  stop_SF <- 1e-6
  min_tail_alpha <- 1.01

  t_PP <- numeric(2000)
  alpha_PP <- numeric(2000)
  SF <- numeric(2000)
  xi <- numeric(2000)
  zeta <- numeric(2000)
  rel_error_lb <- numeric(2000)
  rel_error_ub <- numeric(2000)

  nalpha <- function(x) {-alpha(x)}
  alpha_by_x <- function(x) {alpha(x)/x}

  if (t == 0) {
    min_i <- -33
    max_i <- 3
    for (i in min_i:max_i) {
      res <- NULL
      try(res <- integrate(alpha_by_x, lower = 0, upper = 2^i), silent = T)
      if (!is.numeric(res$value)) {
        return("alpha not admissible")
      }
      if (res$value > -log(0.99)) {
        break
      }
    }
    t_1 <- 2^(i-1)
    SF_t_1 <- exp(- integrate(alpha_by_x, lower = 0, upper = t_1)$value)
    rate <- -log(SF_t_1) / t_1

    t_0 <- t_1 / 1000
    n <- 30
    t_PP[1:(n+1)] <- exp(log(1/1000) / n * n:0) * t_1
    SF[1] <- exp(- integrate(alpha_by_x, lower = 0, upper = t_PP[1])$value)
    for (i in 2:(n+1)) {
      SF[i] <- exp(- integrate(alpha_by_x, lower = t_PP[i-1], upper = t_PP[i])$value) * SF[i-1]
    }
    alpha_PP[1] <- log(1 / SF[2]) / log(t_PP[2] / t_PP[1])
    for (i in 2:n) {
      alpha_PP[i] <- log(SF[i] / SF[i+1]) / log(t_PP[i+1] / t_PP[i])
    }
  } else {
    n <- 0
    t_PP[1] <- t
    SF[1] <- 1
  }




  for (i in (n+1):1999) {
    max_step_temp <- max_step
    max_alpha <- - stats::optim(t_PP[i], nalpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * max_step_temp)$value
    min_alpha <- stats::optim(t_PP[i], alpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * max_step_temp)$value
    while (max_step_temp > min_step & max_alpha - min_alpha > tol_alpha) {
      max_step_temp <- max(min_step,  1 + (max_step_temp - 1) / 2)
      max_alpha <- - stats::optim(t_PP[i], nalpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * max_step_temp)$value
      min_alpha <- stats::optim(t_PP[i], alpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * max_step_temp)$value
    }
    if (max_alpha == min_alpha) {
      t_PP[i+1] <- t_PP[i] * max_step_temp
      alpha_PP[i] <- max_alpha
      SF[i+1] <- SF[i] * max_step_temp^(-max_alpha)

      xi[i] <- (t_PP[i] + t_PP[i+1]) / 2
      zeta[i] <- (t_PP[i] + t_PP[i+1]) / 2
      rel_error_lb[i] <- 0
      rel_error_ub[i] <- 0
    } else {   # then max_alpha > 0
      step <- min(target_prob ^ (- 1 / max_alpha), max_step_temp)
      t_PP[i+1] <- t_PP[i] * step
      factor_SF <- NULL
      try(factor_SF <- exp(-integrate(alpha_by_x, lower = t_PP[i], upper = t_PP[i+1])$value), silent = T)
      if (is.null(factor_SF)) {
        alpha_PP[i] <- alpha_PP[i-1]
        SF[i+1] <- SF[i] * (t_PP[i] / t_PP[i+1])^alpha_PP[i]
        warning("stop criterium not reached!")
        break
      }
      SF[i+1] <- SF[i] * factor_SF
      alpha_PP[i] <- log(SF[i+1] / SF[i]) / log(t_PP[i] / t_PP[i+1])

      xi[i] <- t_PP[i] ^ ((alpha_PP[i] - max_alpha) / (min_alpha - max_alpha)) * t_PP[i+1] ^ ((min_alpha - alpha_PP[i]) / (min_alpha - max_alpha))
      zeta[i] <- t_PP[i] ^ ((alpha_PP[i] - min_alpha) / (max_alpha - min_alpha)) * t_PP[i+1] ^ ((max_alpha - alpha_PP[i]) / (max_alpha - min_alpha))
      rel_error_lb[i] <- (t_PP[i] / zeta[i]) ^ (alpha_PP[i] - min_alpha) - 1
      rel_error_ub[i] <- (xi[i] / t_PP[i]) ^ (max_alpha - alpha_PP[i]) - 1
    }

    # if ((i %% 100) == 0) {
    #   if (exp(1000 / i * log(SF[i+1])) > stop_SF) {
    #     min_step <- 1 + (min_step - 1) * 2
    #     max_step <- min_step * 2
    #   }
    # }

    if (SF[i+1] < stop_SF) {
      break
    }
  }
  t_PP <- t_PP[1:(i+1)]
  alpha_PP <- alpha_PP[1:(i+1)]
  SF <- SF[1:(i+1)]
  rel_error_lb <- rel_error_lb[1:i]
  rel_error_ub <- rel_error_ub[1:i]

  alpha_PP[i+1] <- max(min_tail_alpha, alpha_PP[i])
  index <- alpha_PP != c(-1, alpha_PP[1:i])
  t_PP <- t_PP[index]
  alpha_PP <- alpha_PP[index]
  SF <- SF[index]
  res <- list(t = t_PP, alpha = alpha_PP, SF = SF, original_steps=i+1, min_step = min_step, rel_error_lb = rel_error_lb, rel_error_ub = rel_error_ub)
  return(res)


}

find_t <- function(alpha) {
  nalpha <- function(x) {-alpha(x)}
  min_i <- -10
  max_i <- 20
  for (i in min_i:max_i) {
    res <- optim(10^i / 2, nalpha, lower = 0, upper = 10^i, method = "Brent")
    if (res$value < 0) {
      break
    }
  }
  if (i == min_i) {
    t <- 0
  } else if (i == max_i & res$value < 0) {
    t <- Inf
  } else {
    ub <- 10^i
    lb <- 10^(i-1)
    while (ub/lb > 1.000001) {
      res <- optim((lb + ub) / 2, nalpha, lower = lb, upper = (lb + ub) / 2, method = "Brent")
      if (res$value < 0) {
        ub <- (ub + lb) / 2
      } else {
        lb <- (ub + lb) / 2
      }
    }
    t <- lb
  }
  return(t)
}






#' Simulation of the local Pareto distribution
#'
#' @description Generates random deviates of a local Pareto distribution
#'
#' @param n Numeric. Number of simulations
#' @param alpha Function. \code{alpha(x))} is the local Pareto alpha at \code{x}.
#' @param t Numeric. Threshold of the local Pareto distribution \code{alpha(x)} is used for \code{x>t}. If NULL then \code{alpha}
#'          is used on \code{[0,infty)}
#'
#' @return A vector of \code{n} samples from the local Pareto distribution with threshold \code{t} and parameter function \code{alpha}
#'
#' @examples
#' t <- 1000
#' alpha <- function(x) {2-1000/x}
#' rLocalPareto(100, alpha, t)
#'
#' @export


rLocalPareto <- function(n, alpha, t = NULL) {

  Parameters_PP <- LocalPareto_2_PiecewisePareto(alpha, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }
  return(rPP(n, Parameters_PP))
}








#' CDF of the local Pareto Distribution
#'
#' @description Calculates the CDF of a local Pareto Distribution
#'
#' @param alpha Function. \code{alpha(x)} is the local Pareto alpha at \code{x}.
#' @param t Numeric. Threshold of the local Pareto distribution \code{alpha(x)} is used for \code{x>t}. If NULL then \code{alpha}
#'          is used on \code{[0,infty)}
#'
#' @return CDF of the piecewise Pareto distribution with threshold \code{t} and parameter function \code{alpha} evaluated at \code{x}
#'
#' @examples
#' t <- 1000
#' alpha <- function(x) {2-1000/x}
#' x <- 0:10 * 1000
#' pLocalPareto(x, alpha, t)
#'
#' @export

pLocalPareto <- function(x, alpha, t = NULL) {
  if (is.null(x) || (is.atomic(x) && length(x) == 0)) {
    return(numeric())
  }

  Parameters_PP <- LocalPareto_2_PiecewisePareto(alpha, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }

  return(pPP(x, Parameters_PP))
}



#' PDF of the local Pareto Distribution
#'
#' @description Calculates the PDF of a local Pareto Distribution
#'
#' @param alpha Function. \code{alpha(x)} is the local Pareto alpha at \code{x}.
#' @param t Numeric. Threshold of the local Pareto distribution \code{alpha(x)} is used for \code{x>t}. If NULL then \code{alpha}
#'          is used on \code{[0,infty)}
#'
#' @return PDF of the piecewise Pareto distribution with threshold \code{t} and parameter function \code{alpha} evaluated at \code{x}
#'
#' @examples
#' t <- 1000
#' alpha <- function(x) {2-1000/x}
#' x <- 0:10 * 1000
#' dLocalPareto(x, alpha, t)
#'
#' @export

dLocalPareto <- function(x, alpha, t = NULL) {
  if (is.null(x) || (is.atomic(x) && length(x) == 0)) {
    return(numeric())
  }

  Parameters_PP <- LocalPareto_2_PiecewisePareto(alpha, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }

  return(dPP(x, Parameters_PP))

}








#' Maximum Likelihood Estimation for the local Pareto distribution
#'
#' @description Calculates the maximum likelihood estimator for the alpha_0 if alpha(x) is of the form alpha(x) = alpha_0 * beta(x)
#' with a known beta
#'
#' @param losses Numeric vector. Losses that are used for the ML estimation.
#' @param t Numeric. Threshold of the local Pareto distribution.
#' @param beta Function. Function of local Pareto alphas.
#' @param weights Numeric vector or NULL. Weights for the losses. For instance \code{weights[i] = 2} and \code{weights[j] = 1} for \code{j != i} has the same effect as adding another loss of size \code{loss[i]}.
#'
#' @return Maximum likelihood estimator for the parameters \code{alpha_0} of a local Pareto distribution with alpha(x) = alpha_0 * beta(x)
#'
#' @examples
#' t <- 1000
#' beta <- function(x) {2-1000/x}
#' alpha_0 <- 2
#' alpha <- function(x) {alpha_0 * beta(x)}
#' losses <- rLocalPareto(10000, alpha, t)
#' LocalPareto_ML_Estimator_alpha_0(losses, t, beta)
#'
#' @export

LocalPareto_ML_Estimator_alpha_0 <- function(losses, t, beta, weights = NULL) {
  if (!is.nonnegative.finite.vector(losses)) {
    warning("losses must be non-negative.")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(t)) {
    warning("t must be a non-negative number.")
    return(NaN)
  }
  Parameters_PP <- LocalPareto_2_PiecewisePareto(beta, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }

  if (is.null(weights)) weights <- rep(1, length(losses))
  if (!is.positive.finite.vector(weights)) {
    warning("weights must NULL or positive.")
    return(NaN)
  }
  if (length(weights) != length(losses)) {
    warning("weights must have the same length as losses.")
    return(NaN)
  }

  t_PP <- Parameters_PP$t

  index <- (losses > t & losses > t_PP[1])
  if (sum(index) == 0) {
    warning("No loss is larger than t and the specific reporting_threshold.")
    return(NaN)
  }
  losses <- losses[index]
  weights <- weights[index]

  n <- length(losses)

  SF <- 1 - pPP(losses, Parameters_PP)
  alpha_0_hat <- n / sum(-log(SF))

  alpha_hat <- function(x) {ifelse(x>=t, alpha_0_hat * beta(x), 0)}
  MLE <- list(alpha_0_hat = alpha_0_hat, FisherInformation = n / alpha_0_hat^2, ApproxVariance = alpha_0_hat^2 / n, alpha_hat = alpha_hat)

  return(MLE)

}





#' Layer mean of the local Pareto Distribution

#' @description  Calculates the expected loss of a local Pareto
#'               distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha Function. \code{alpha(x)} is the local Pareto alpha at \code{x}.
#' @param t Numeric. Threshold of the local Pareto distribution \code{alpha(x)} is used for \code{x>t}. If NULL then \code{alpha}
#'          is used on \code{[0,infty)}
#'
#' @return The expected loss of the local Pareto distribution with threshold \code{t} and parameter function \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' t <- 1000
#' alpha <- function(x) {2-1000/x}
#' Cover <- c(1000, Inf)
#' AttachmentPoint <- c(1000, 2000)
#' LocalPareto_Layer_Mean(Cover, AttachmentPoint, alpha, t = 1000)
#'
#' @export

LocalPareto_Layer_Mean <- function(Cover, AttachmentPoint, alpha, t = NULL) {
  if (!is.nonnegative.vector(Cover)) {
    warning("Cover must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.nonnegative.finite.vector(AttachmentPoint)) {
    warning("AttachmentPoint must be a nonnegative vector.")
    return(NaN)
  }

  Parameters_PP <- LocalPareto_2_PiecewisePareto(alpha, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }

  return(PP_Layer_Mean(Cover, AttachmentPoint, Parameters_PP))

}



#' Layer variance of the local Pareto Distribution

#' @description  Calculates the variance of a local Pareto
#'               distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha Function. \code{alpha(x)} is the local Pareto alpha at \code{x}.
#' @param t Numeric. Threshold of the local Pareto distribution \code{alpha(x)} is used for \code{x>t}. If NULL then \code{alpha}
#'          is used on \code{[0,infty)}
#'
#' @return The variance of the local Pareto distribution with threshold \code{t} and parameter function \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' t <- 1000
#' alpha <- function(x) {2.5 - 1000/x}
#' Cover <- c(1000, Inf)
#' AttachmentPoint <- c(1000, 2000)
#' LocalPareto_Layer_Var(Cover, AttachmentPoint, alpha, t = 1000)
#'
#' @export

LocalPareto_Layer_Var <- function(Cover, AttachmentPoint, alpha, t = NULL) {
  if (!is.nonnegative.vector(Cover)) {
    warning("Cover must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.nonnegative.finite.vector(AttachmentPoint)) {
    warning("AttachmentPoint must be a nonnegative vector.")
    return(NaN)
  }

  Parameters_PP <- LocalPareto_2_PiecewisePareto(alpha, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }

  return(PP_Layer_Var(Cover, AttachmentPoint, Parameters_PP))

}


#' Second layer moment of the local Pareto Distribution

#' @description  Calculates the second moment of a local Pareto
#'               distribution in a reinsurance layer
#'
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#' @param alpha Function. \code{alpha(x)} is the local Pareto alpha at \code{x}.
#' @param t Numeric. Threshold of the local Pareto distribution \code{alpha(x)} is used for \code{x>t}. If NULL then \code{alpha}
#'          is used on \code{[0,infty)}
#'
#' @return The second moment of the local Pareto distribution with threshold \code{t} and parameter function \code{alpha} in the layer
#'         \code{Cover} xs \code{AttachmentPoint}
#'
#' @examples
#' t <- 1000
#' alpha <- function(x) {2.5 - 1000/x}
#' Cover <- c(1000, Inf)
#' AttachmentPoint <- c(1000, 2000)
#' LocalPareto_Layer_SM(Cover, AttachmentPoint, alpha, t = 1000)
#'
#' @export

LocalPareto_Layer_SM <- function(Cover, AttachmentPoint, alpha, t = NULL) {
  if (!is.nonnegative.vector(Cover)) {
    warning("Cover must be a nonnegative vector.")
    return(NaN)
  }
  if (!is.nonnegative.finite.vector(AttachmentPoint)) {
    warning("AttachmentPoint must be a nonnegative vector.")
    return(NaN)
  }

  Parameters_PP <- LocalPareto_2_PiecewisePareto(alpha, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }

  return(PP_Layer_SM(Cover, AttachmentPoint, Parameters_PP))

}





#' Approximate local Pareto distribution by a piecewise Pareto distribution
#'
#' @description  Calculates the parameters of an approximating piecewise Pareto distribution
#'
#' @param alpha Function.
#' @param t Numeric. Threshold of the distribution. If \code{t = 0} then \code{t} is derived from \code{alpha}.
#' @param rel.tolerance Numeric. Maximum of the relative error of the approximation for the survival function.
#' @param stop.EP Numeric. Algorithm stops if the survival function is below \code{stop.EP}.
#'
#' @return List with parameters and infos.
#'
#' @examples
#' alpha_0 <- 1.5
#' t <- 1000
#' gamma <- 0.5
#' alpha <- function(x) alpha_0 + alpha_0 * gamma * log(x/t)
#' LocalPareto_2_PiecewisePareto(t, alpha)
#'
#' @export


LocalPareto_2_PiecewisePareto <- function(alpha, t = NULL, rel.tolerance = 1e-4, stop.EP = 1e-6) {

  # check if alpha is a vectorized function
  res <- NULL
  try(res <- alpha(c(0,1)), silent = T)
  if (!is.numeric(res)) {
    alpha_s <- alpha
    alpha <- Vectorize(alpha_s)
    res <- NULL
    try(res <- alpha(c(0,1)), silent = T)
    if (!is.numeric(res)) {
      return("alpha is not a vectorized function")
    }
  }
  if (is.positive.finite.number(t)) {
    x_test <- t * 1.9^(0:20)
  } else {
    x_test <- 1.9^((-20):20)
  }
  if (min(alpha(x_test)) < 0) {
    return("alpha is not non-negative")
  }

  if (is.null(t)) {
    t <- find_t(alpha)
  } else{
    if (!is.nonnegative.finite.number(t)) {
      return("t must be a non-negative number")
    }
  }

  epsilon <- rel.tolerance #min(rel.tolerance, abs.tolerance)
  phi_max <- 2
  min_tail_alpha <- 1.01


  t_PP <- numeric(5000)
  alpha_PP <- numeric(5000)
  SF <- numeric(5000)
  xi <- numeric(5000)
  zeta <- numeric(5000)
  rel_error_lb  <- numeric(5000)
  rel_error_ub  <- numeric(5000)

  nalpha <- function(x) {-alpha(x)}
  alpha_by_x <- function(x) {alpha(x)/x}

  if (t == 0) {
    min_i <- -33
    max_i <- 3
    for (i in min_i:max_i) {
      res <- NULL
      try(res <- integrate(alpha_by_x, lower = 0, upper = 2^i), silent = T)
      if (!is.numeric(res$value)) {
        return("alpha not admissible")
      }
      if (res$value > -log(0.99)) {
        break
      }
    }
    t_1 <- 2^(i-1)
    SF_t_1 <- exp(- integrate(alpha_by_x, lower = 0, upper = t_1)$value)
    rate <- -log(SF_t_1) / t_1

    t_0 <- t_1 / 1000
    n <- 30
    t_PP[1:(n+1)] <- exp(log(1/1000) / n * n:0) * t_1
    SF[1] <- exp(- integrate(alpha_by_x, lower = 0, upper = t_PP[1])$value)
    for (i in 2:(n+1)) {
      SF[i] <- exp(- integrate(alpha_by_x, lower = t_PP[i-1], upper = t_PP[i])$value) * SF[i-1]
    }
    alpha_PP[1] <- log(1 / SF[2]) / log(t_PP[2] / t_PP[1])
    for (i in 2:n) {
      alpha_PP[i] <- log(SF[i] / SF[i+1]) / log(t_PP[i+1] / t_PP[i])
    }
  } else {
    n <- 0
    t_PP[1] <- t
    SF[1] <- 1
  }



  phi <- phi_max
  for (i in (n+1):4999) {
    #max_step_temp <- max_step
    phi <- min(phi_max, 1 + (phi - 1) * 1.5)
    max_alpha <- - stats::optim(t_PP[i], nalpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * phi)$value
    min_alpha <- stats::optim(t_PP[i], alpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * phi)$value
    while (log(1 + epsilon) / log(phi) < max_alpha - min_alpha) {
      phi <- 1 + (phi - 1) / 1.5  #max(min_step,  1 + (max_step_temp - 1) / 2)
      max_alpha <- - stats::optim(t_PP[i], nalpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * phi)$value
      min_alpha <- stats::optim(t_PP[i], alpha, method = "Brent", lower = t_PP[i], upper = t_PP[i] * phi)$value
    }
    if (max_alpha == min_alpha) {
      t_PP[i+1] <- t_PP[i] * phi
      alpha_PP[i] <- max_alpha
      SF[i+1] <- SF[i] * phi^(-max_alpha)

      xi[i] <- (t_PP[i] + t_PP[i+1]) / 2
      zeta[i] <- (t_PP[i] + t_PP[i+1]) / 2
      rel_error_lb[i] <- 0
      rel_error_ub[i] <- 0
    } else {   # then max_alpha > 0
      #step <- min(target_prob ^ (- 1 / max_alpha), max_step_temp)
      t_PP[i+1] <- t_PP[i] * phi
      factor_SF <- NULL
      try(factor_SF <- exp(-integrate(alpha_by_x, lower = t_PP[i], upper = t_PP[i+1])$value), silent = T)
      if (is.null(factor_SF)) {
        alpha_PP[i] <- alpha_PP[i-1]
        SF[i+1] <- SF[i] * (t_PP[i] / t_PP[i+1])^alpha_PP[i]
        warning("stop criterium not reached!")
        break
      }
      SF[i+1] <- SF[i] * factor_SF
      alpha_PP[i] <- log(SF[i+1] / SF[i]) / log(t_PP[i] / t_PP[i+1])

      xi[i] <- t_PP[i] ^ ((alpha_PP[i] - max_alpha) / (min_alpha - max_alpha)) * t_PP[i+1] ^ ((min_alpha - alpha_PP[i]) / (min_alpha - max_alpha))
      zeta[i] <- t_PP[i] ^ ((alpha_PP[i] - min_alpha) / (max_alpha - min_alpha)) * t_PP[i+1] ^ ((max_alpha - alpha_PP[i]) / (max_alpha - min_alpha))
      rel_error_lb[i] <- (t_PP[i] / zeta[i]) ^ (alpha_PP[i] - min_alpha) - 1
      rel_error_ub[i] <- (xi[i] / t_PP[i]) ^ (max_alpha - alpha_PP[i]) - 1
    }

    if (SF[i+1] < stop.EP) {
      break
    }
  }
  t_PP <- t_PP[1:(i+1)]
  alpha_PP <- alpha_PP[1:(i+1)]
  SF <- SF[1:(i+1)]
  rel_error_lb <- rel_error_lb[1:i]
  rel_error_ub <- rel_error_ub[1:i]

  ##################

  factor_SF <- NULL
  try(factor_SF <- exp(-integrate(alpha_by_x, lower = t_PP[i+1], upper = t_PP[i+1] * phi_max)$value), silent = T)
  if (is.null(factor_SF)) {
    alpha_PP[i+1] <- max(min_tail_alpha, alpha_PP[i])
  } else {
    alpha_PP[i+1] <- - log(factor_SF) / log(phi_max)
  }
  ##############

  #alpha_PP[i+1] <- max(min_tail_alpha, alpha_PP[i])

  index <- alpha_PP != c(-1, alpha_PP[1:i])
  t_PP <- t_PP[index]
  alpha_PP <- alpha_PP[index]
  SF <- SF[index]
  res <- list(t = t_PP, alpha = alpha_PP, SF = SF, original_steps=i+1, epsilon = epsilon, rel_error_lb = rel_error_lb, rel_error_ub = rel_error_ub)
  return(res)


}

