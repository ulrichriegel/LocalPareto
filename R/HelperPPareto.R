

pPP <- function(x, Parameters_PP) {
  t_PP <- Parameters_PP$t
  alpha_PP <- Parameters_PP$alpha
  excess_prob <- Parameters_PP$SF

  k <- length(t_PP)

  index1 <- rep(F, length(x))
  index2 <- index1
  EP <- rep(NaN, length(x))

  for (i in 1:k) {
    index1 <- x <= t_PP[i]
    diff_index <- index1 & !index2
    if (i == 1) {
      EP[diff_index] <- 1
    } else {
      EP[diff_index] <- (t_PP[i-1] / x[diff_index])^alpha_PP[i-1] * excess_prob[i-1]
    }
    index2 <- index1
  }
  EP[!index1] <- (t_PP[k] / x[!index1])^alpha_PP[k] * excess_prob[k]

  return(1 - EP)

}


dPP <- function(x, Parameters_PP) {
  t_PP <- Parameters_PP$t
  alpha_PP <- Parameters_PP$alpha
  excess_prob <- Parameters_PP$SF

  k <- length(t_PP)
  index1 <- rep(F, length(x))
  index2 <- index1
  Dens <- rep(NaN, length(x))

  for (i in 1:k) {
    index1 <- x < t_PP[i]
    diff_index <- index1 & !index2
    if (i == 1) {
      Dens[diff_index] <- 0
    } else {
      Dens[diff_index] <- alpha_PP[i-1] * t_PP[i-1]^alpha_PP[i-1] * x[diff_index]^(-alpha_PP[i-1]-1) * excess_prob[i-1]
    }
    index2 <- index1
  }
  Dens[!index1] <- alpha_PP[k] * t_PP[k]^alpha_PP[k] * x[!index1]^(-alpha_PP[k]-1) * excess_prob[k]

  return(Dens)
}



qPP <- function(y, Parameters_PP) {
  t_PP <- Parameters_PP$t
  alpha_PP <- Parameters_PP$alpha
  excess_prob <- Parameters_PP$SF

  k <- length(t_PP)
  CDF_at_t <- 1 - excess_prob # pPP(t_PP, Parameters_PP)

  index1 <- rep(F, length(y))
  index2 <- index1
  Amount <- rep(NaN, length(y))

  for (i in 1:k) {
    index1 <- y <= CDF_at_t[i]
    diff_index <- index1 & !index2
    if (i == 1) {
      Amount[diff_index] <- t_PP[1]
    } else {
      Amount[diff_index] <- t_PP[i-1] / ((1 - y[diff_index]) / (1 - CDF_at_t[i-1]))^(1 / alpha_PP[i-1])
    }
    index2 <- index1
  }
  Amount[!index1] <- t_PP[k] / ((1 - y[!index1]) / (1 - CDF_at_t[k]))^(1 / alpha_PP[k])

  return(Amount)
}






PP_Layer_Mean <- function(Cover, AttachmentPoint, Parameters_PP) {
  VecFun <- Vectorize(PP_Layer_Mean_s, c("Cover", "AttachmentPoint"))
  return(VecFun(Cover, AttachmentPoint, Parameters_PP))
}

PP_Layer_Mean_s <- function(Cover, AttachmentPoint, Parameters_PP) {
  t <- Parameters_PP$t
  alpha <- Parameters_PP$alpha
  excess_prob <- Parameters_PP$SF

  if (!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }

  n <- length(t)
  if (n == 1) {
    return(P_Layer_Mean_s(Cover, AttachmentPoint, alpha, t))
  }
  if (Cover == 0) {
    return(0)
  }
  if (is.infinite(Cover) && alpha[n] <= 1) {
    return(Inf)
  }

  k1 <- sum(t <= AttachmentPoint)
  if (Cover == Inf) {
    k2 <- n
  } else {
    k2 <- sum(t < AttachmentPoint + Cover)
  }

  Cover_orig <- Cover
  if (k1 == 0 && k2 == 0) {
    return(Cover)
  } else if (k1 == 0) {
    Result <- t[1] - AttachmentPoint
    AttachmentPoint <- t[1]
    if (is.numeric(Cover)) {
      Cover <- Cover - Result
    }
    k1 <- 1
  } else {
    Result <- 0
  }

  Att <- numeric(n)
  Exit <- numeric(n)

  Att[k1:k2] <- pmax(t[k1:k2], AttachmentPoint)
  if (k2 < n) {
    Exit[k1:k2] <- pmin(t[(k1+1):(k2+1)], AttachmentPoint + Cover)
  } else if (k1 < k2) {
    if (is.numeric(Cover)) {
      Exit[k1:k2] <- c(t[(k1+1):k2], AttachmentPoint + Cover)
    } else {
      Exit[k1:k2] <- c(t[(k1+1):k2], 0)
    }
  } else {
    if (is.numeric(Cover)) {
      Exit[k1] <- AttachmentPoint + Cover
    } else {
      Exit[k1] <- 0
    }
  }

  if (!is.infinite(Cover)) {
    for (i in k1:k2) {
      Result <- Result + P_Layer_Mean_s(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
    }
  } else if (is.infinite(Cover)) {
    if (k1 < k2) {
      for (i in k1:(k2-1)) {
        Result <- Result + P_Layer_Mean_s(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
      }
    }
    Result <- Result + P_Layer_Mean_s(Inf, Att[n], alpha[n], t[n]) * excess_prob[n]
  }

  return(Result)
}




PP_Layer_SM <- function(Cover, AttachmentPoint, Parameters_PP) {
  VecFun <- Vectorize(PP_Layer_SM_s, c("Cover", "AttachmentPoint"))
  return(VecFun(Cover, AttachmentPoint, Parameters_PP))
}


PP_Layer_SM_s <- function(Cover, AttachmentPoint, Parameters_PP) {
  t <- Parameters_PP$t
  alpha <- Parameters_PP$alpha
  excess_prob <- Parameters_PP$SF

  n <- length(t)

  if (!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }


  if (n == 1) {
    return(P_Layer_SM_s(Cover, AttachmentPoint, alpha, t))
  }

  if (Cover == 0) {
    return(0)
  }
  if (is.infinite(Cover) && alpha[n] <= 2) {
    return(Inf)
  }


  prob <- c(- diff(excess_prob), excess_prob[n])

  k1 <- sum(t <= AttachmentPoint)
  if (Cover == Inf) {
    k2 <- n
  } else {
    k2 <- sum(t < AttachmentPoint + Cover)
  }

  AttachmentPoint_orig <- AttachmentPoint
  Cover_orig <- Cover

  if (k1 == 0 && k2 == 0) {
    return(Cover^2)
  } else if (k1 == 0) {
    AttachmentPoint <- t[1]
    if (is.numeric(Cover)) {
      Cover <- Cover - (AttachmentPoint - AttachmentPoint_orig)
    }
    k1 <- 1
  }

  Result <- 0


  Att <- numeric(n)
  Exit <- numeric(n)

  Att[k1:k2] <- pmax(t[k1:k2], AttachmentPoint)
  if (k2 < n) {
    Exit[k1:k2] <- pmin(t[(k1+1):(k2+1)], AttachmentPoint + Cover)
  } else if (k1 < k2) {
    Exit[k1:k2] <- c(t[(k1+1):k2], AttachmentPoint + Cover)
  } else {
    Exit[k1] <- AttachmentPoint + Cover
  }

  for (i in k1:k2) {
    if (i == k2) {
      Result <- Result + P_Layer_SM_s(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
      Result <- Result + (Att[i] - AttachmentPoint_orig)^2 * excess_prob[i] * (t[i] / Att[i])^alpha[i]
      Result <- Result + 2 * (Att[i] - AttachmentPoint_orig) * P_Layer_Mean_s(Exit[i]-Att[i], Att[i], alpha[i], t[i]) * excess_prob[i]
    } else {
      Result <- Result + P_Layer_SM_s(Exit[i]-Att[i], Att[i], alpha[i], t[i], truncation = Exit[i]) * prob[i]
      Result <- Result + (Att[i] - AttachmentPoint_orig)^2 * prob[i] * (t[i] / Att[i])^alpha[i]
      Result <- Result + 2 * (Att[i] - AttachmentPoint_orig) * P_Layer_Mean_s(Exit[i]-Att[i], Att[i], alpha[i], t[i], truncation = Exit[i]) * prob[i]
    }
  }


  return(Result)
}






PP_Layer_Var <- function(Cover, AttachmentPoint, Parameters_PP) {
  PP_Layer_SM(Cover, AttachmentPoint, Parameters_PP) - PP_Layer_Mean(Cover, AttachmentPoint, Parameters_PP)^2
}




rPP <- function(n, Parameters_PP) {

  t_PP <- Parameters_PP$t
  alpha_PP <- Parameters_PP$alpha

  k <- length(t_PP)
  excess_prob <- Parameters_PP$SF

  prob_for_pieces <- c(-diff(excess_prob), excess_prob[k])


  Simulated_Pieces <- c(sample(1:k, n, replace = T, prob = prob_for_pieces))
  NumberOfSimulations_for_Pieces <- numeric(k)
  for (i in 1:k) {
    NumberOfSimulations_for_Pieces[i] <- sum(Simulated_Pieces == i)
  }

  FinvPareto <- function(x,t,alpha) {
    return(t/(1-x)^(1/alpha))
  }

  Result <- numeric(n)

  for (i in 1:k) {
    if (i == k) {
      CDF <- 1
    } else {
      CDF <- 1 - (t_PP[i] / t_PP[i+1])^alpha_PP[i]
    }
    Result[Simulated_Pieces == i] <- FinvPareto(stats::runif(NumberOfSimulations_for_Pieces[i], min = 0, max = CDF), t_PP[i], alpha_PP[i])
  }

  return(Result)

}



# Pareto functions





P_Layer_Mean_s <- function(Cover, AttachmentPoint, alpha, t=NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point in zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  if (is.infinite(Cover)) {
    if (alpha <= 1) {
      return(Inf)
    } else if (t <= AttachmentPoint) {
      EP <- -(t / AttachmentPoint)^alpha / (1 - alpha) * AttachmentPoint
    } else {
      EP <- t - AttachmentPoint
      EP <- EP - t / (1 - alpha)
    }
    return(EP)
  } else {
    # Calculation ignoring truncation
    if (t <= AttachmentPoint) {
      if (alpha == 0) {
        EP <- Cover
      } else if (alpha == 1) {
        EP <- t * (log(Cover + AttachmentPoint) - log(AttachmentPoint))
      } else {
        EP <- t / (1 - alpha) * (((Cover + AttachmentPoint) / t)^(1 - alpha) - (AttachmentPoint / t)^(1 - alpha))
      }
    } else if (t >= AttachmentPoint + Cover) {
      EP <- Cover
    } else {
      EP <- t - AttachmentPoint
      if (alpha == 0) {
        EP <- Cover
      } else if (alpha == 1) {
        EP <- EP + t * (log(Cover + AttachmentPoint) - log(t))
      } else {
        EP <- EP + t / (1 - alpha) * (((Cover + AttachmentPoint) / t)^(1 - alpha) - 1)
      }
    }

    if (is.positive.finite.number(truncation)) {
      # then Cover + AttachmentPoint <= truncation
      if (alpha < 1e-6) {
        IndefitineIntegral <- function(x) {
          return(x - 1 / log(truncation / t) * (x * log(x / t) - x))
        }
        if (t <= AttachmentPoint) {
          EP <- IndefitineIntegral(Cover + AttachmentPoint) - IndefitineIntegral(AttachmentPoint)
        } else if (t >= Cover + AttachmentPoint) {
          EP <- Cover
        } else {
          EP <- t - AttachmentPoint + IndefitineIntegral(Cover + AttachmentPoint) - IndefitineIntegral(t)
        }
      } else {
        # Adjustment for truncation
        FQ_at_truncation <- (t / truncation)^alpha
        EP <- (EP - FQ_at_truncation * Cover) / (1 - FQ_at_truncation)
      }
    }
    return(EP)
  }

}

P_Layer_Second_Moment_simple <- function(Cover, AttachmentPoint, alpha) {
  if (!is.positive.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a positive number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (Cover == 0) {
    return(0)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }

  if (is.infinite(Cover)) {
    if (alpha <= 2) {
      # warning("alpha must be > 2 for unlimited covers!")
      return(Inf)
    } else {
      SM <- 2 * AttachmentPoint^2 * (1/(alpha-2) - 1/(alpha-1))
    }
    return(SM)
  } else {
    if (alpha == 0) {
      SM <- Cover^2
    } else if (alpha == 1) {
      SM <- 2 * AttachmentPoint^2 * (Cover/AttachmentPoint - log(1+Cover/AttachmentPoint))
    } else if (alpha == 2) {
      SM <- 2 * AttachmentPoint^2 * (-Cover/(Cover+AttachmentPoint) + log(1+Cover/AttachmentPoint))
    } else {
      SM <- 2 * AttachmentPoint^2 * (((1+Cover/AttachmentPoint)^(2-alpha)-1) / (2-alpha) - ((1+Cover/AttachmentPoint)^(1-alpha)-1) / (1-alpha))
    }
    return(SM)
  }

}





P_Layer_Var_s <- function(Cover, AttachmentPoint, alpha, t = NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point is zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  ExitPoint <- Cover + AttachmentPoint
  if (t >= ExitPoint) {
    return(0)
  }
  AttachmentPoint <- max(AttachmentPoint, t)
  Cover <- ExitPoint - AttachmentPoint

  if (is.infinite(Cover) && alpha <= 2) {
    return(Inf)
  }

  # Second moment of layer loss if t = AttachmentPoint
  SM <- P_Layer_Second_Moment_simple(Cover, AttachmentPoint, alpha)
  if (is.positive.finite.number(truncation) && alpha > 0) {
    # probability of truncation if t = AttachmentPoint
    p <- 1- pP_s(truncation, AttachmentPoint, alpha)
    # consider truncation in second moment
    SM <- (SM - p * Cover^2) / (1 - p)
  }
  # consider thresholds t < AttachmentPoint
  p <- 1 - pP_s(AttachmentPoint, t, alpha, truncation = truncation)
  SM <- p * SM

  if (is.positive.finite.number(truncation) && alpha == 0) {
    IndefiniteIntegral <- function(x) {
      if (x <= t) {
        return(x^2 + 0.5 * t^2 / log(truncation / t))
      } else {
        return(x^2 - x^2 / log(truncation / t) * log(x / t) + 0.5 * x^2 / log(truncation / t))
      }
    }
    SM <- IndefiniteIntegral(Cover + AttachmentPoint) - IndefiniteIntegral(AttachmentPoint) - 2 * AttachmentPoint * P_Layer_Mean_s(Cover, AttachmentPoint, alpha, t, truncation = truncation)
  }


  # calculate Variance
  Result <- SM - P_Layer_Mean_s(Cover, AttachmentPoint, alpha, t, truncation = truncation)^2
  return(Result)
}



P_Layer_SM_s <- function(Cover, AttachmentPoint, alpha, t=NULL, truncation = NULL) {
  if (!is.nonnegative.finite.number(AttachmentPoint)) {
    warning("AttachmentPoint must be a non-negative number.")
    return(NaN)
  }
  if(!is.nonnegative.number(Cover)) {
    warning("Cover must be a non-negative number ('Inf' allowed).")
    return(NaN)
  }
  if (!is.nonnegative.finite.number(alpha)) {
    warning("alpha must be a non-negative number.")
    return(NaN)
  }
  if (is.null(t)) {
    if (AttachmentPoint == 0) {
      warning("If Attachment Point is zero, then a t>0 has to be entered.")
      return(NaN)
    }
    t <- AttachmentPoint
  }
  if (!is.positive.finite.number(t)) {
    warning("t must be a positive number.")
    return(NaN)
  }
  if (!is.null(truncation)) {
    if (!is.positive.number(truncation)) {
      warning("truncation must be NULL or a positive number ('Inf' allowed).")
      return(NaN)
    }
    if (truncation <= t) {
      warning("truncation must be larger than t.")
      return(NaN)
    }
    if (truncation <= AttachmentPoint) {
      return(0)
    }
    if (AttachmentPoint + Cover > truncation) {
      Cover <- truncation - AttachmentPoint
    }
  }

  Var <- P_Layer_Var_s(Cover, AttachmentPoint, alpha, t, truncation = truncation)
  if (is.infinite(Var)) {
    return(Inf)
  }
  Result <- Var + P_Layer_Mean_s(Cover, AttachmentPoint, alpha, t, truncation = truncation)^2
  return(Result)
}



pP_s <- function(x, t, alpha, truncation = NULL) {
  if (!is.number(x)) {
    warning("x must be a number ('Inf' allowed).")
    return(NaN)
  }

  if (alpha <= 1e-6) {
    if (is.positive.finite.number(truncation)) {
      if (x <= t) {
        return(0)
      } else if (x < truncation) {
        return(log(x/t) / log(truncation / t))
      } else {
        return(1)
      }
    } else {
      if (is.infinite(x) && is.positive.number(x)) {
        return(1)
      } else {
        return(0)
      }
    }

  }

  if (x <= t) {
    return(0)
  } else if (is.null(truncation)) {
    Result <- 1 - (t / x)^alpha
    return(Result)
  } else if (x >= truncation) {
    return(1)
  } else {
    Result <- (1 - (t / x)^alpha) / (1 - (t / truncation)^alpha)
    return(Result)
  }
}


