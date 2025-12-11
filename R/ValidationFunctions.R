is.vector <- function(x) {
  if(!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    return(FALSE)
  }
  return(TRUE)
}

is.positive.vector <- function(x) {
  if(!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    return(FALSE)
  }
  k <- length(x)
  if (sum(x > 0, na.rm = TRUE) < k) {
    return(FALSE)
  }
  return(TRUE)
}


is.positive_or_NA.finite.vector <- function(x) {
  if(!is.atomic(x) || length(x) < 1) {
    return(FALSE)
  }
  if (length(x[!is.na(x)]) > 0 && !is.positive.finite.vector(x[!is.na(x)])) {
    return(FALSE)
  }
  return(TRUE)
}


is.nonnegative_or_NA.finite.vector <- function(x) {
  if(!is.atomic(x) || length(x) < 1) {
    return(FALSE)
  }
  if (length(x[!is.na(x)]) > 0 && !is.nonnegative.finite.vector(x[!is.na(x)])) {
    return(FALSE)
  }
  return(TRUE)
}

is.NA.vector <- function(x) {
  if(!is.atomic(x) || length(x) < 1) {
    return(FALSE)
  }
  if (length(x[!is.na(x)]) > 0) {
    return(FALSE)
  }
  return(TRUE)
}



is.nonnegative.vector <- function(x) {
  if(!is.atomic(x) || !is.numeric(x) || length(x) < 1) {
    return(FALSE)
  }
  k <- length(x)
  if (sum(x >= 0, na.rm = TRUE) < k) {
    return(FALSE)
  }
  return(TRUE)
}


is.positive.finite.vector <- function(x) {
  if(!is.positive.vector(x)) {
    return(FALSE)
  }
  if (sum(is.infinite(x)) > 0) {
    return(FALSE)
  }
  return(TRUE)
}


is.nonnegative.finite.vector <- function(x) {
  if(!is.nonnegative.vector(x)) {
    return(FALSE)
  }
  if (sum(is.infinite(x)) > 0) {
    return(FALSE)
  }
  return(TRUE)
}

is.number <- function(x) {
  if(!is.atomic(x) || !is.numeric(x) || length(x) != 1 || is.na(x)) {
    return(FALSE)
  }
  return(TRUE)
}

is.positive.number <- function(x) {
  if (!is.positive.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.positive.finite.number <- function(x) {
  if (!is.positive.finite.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.nonnegative.number <- function(x) {
  if (!is.nonnegative.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.nonnegative.finite.number <- function(x) {
  if (!is.nonnegative.finite.vector(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.string <- function(x) {
  if (!is.character(x) || !is.atomic(x) || length(x) != 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

is.TRUEorFALSE <- function(x) {
  if (isTRUE(x) || isFALSE(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


is.TRUEorFALSE.vector <- function(x) {
  if(!is.atomic(x) || !is.logical(x) || length(x) < 1) {
    return(FALSE)
  }
  k <- length(x)
  if (sum(x, na.rm = T) + sum(!x, na.rm = T) < k) {
    return(FALSE)
  }
  return(TRUE)
}






valid.parameters.LALocPareto <- function(t, alpha_0, gamma, delta) {

  if (!is.positive.finite.number(t) || !is.positive.finite.number(alpha_0)) {
    return("t and alpha must be positive numbers.")
    return(FALSE)
  }
  if (!is.positive.finite.number(gamma) && !is.positive.finite.number(delta)) {
    return("Either gamma or delta must be a positive number.")
    return(FALSE)
  }
  if (is.positive.finite.number(gamma) && is.positive.finite.number(delta)) {
    if (round(gamma * alpha_0 * log(2), 8) != round(delta, 8)) {
      return("gamma * alpha_0 * log(2) must equal delta.")
      return(FALSE)
    }
  }
  return(TRUE)
}


