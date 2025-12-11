#' PLP_Model (Collective Panjer & Local Pareto Model) Object
#'
#' @description Constructor function for the PLP_Model object
#'
#' @param FQ Numerical. Expected claim count of the collective model.
#' @param alpha Function. \code{alpha(x)} is the local Pareto alpha at \code{x}.
#' @param t Numerical. Optional thresholds of the local Pareto distribution. If \code{t} is not \code{NULL} then alpha is only used above \code{t}.
#' @param dispersion Numerical. Dispersion of the Panjer distribution (i.e. variance to mean ratio).
#' @param Status Numerical indicator if a function returns a PPP_Model object: 0 = success, 1 = some information has been ignored, 2 = no solution found
#' @param Comment Charakter. An optional comment.

#' @examples
#' t <- 1000
#' alpha <- function(x) {2-1000/x}
#' PLPM <- PLP_Model(2, alpha, t, dispersion = 2)
#' PLPM
#'
#' @export

PLP_Model <- function(FQ = NULL, alpha, t = NULL, dispersion = 1, Status = 0, Comment = "OK") {

  Parameters_PP <- LocalPareto_2_PiecewisePareto(alpha, t)
  if (is.character(Parameters_PP)) {
    stop(Parameters_PP)
  }

  obj <- list(FQ = FQ, alpha = alpha, t = t, Parameters_PP = Parameters_PP, dispersion = dispersion, Status = Status, Comment = Comment)
  class(obj) <- "PLP_Model"

  if (!is.valid.PLP_Model(obj)) {
    obj$Status <- 2
    obj$Comment = is.valid.PPP_Model(obj, comment = TRUE)
  }
  return(obj)

}

#' Print a PLP_Model Object (Collective Panjer & Local Pareto Model) Object
#'
#' @description Print method for PLP_Model objects
#'
#' @param x PLP_Model object.
#' @param ... Other arguments, all currently ignored.
#'
#' @export

print.PLP_Model <- function(x, ...) {
  if (!is.positive.finite.number(x$dispersion)) {
    fq_dist <- "Panjer"
  } else if (x$dispersion == 1) {
    fq_dist <- "Poisson"
  } else if (x$dispersion > 1) {
    fq_dist <- "Negative Binomial"
  } else {
    fq_dist <- "Binomial"
  }
  cat("\nPanjer & Local Pareto model\n\n")
  cat("Collective model with a ", fq_dist, " distribution for the claim count and a local Pareto distributed severity.", sep = "")
  cat("\n\n", fq_dist, " Distribution:", sep = "")

  cat("\nExpected Frequency:   ", x$FQ, sep = "")
  if (is.positive.finite.number(x$dispersion)  && x$dispersion != 1) {
    cat("\nDispersion:           ", x$dispersion, sep = "")
    if (is.positive.finite.number(x$FQ) && x$dispersion > 1) {
      cat(" (i.e. contagion = ", (x$dispersion - 1)/x$FQ, ")", sep = "")
    }
  }
  ParPP <- x$Parameters_PP
  nPP <-  length(ParPP$t)
  cat("\n\nLocal Pareto Distribution:")
  cat("\nalpha:           ", paste0("alpha <- function(",names(formals(alpha)),") "), deparse(body(alpha)), sep = " ")
  cat("\nApproximated with", paste0(nPP), "Pareto pieces.")

  cat("\n\nStatus:              ", x$Status)
  cat("\nComments:            ", x$Comment)
  if (!is.valid.PLP_Model(x)) {
    cat("\n\nThe model is not valid.\n")
    cat(is.valid.PLP_Model(x, comment = TRUE))
  }
  cat("\n\n")

}

#' Check if an object is a PLP_Model
#'
#' @description Checks if the class of an object is 'PLP_Model'
#'
#' @param x Object to be checked.

#' @examples
#' t <- 1000
#' alpha <- function(x) {2-1000/x}
#' PLPM <- PLP_Model(2, alpha, t, dispersion = 2)
#' PLPM
#' is.PLP_Model(PLPM)
#'
#' @export

is.PLP_Model <- function(x) {
  if (inherits(x, "PLP_Model")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


#' Check if an object is a valid PPP_Model
#'
#' @description Checks if an object is a PPP_Model object and whether it is valid for the use in functions like \code{Layer_Mean}
#'
#' @param x Object to be checked.
#' @param comment If FALSE then the function returns a boolean indicating whether \code{x} is a valid PPP_Model. If TRUE then the function returns a comment instead.

#' @examples
#' t <- 1000
#' alpha <- function(x) {2-1000/x}
#' PLPM <- PLP_Model(2, alpha, t, dispersion = 2)
#' PLPM
#' is.valid.PLP_Model(PLPM)
#' is.valid.PLP_Model(PLPM, comment = TRUE)
#'
#' @export

is.valid.PLP_Model <- function(x, comment = FALSE) {
  if (!inherits(x, "PLP_Model") || typeof(x) != "list") {
    if (!comment) {
      return(FALSE)
    } else {
      return("Object does not have class PLP_Model.")
    }
  }
  required_elements <- c("FQ", "alpha", "t", "Parameters_PP", "dispersion", "Status", "Comment")
  available <- required_elements %in% names(x)
  if (sum(!available) > 0) {
    if (!comment) {
      return(FALSE)
    } else {
      return(paste("Not all required list elements available. Missing elements:", paste(required_elements[!available], collapse = ", ")))
    }
  }


  if (!is.nonnegative.finite.number(x$FQ)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("FQ must be a nonnegative number.")
    }
  }
  # if (!valid.parameters.PiecewisePareto()) {
  #   if (!comment) {
  #     return(FALSE)
  #   } else {
  #     return(valid.parameters.PiecewisePareto(x$t, x$alpha, x$truncation, x$truncation_type, comment = TRUE))
  #   }
  # }
  if (!is.positive.finite.number(x$dispersion)) {
    if (!comment) {
      return(FALSE)
    } else {
      return("dispersion must be a positive number.")
    }
  }


  if (!comment) {
    return(TRUE)
  } else {
    return("OK")
  }
}


