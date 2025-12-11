#' Expected Loss of a Reinsurance Layer
#'
#' @description  Calculates the expected loss of a reinsurance layer for a collective model
#'
#' @param CollectiveModel A collective model object. Currently only \code{PLP_Models} are handled.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The expected loss of the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @examples
#' t <- 1000
#' alpha <- function(x) 2 - 1000 / x
#' PLPM <- PLP_Model(2, alpha, t, dispersion = 2)
#' Layer_Mean(PLPM, 2000, 2000)
#'
#' @export

Layer_Mean <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) UseMethod("Layer_Mean")



#' Expected Loss of a Reinsurance Layer
#'
#' @description  Calculates the expected loss of a reinsurance layer for a PPP_Model
#'
#' @param CollectiveModel PPP_Model object.
#' @param Cover Numeric. Cover of the reinsurance layer. Use \code{Inf} for unlimited layers.
#' @param AttachmentPoint Numeric. Attachment point of the reinsurance layer.
#'
#' @return The expected loss of the layer \code{Cover} xs \code{AttachmentPoint} for the given \code{CollectiveModel}
#'
#' @export

Layer_Mean.PLP_Model <- function(CollectiveModel, Cover = Inf, AttachmentPoint = 0) {
  PP_Layer_Mean(Cover, AttachmentPoint, CollectiveModel$Parameters_PP) * CollectiveModel$FQ
}


