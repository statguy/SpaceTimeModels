#' @export findScale
#' @keywords internal
findScale <- function(x) max(10^floor(log10(abs(x))), 1)

#' @export nullScale
#' @keywords internal
nullScale <- function(x, y) {
  if (is.null(x)) return(NULL)
  return(x / y)
}

#' @export getINLAModelMatrix
#' @keywords internal
getINLAModelMatrix = function(covariatesModel, covariates) {
  if (missing(covariatesModel) || is.null(covariatesModel))
    stop("Required argument 'covariatesModel' missing.")
  
  x <- if (missing(covariates) || is.null(covariates)) terms(covariatesModel)
  else terms(covariatesModel, data=covariates)
  
  if (length(attr(x, "term.labels")) > 0) {
    if (missing(covariates) || is.null(covariates))
      stop("Covariates data do not match with covariates model.")
    
    modelMatrix <- as.data.frame(model.matrix(covariatesModel, data=covariates))
    terms <- colnames(modelMatrix)
    interceptIndex <- terms %in% "(Intercept)"
    if (any(interceptIndex)) {
      terms <- terms[!interceptIndex]
      modelMatrix <- modelMatrix[,!interceptIndex]
    }
    
    if (any(is.na(modelMatrix)))
      stop("Covariates contain missing values which are not allowed.")
    
    return(modelMatrix)
  }
  else return(NULL)
}

#' @export summaryINLAParameter
#' @keywords internal
summaryINLAParameter <- function(marginal, fun=identity, coordsScale=1) {
  m <- inla.tmarginal(function(x) fun(x) * coordsScale, marginal)
  e <- inla.emarginal(function(x) x, m)
  e2 <- inla.emarginal(function(x) x^2, m)
  sd <- sqrt(e2-e^2)
  q <- inla.qmarginal(c(0.025, 0.5, 0.975), m)
  mode <- inla.mmarginal(m)
  x <- data.frame(e=e, sd=sd, q1=q[1], q2=q[2], q3=q[3], mode=mode)
  colnames(x) <- c("mean", "sd", "0.025quant","0.5quant","0.975quant", "mode")
  return(x)
}
