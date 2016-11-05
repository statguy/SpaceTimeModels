#' @export findScale
#' @keywords internal
findScale <- function(x) max(10^floor(log10(abs(x))), 1)

#' @export nullScale
#' @keywords internal
nullScale <- function(x, y) {
  if (is.null(x)) return(NULL)
  return(x / y)
}

#' @export getCovariateNames
#' @keywords internal
getCovariateNames = function(covariatesModel, covariates) {
  x <- if (missing(covariates) || is.null(covariates)) terms(covariatesModel)
  else terms(covariatesModel, data = covariates)
  return(attr(x, "term.labels"))
}

#' @export getINLAModelMatrix
#' @keywords internal
getINLAModelMatrix = function(covariatesModel, covariates) {
  if (missing(covariatesModel) || is.null(covariatesModel))
    stop("Required argument 'covariatesModel' missing.")
  
  covariateNames <- SpaceTimeModels::getCovariateNames(covariatesModel, covariates)
  if (length(covariateNames) > 0) {
    if (missing(covariates) || is.null(covariates))
      stop("Covariates data do not match with covariates model.")
    
    modelMatrix <- as.data.frame(model.matrix(covariatesModel, data = covariates))
    terms <- colnames(modelMatrix)
    interceptIndex <- terms %in% "(Intercept)"
    if (any(interceptIndex)) {
      terms <- terms[!interceptIndex]
      modelMatrix <- modelMatrix[,!interceptIndex, drop = F]
    }
    
    if (any(is.na(modelMatrix)))
      stop("Covariates contain missing values which are not allowed.")
    
    return(modelMatrix)
  }
  else return(NULL)
}

#' @export summaryINLAParameter
#' @keywords internal
summaryINLAParameter <- function(marginal, fun = identity, coordinatesScale = 1) {
  m <- inla.tmarginal(function(x) fun(x) * coordinatesScale, marginal)
  e <- inla.emarginal(function(x) x, m)
  e2 <- inla.emarginal(function(x) x^2, m)
  sd <- sqrt(e2 - e^2)
  q <- inla.qmarginal(c(0.025, 0.5, 0.975), m)
  mode <- inla.mmarginal(m)
  x <- data.frame(e = e, sd = sd, q1 = q[1], q2 = q[2], q3 = q[3], mode = mode)
  colnames(x) <- c("mean", "sd", "0.025quant","0.5quant","0.975quant", "mode")
  return(x)
}

#' @export assertCompleteCovariates
#' @keywords internal
assertCompleteCovariates <- function(covariatesModel, covariates) {
  #if (inherits(covariates, "STIDF") == FALSE)
  #  stop("Covariates must be of class STIDF or a subclass.")
  x <- terms(covariatesModel, data = covariates)
  #complete <- complete.cases(covariates@data[,attr(x, "term.labels"), drop = F])
  complete <- complete.cases(covariates[,attr(x, "term.labels"), drop = F])
  if (any(complete == FALSE))
    stop("Covariates cannot contain missing data.")
}

#' @export theme_raster
theme_raster <- function(base_size = 12, base_family = "", ...) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.border = element_blank(),
      panel.margin = unit(0, "lines"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      #plot.margin = unit(c(0, 0, -1, -1), "lines"),
      #plot.margin = unit(c(0, 0, -.5, -.5), "lines"),
      plot.margin = unit(c(0, 0, 0, 0), "lines"),
      #axis.ticks.length = unit(0, "lines"), axis.ticks.margin = unit(0, "lines"),
      #plot.margin = rep(unit(0, "null"), 4), panel.margin = unit(0, "null"), axis.ticks.length = unit(0, "null"), axis.ticks.margin = unit(0, "null"),
      legend.position = "none",
      ...
    )
}

# Taken from https://groups.google.com/forum/#!topic/r-inla-discussion-group/cPU0iJA2UqY
#' @export local.inla.spde2.matern.new
#' @keywords internal
# - YOU MUST SET rho0 , e.g. to half the length/width of your space
# - You may reduce sig0 to have less spatial effect
# - You may increase rho0 additionally, to make the spatial effect smoother
local.inla.spde2.matern.new <- function(mesh, alpha=2, prior.pc.rho, prior.pc.sig) {
  # Call inla.spde2.matern with range and standard deviation parametrization
  d = INLA:::inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
  nu = alpha-d/2
  kappa0 = log(8*nu)/2
  tau0   = 0.5*(lgamma(nu)-lgamma(nu+d/2)-d/2*log(4*pi))-nu*kappa0
  spde   = inla.spde2.matern(mesh = mesh,
                             B.tau   = cbind(tau0,   nu,  -1),
                             B.kappa = cbind(kappa0, -1, 0))
  
  # Change prior information
  param = c(prior.pc.rho, prior.pc.sig)
  spde$f$hyper.default$theta1$prior = "pcspdega"
  spde$f$hyper.default$theta1$param = param
  spde$f$hyper.default$theta1$initial = log(prior.pc.rho[1])+1
  spde$f$hyper.default$theta2$initial = log(prior.pc.sig[1])-1
  
  # End and return
  return(invisible(spde))  
}
