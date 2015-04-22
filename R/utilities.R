#' @export findScale
#' @keywords internal
findScale <- function(x) max(10^floor(log10(abs(x))), 1)

#' @export nullScale
#' @keywords internal
nullScale <- function(x, y) {
  if (is.null(x)) return(NULL)
  return(x / y)
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
