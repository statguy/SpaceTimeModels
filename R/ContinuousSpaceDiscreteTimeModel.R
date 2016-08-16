#' @title Continuous space, discrete time model
#' @description Building and estimating discrete time, continuous space models.
#' @references Lindgren, F. & Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software, 63(19).
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo <\email{jvj@@iki.fi}>
#' @exportClass ContinuousSpaceDiscreteTimeModel
#' @export ContinuousSpaceDiscreteTimeModel
ContinuousSpaceDiscreteTimeModel <- R6::R6Class(
  "ContinuousSpaceDiscreteTimeModel",
  lock_objects = FALSE,
  inherit = SpaceTimeModels::ContinuousSpaceTimeModel,
  public = list(
    getRandomEffectTerm = function() {
      return("f(spatial, model=spde, group=spatial.group, control.group=list(model=\"ar1\"))")
    },
    
    addObservationStack = function(sp, response=NA, covariates, offset, tag="obs") {
      # TODO: allow defining link function
      
      if (missing(sp))
        stop("Required argument 'sp' missing.")
      if (!inherits(sp, "STI"))
        stop("Argument 'sp' must be of class 'STI'.")
      
      if (is.null(self$getSpatialMesh()))
        stop("Mesh must be defined first.")
      if (is.null(self$getSPDEObject()))
        stop("Spatial prior must be defined first.")
      if (is.null(self$covariatesModel))
        stop("Covariates model must be defined first.")
      #if (missing(coordinates)) coordinates <- model$getSpatialMesh()$getKnots()
      
      dataList <- list(response=response)
      if (!missing(offset)) dataList$E <- offset / self$getOffsetScale()

      coordinates <- self$scaleCoordinates(sp::coordinates(sp))
      SpaceTimeModels::assertCompleteCovariates(self$covariatesModel, covariates)
      modelMatrix <- SpaceTimeModels::getINLAModelMatrix(self$covariatesModel, covariates)
      self$time <- time <- time(sp)
      self$timeIndex <- timeIndex <- time - min(time) + 1
      self$nTime <- nTime <- length(unique(timeIndex))
      fieldIndex <- inla.spde.make.index("spatial", n.spde=self$getSPDEObject()$n.spde, n.group=nTime)
      A <- inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc=coordinates, group=timeIndex, n.group=nTime)
      
      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept=1))) else list(fieldIndex)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      self$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
    
    addPredictionStack = function(sp, response=NA, covariates, tag="pred") {
      # TODO: finish this function
      
      if (missing(sp))
        stop("Required argument 'sp' missing.")
      if (!inherits(sp, "STI"))
        stop("Argument 'sp' must be of class 'STI'.")
      coordinates <- sp::coordinates(sp)
      
      effects <- if (self$hasIntercept()) list(c(fieldIndex, coordinates, list(intercept=1))) else list(c(index, coordinates))
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      self$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
    
    getSPDEResult = function() {
      if (is.null(self$result) || is.null(self$spde))
        stop("The model has not been estimated.")
      return(inla.spde2.result(self$result, "spatial", self$spde))
    },
    
    summarySpatialParameters = function() {
      spdeResult <- self$getSPDEResult()
      range <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.range.nominal[[1]], coordinatesScale=self$getCoordinatesScale())
      variance <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.variance.nominal[[1]])
      kappa <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.kappa[[1]], coordinatesScale=1/self$getCoordinatesScale())
      tau <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.tau[[1]])
      x <- rbind(kappa=kappa, tau=tau, range=range, variance=variance)
      colnames(x) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
      print(x)
      return(invisible(self))
    },
    
    summary = function() {
      if (is.null(self$result))
        stop("The model has not been estimated.")
      print(summary(self$result))
      return(invisible(self))
    },
    
    #getFitted = function() {
    #  indexObserved <- inla.stack.index(self$fullStack, "obs")$data
    #  fitted <- self$result$summary.fitted.values$mean[indexObserved] * self$offsetScale
    #  return(fitted)      
    #},
    
    getFittedResponse = function(tag="obs") {
      index <- self$getIndex(tag)
      offset <- self$getOffset(tag)
      nodes <- self$getSpatialMesh()$getNumNodes()
      data <- list()
      data$responseMean <- inla.vector2matrix(self$getResult()$summary.fitted.values$mean[index] * offset, nrow = nodes, ncol = self$nTime)
      colnames(data$responseMean) <- unique(self$time)
      data$responseSd <- inla.vector2matrix(self$getResult()$summary.fitted.values$sd[index] * offset, nrow = nodes, ncol = self$nTime)
      colnames(data$responseMean) <- unique(self$time)
      return(data)
    },
    
    getFittedLinearPredictor = function(tag="obs") {
      index <- self$getIndex(tag)
      data <- list()
      return(data) 
    },
    
    summaryTemporalVariation = function() {
      observed <- self$getObserved()
      fitted <- self$getFittedResponse()$responseMean
      x <- data.frame(time=self$time, observed=observed, fitted=fitted)
      x <- ddply(x, .(time), function(x) data.frame(observed=sum(x$observed), fitted=sum(x$fitted)))
      print(x)
      return(invisible(self))
    },
    
    plotTemporalVariation = function() {
      observed <- self$getObserved()
      fitted <- self$getFittedResponse()$responseMean
      x <- data.frame(time=self$time, observed=observed, fitted=fitted)
      x <- ddply(x, .(time), function(x) data.frame(observed=sum(x$observed), fitted=sum(x$fitted)))
      x <- melt(x, id.vars="time", measure.vars=c("observed", "fitted"))
      p <- ggplot(x, aes(time, value, colour=variable)) + geom_line()
      print(p)
      return(invisible(self))
    },
    
    plotSpatialVariation = function(timeIndex) {
      
    }
  )
)
