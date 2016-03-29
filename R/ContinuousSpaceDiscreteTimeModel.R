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
  inherit = SpaceTimeModels::ContinuousSpaceTimeModel,
  private = list(
    getRandomEffectTerm = function() {
      return("f(spatial, model=spde, group=spatial.group, control.group=list(model=\"ar1\"))")
    }
  ),
  public = list(
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
      if (is.null(private$covariatesModel))
        stop("Covariates model must be defined first.")
      #if (missing(coordinates)) coordinates <- model$getSpatialMesh()$getKnots()
      
      dataList <- list(response=response)
      if (!missing(offset)) dataList$E <- offset / self$getOffsetScale()

      coordinates <- private$scaleCoordinates(sp::coordinates(sp))
      SpaceTimeModels::assertCompleteCovariates(private$covariatesModel, covariates)
      modelMatrix <- SpaceTimeModels::getINLAModelMatrix(private$covariatesModel, covariates)
      time <- time(sp)
      timeIndex <- time - min(time) + 1
      nTime <- length(unique(timeIndex))
      fieldIndex <- inla.spde.make.index("spatial", n.spde=self$getSPDEObject()$n.spde, n.group=nTime)
      A <- inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc=coordinates, group=timeIndex, n.group=nTime)
      
      effects <- if (private$hasIntercept()) list(c(fieldIndex, list(intercept=1))) else list(fieldIndex)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      private$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
    
    addPredictionStack = function(sp, response=NA, covariates, tag="pred") {
      # TODO: finish this function
      
      if (missing(sp))
        stop("Required argument 'sp' missing.")
      if (!inherits(sp, "STI"))
        stop("Argument 'sp' must be of class 'STI'.")
      coordinates <- sp::coordinates(sp)
      
      effects <- if (private$hasIntercept()) list(c(fieldIndex, coordinates, list(intercept=1))) else list(c(index, coordinates))
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      private$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
    
    getSPDEResult = function() {
      if (is.null(private$result) || is.null(private$spde))
        stop("The model has not been estimated.")
      return(inla.spde2.result(private$result, "spatial", private$spde))
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
      if (is.null(private$result))
        stop("The model has not been estimated.")
      print(summary(private$result))
      return(invisible(self))
    },
    
    #getFitted = function() {
    #  indexObserved <- inla.stack.index(private$fullStack, "obs")$data
    #  fitted <- private$result$summary.fitted.values$mean[indexObserved] * private$offsetScale
    #  return(fitted)      
    #},
    
    summaryTemporalVariation = function() {
      observed <- self$getObserved()
      fitted <- self$getFittedResponse()$responseMean
      x <- data.frame(time=private$time, observed=observed, fitted=fitted)
      x <- ddply(x, .(time), function(x) data.frame(observed=sum(x$observed), fitted=sum(x$fitted)))
      print(x)
      return(invisible(self))
    },
    
    plotTemporalVariation = function() {
      observed <- self$getObserved()
      fitted <- self$getFittedResponse()$responseMean
      x <- data.frame(time=private$time, observed=observed, fitted=fitted)
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
