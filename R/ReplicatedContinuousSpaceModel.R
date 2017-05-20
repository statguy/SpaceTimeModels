#' @title Replicated continuous space model
#' @description Continuous space model for multiple inpedendent observations at the same location.
#' @references Lindgren, F. & Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software, 63(19).
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo <\email{jvj@@iki.fi}>
#' @exportClass ReplicatedContinuousSpaceModel
#' @export ReplicatedContinuousSpaceModel
ReplicatedContinuousSpaceModel <- R6::R6Class(
  "ContinuousSpaceDiscreteTimeModel",
  lock_objects = FALSE,
  inherit = SpaceTimeModels::ContinuousSpaceTimeModel,
  public = list(
    initialize = function() {
    },
    
    getRandomEffectTerm = function() {
      return("f(spatial, model=spde, replicate=spatial.repl)")
    },
    
    addObservationStack = function(sp, index = as.numeric(sp@data$id), response = NA, covariates, offset, tag = "obs") {
      if (missing(sp))
        stop("Required argument 'sp' missing.")
      if (!inherits(sp, "STI"))
        stop("Argument 'sp' must be of class 'STI'.")
      
      if (is.null(self$getSpatialMesh()))
        stop("Mesh must be defined first.")
      if (is.null(self$getSPDE()))
        stop("Spatial prior must be defined first.")
      if (is.null(self$covariatesModel))
        stop("Covariates model must be defined first.")
      #if (missing(coordinates)) coordinates <- model$getSpatialMesh()$getKnots()
      
      dataList <- list(response = response)
      if (!missing(offset)) dataList$E <- offset / self$getOffsetScale()
      if (!is.null(self$getLinkFunction())) dataList$link <- self$getLinkFunction()
      
      coordinates <- self$scaleCoordinates(sp::coordinates(sp))
      if (!missing(covariates)) SpaceTimeModels::assertCompleteCovariates(self$covariatesModel, covariates)
      if (length(SpaceTimeModels::getCovariateNames(self$covariatesModel)) > 0 && missing(covariates))
        stop("Covariates specified in the model but argument 'covariates' missing.")
      modelMatrix <- SpaceTimeModels::getINLAModelMatrix(self$covariatesModel, covariates)
      
      time <- time(sp)
      timeIndex <- time - min(time) + 1
      nTime <- length(unique(timeIndex))
      
      n <- self$getSpatialMesh()$getNumNodes()
      fieldIndex <- INLA::inla.spde.make.index("spatial", n.spde = self$getSPDE()$n.spde, n.repl = nTime)
      A <- INLA::inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc = coordinates, index = index, repl = timeIndex)
      
      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept = 1))) else list(fieldIndex)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      self$addStack(data = dataList, A = AList, effects = effects, tag = tag)
      
      return(invisible(self))
    },
    
    addPredictionStack = function(sp, tag = "pred") {
      if (missing(sp))
        stop("Required argument 'sp' missing.")
      if (!inherits(sp, "STI"))
        stop("Argument 'sp' must be of class 'STI'.")
      
      dataList <- list(response = NA)
      if (!is.null(self$getLinkFunction())) dataList$link <- self$getLinkFunction()
      
      coordinates <- self$scaleCoordinates(sp::coordinates(sp))
      nTime <- length(unique(time(sp)))
      fieldIndex <- INLA::inla.spde.make.index("spatial", n.spde = self$getSPDE()$n.spde, n.repl = nTime)
      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept = 1))) else list(c(fieldIndex))
      AList <- list(1)
      
      self$addStack(data = dataList, A = AList, effects = effects, tag = tag)
      
      return(invisible(self))
    },
    
    summary = function() {
      if (is.null(self$result))
        stop("The model has not been estimated.")
      summary(self$result)
    }
    
  )
)