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
    initialize = function() {
      self$temporalModel <- "ar1"
    },
 
    getRandomEffectTerm = function() {
      model <- if (is.null(self$temporalPrior))
        "f(spatial, model=spde, group=spatial.group, control.group=c(list(model=self$temporalModel)"
      else
        "f(spatial, model=spde, group=spatial.group, control.group=c(list(model=self$temporalModel, hyper=self$temporalPrior)"
      if (!is.null(self$temporalParams)) model <- paste0(model, ", self$temporalParams")
      model <- paste0(model, "))")
      return(model)
    },
    
    addObservationStack = function(sp, response = NA, covariates, offset, tag = "obs") {
      # TODO: allow defining link function
      
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
      
      fieldIndex <- INLA::inla.spde.make.index("spatial", n.spde = self$getSPDE()$n.spde, n.group = nTime)
      A <- INLA::inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc = coordinates, group = timeIndex, n.group = nTime)
      
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
      fieldIndex <- inla.spde.make.index("spatial", n.spde = self$getSPDE()$n.spde, n.group = nTime)
      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept = 1))) else list(c(fieldIndex))
      AList <- list(1)
      
      self$addStack(data = dataList, A = AList, effects = effects, tag = tag)
      
      return(invisible(self))
    },
    
    summary = function() {
      if (is.null(self$result))
        stop("The model has not been estimated.")
      summary(self$result)
    },

    summaryTemporalVariation = function(variable = "mean", timeIndex, summariseFun = sum, tag = "obs", ...) {
      observed <- self$getObserved(tag = tag)
      fitted <- self$getFittedResponse(variable = variable, withOffset = TRUE, tag = tag)
      
      timeIndex <- if (missing(timeIndex)) {
        index <- self$getIndex(tag = tag)
        INLA::inla.stack.RHS(self$getFullStack())$spatial.group[index]
      }
      else timeIndex
      
      offset <- self$getOffset(tag)
      x <- data.frame(time = timeIndex, observed = observed, fitted = fitted, offset = offset)
      df <- x %>% dplyr::group_by(time) %>% 
        dplyr::summarise(observed = summariseFun(observed, ...), fitted = summariseFun(fitted, ...))
      return(df)
    },
    
    plotTemporalVariation = function(timeIndex, tag = "obs") {
      x <- self$summaryTemporalVariation(timeIndex = timeIndex, tag = tag)
      p <- x %>% tidyr::gather(observed, fitted, value = "value", key = "variable") %>% 
        ggplot2::ggplot(aes(time, value, colour = variable)) + ggplot2::geom_line()
      return(p)
    },
    
    getSpatialVariationRaster = function(variable = "mean", timeIndex, timeLabels, template = self$getSpatialMesh()$getKnots(), height = 100, width = 100, crs = self$getSpatialMesh()$getCRS(), tag = "pred") {
      predictedValues <- self$getFittedResponse(variable = variable, tag = tag)
      meshNodes <- self$getSpatialMesh()$getINLAMesh()$n
      maxTimeIndex <- length(na.omit(unique(INLA::inla.stack.data(self$getFullStack())$spatial.group)))
      predictions <- INLA::inla.vector2matrix(predictedValues, nrow = meshNodes, ncol = maxTimeIndex)
      if (!missing(timeIndex)) predictions <- predictions[,timeIndex, drop = F]
      
      r <- SpaceTimeModels::SpaceTimeRaster$new(x = template, height = height, width = width, crs = crs)
      r$project(self$getSpatialMesh(), predictions, timeLabels = timeLabels)
      return(r)
    },
    
    plotSpatialVariation = function(variable = "mean", timeIndex, xlim, ylim, dims, tag = "pred") {
      str <- self$getSpatialVariationRaster(variable = variable, timeIndex = timeIndex, tag = tag)
      p <- rasterVis::gplot(str$getLayer(1)) + ggplot2::geom_raster(aes(fill = value))
      return(p)
    }
  )
)
