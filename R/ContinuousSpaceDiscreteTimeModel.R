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
    
    addObservationStack = function(sp, response = NA, covariates, offset, tag = "obs") {
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
      
      fieldIndex <- inla.spde.make.index("spatial", n.spde = self$getSPDEObject()$n.spde, n.group = nTime)
      A <- inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc = coordinates, group = timeIndex, n.group = nTime)
      
      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept = 1))) else list(fieldIndex)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      self$addStack(data = dataList, A = AList, effects = effects, tag = tag)
      
      return(invisible(self))
    },
    
    addPredictionStack = function(sp, predTimeIndex, tag = "pred") {
      if (missing(sp))
        stop("Required argument 'sp' missing.")
      if (!inherits(sp, "STI"))
        stop("Argument 'sp' must be of class 'STI'.")
      
      dataList <- list(response=NA)
      if (!is.null(self$getLinkFunction())) dataList$link <- self$getLinkFunction()
      
      coordinates <- self$scaleCoordinates(sp::coordinates(sp)) # TODO: sp not required parameter
      time <- time(sp)
      timeIndex <- time - min(time) + 1
      nTime <- length(unique(timeIndex))
      fieldIndex <- inla.spde.make.index("spatial", n.spde = self$getSPDEObject()$n.spde, n.group = nTime) # TODO: get index from obs stack

      effects <- if (missing(predTimeIndex)) {
        if (self$hasIntercept()) list(c(fieldIndex, list(intercept = 1))) else list(c(fieldIndex))
      }
      else {
        if (self$hasIntercept()) list(c(fieldIndex, coordinates, list(intercept=1))) else list(c(fieldIndex, coordinates))
      }
      AList <- if (missing(predTimeIndex)) {
        list(1)
      }
      else {
        A <- inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), group = predTimeIndex, n.group = nTime)
        list(A)
      }

      self$addStack(data = dataList, A = AList, effects = effects, tag = tag)
      
      #st.est = inla.stack(data = list(y = y),
      #                    + A = list(A, 1),
      #                    + effects = list(c(mesh.index, list(intercept = 1)),
      #                                     + list(cov = covariate)),
      #                    + tag = "est")
      #st.pred = inla.stack(data = list(y = NA),
      #                        + A = list(1),
      #                        + effects = list(c(mesh.index, list(intercept = 1))),
      #                        + tag = "pred")
      
      
      return(invisible(self))
    },
    
    getSPDEResult = function() {
      if (is.null(self$result) || is.null(self$spde))
        stop("The model has not been estimated.")
      return(inla.spde2.result(self$result, "spatial", self$spde))
    },
    
    summarySpatialParameters = function() {
      spdeResult <- self$getSPDEResult()
      range <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.range.nominal[[1]], coordinatesScale = self$getSpatialMesh()$getScale())
      variance <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.variance.nominal[[1]])
      kappa <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.kappa[[1]], coordinatesScale = 1/self$getSpatialMesh()$getScale())
      tau <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.tau[[1]])
      x <- rbind(kappa = kappa, tau = tau, range = range, variance = variance)
      colnames(x) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
      x
    },
    
    summary = function() {
      if (is.null(self$result))
        stop("The model has not been estimated.")
      summary(self$result)
    },

    #getFittedResponse = function(tag="obs") {
    #  data <- super$getFittedResponse(tag=tag)
    #  #data$time <- self$time
    #  return(data)
    #},
    
    #getFittedLinearPredictor = function(tag="obs") {
    #  data <- super$getFittedLinearPredictor(tag=tag)
    #  #data$time <- self$time
    #  return(data)
    #},
    
    #getFittedSpatialEffect = function() {
    #  data <- super$getFittedSpatialEffect()
    #  return(data)
    #}
    
    summaryTemporalVariation = function(variable = "mean", timeIndex, tag = "obs") {
      observed <- self$getObserved(tag = tag)
      fitted <- self$getFittedResponse(variable = variable, tag = tag)
      
      timeIndex <- if (missing(timeIndex)) {
        index <- self$getIndex(tag = tag)
        inla.stack.RHS(self$getFullStack())$spatial.group[index]
      }
      else timeIndex
      
      x <- data.frame(time = timeIndex, observed = observed, fitted = fitted)
      df <- x %>% group_by(time) %>% summarise(observed = sum(observed, na.rm = T), fitted = sum(fitted, na.rm = T))
      return(df)
    },
    
    plotTemporalVariation = function(timeIndex, tag = "obs") {
      x <- self$summaryTemporalVariation(timeIndex = timeIndex, tag = tag)
      p <- x %>% gather(observed, fitted, value = "value", key = "variable") %>% 
        ggplot(aes(time, value, colour = variable)) + geom_line()
      return(p)
    },
    
    plotSpatialVariation = function(variable = "mean", timeIndex, xlim, ylim, dims, tag = "pred") {
      str <- self$getSpatialVariationRaster(variable = variable, timeIndex = timeIndex, tag = tag)
      p <- gplot(str$getLayer(1)) + geom_raster(aes(fill = value))
      return(p)
    },
    
    getSpatialVariationRaster = function(variable = "mean", timeIndex, timeLabels, height = 100, width = 100, tag = "pred") {
      predictedValues <- self$getFittedResponse(variable = variable, tag = tag)
      meshNodes <- self$getSpatialMesh()$getINLAMesh()$n
      maxTimeIndex <- length(na.omit(unique(INLA::inla.stack.data(self$getFullStack())$spatial.group)))
      predictions <- INLA::inla.vector2matrix(predictedValues, nrow = meshNodes, ncol = maxTimeIndex)
      if (!missing(timeIndex)) predictions <- predictions[,timeIndex, drop=F]
      
      str <- SpaceTimeRaster$new(sp = self$getSpatialMesh()$getKnots(), height = height, width = width)
      str$project(self$getSpatialMesh(), predictions, timeLabels=timeLabels)
      return(str)
    }
  )
)
