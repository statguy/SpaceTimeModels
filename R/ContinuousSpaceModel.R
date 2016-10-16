#' @title Continuous space model
#' @description Building and estimating continuous space models.
#' @references Lindgren, F. & Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software, 63(19).
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass ContinuousSpaceModel
#' @export ContinuousSpaceModel
#' @keywords internal
ContinuousSpaceModel <- R6::R6Class(
  "ContinuousSpaceModel",
  lock_objects = FALSE,
  inherit = SpaceTimeModels::SpaceModel,
  public = list(
    spaceMesh = NULL,
    spde = NULL,
    fullStack = NULL,
    
    scaleCoordinates = function(coordinates) {
      return(as.matrix(coordinates) / self$getSpatialMesh()$getScale())
    },
    
    getRandomEffectTerm = function() {
      return("f(spatial, model=spde)")
    },
    
    hasIntercept = function() {
      if (is.null(self$linearModel))
        stop("Linear model must be defined first.")
      return("intercept" %in% attr(terms(self$linearModel), "term.labels"))
    },
    
    addStack = function(data, A, effects, tag) {
      obsStack <- inla.stack(data=data, A=A, effects=effects, tag=tag)
      self$fullStack <- if (is.null(self$fullStack)) inla.stack(obsStack)
      else {
        if (tag %in% names(self$fullStack$data$index))
          stop("Stack with tag '", tag, "' already exists.")
        inla.stack(self$fullStack, obsStack)
      }
      return(invisible(self))
    },
    
    getSpatialMesh = function() return(self$spaceMesh),
    getSPDEObject = function() return(self$spde),
    getFullStack = function() return(self$fullStack),
    
    clearStack = function() {
      self$fullStack <- NULL
      return(invisible(self))
    },
    
    setSpatialMesh = function(mesh) {
      if (!inherits(mesh, "Mesh"))
        stop("Argument 'mesh' must be of class 'Mesh'.")
      if (is.null(mesh$getINLAMesh()))
        stop("Mesh has not been initialized.")
      self$spaceMesh <- mesh
      return(invisible(self))
    },
    
    setCovariatesModel = function(covariatesModel, covariates) {
      if (missing(covariatesModel)) covariatesModel <- ~ 1
      x <- if (missing(covariates)) terms(covariatesModel)
      else {
        SpaceTimeModels::assertCompleteCovariates(covariatesModel, covariates)
        terms(covariatesModel, data=covariates)
      }
      
      if (attr(x, "response") != 0)
        stop("The covariates model formula must be right-sided.")
      
      if (!is.null(self$covariatesModel)) {
        warning("Covariates model has been respecified. To enable reuse of the model object, clearStack() method must be called and the data stack needs to be reconstructed.")
      }
      
      self$covariatesModel <- covariatesModel
      if (length(SpaceTimeModels::getCovariateNames(self$covariatesModel)) > 0 && missing(covariates))
        stop("Covariates specified in the model but argument 'covariates' missing.")
      covariates <- colnames(getINLAModelMatrix(covariatesModel, covariates))
      intercept <- if (attr(x, "intercept")[1] == 0) NULL else "intercept"
      randomEffect <- self$getRandomEffectTerm()
      self$linearModel <- reformulate(termlabels=c(intercept, covariates, randomEffect), response="response", intercept=FALSE)
      
      return(invisible(self))
    },
    
    setSpatialPrior = function(range, rangeFactor=5) {
      if (is.null(self$getSpatialMesh()))
        stop("Mesh must be defined first.")
      
      sigma0 <- 1
      range0 <- if (missing(range)) {
        size <- min(c(diff(base::range(self$getSpatialMesh()$getINLAMesh()$loc[,1])), diff(base::range(self$getSpatialMesh()$getINLAMesh()$loc[,2]))))
        size / rangeFactor
      }
      else range0 <- range
      kappa0 <- sqrt(8) / range0
      tau0 <- 1 / (sqrt(4 * pi) * kappa0 * sigma0)
      self$spde <- inla.spde2.matern(mesh=self$getSpatialMesh()$getINLAMesh(),
                                        B.tau=cbind(log(tau0), -1, +1),
                                        B.kappa=cbind(log(kappa0), 0, -1),
                                        theta.prior.mean=c(0, 0),
                                        theta.prior.prec=c(0.1, 1))
      
      return(invisible(self))        
    },
    
    addObservationStack = function(sp, response, covariates, offset, tag="obs") {
      # TODO: allow defining link function
      
      if (is.null(self$getSpatialMesh()))
        stop("Mesh must be defined first.")
      if (is.null(self$spde))
        stop("Spatial prior must be defined first.")
      if (is.null(self$covariatesModel))
        stop("Covariates model must be defined first.")
      #if (missing(coordinates)) coordinates <- model$getSpatialMesh()$getKnots()
      if (missing(sp))
        stop("Required argument 'sp' must be given.")
      if (!inherits(sp, "SpatialPoints"))
        stop("Argument 'sp' must be of class 'SpatialPoints'.")
      if (missing(response))
        stop("Required argument 'response' must be given.")
      
      dataList <- list(response=response)
      if (!missing(offset)) dataList$E <- offset / self$getOffsetScale()
      if (!is.null(self$getLinkFunction())) dataList$link <- self$getLinkFunction()
      
      coordinates <- self$scaleCoordinates(sp::coordinates(sp))
      SpaceTimeModels::assertCompleteCovariates(self$covariatesModel, covariates)
      modelMatrix <- SpaceTimeModels::getINLAModelMatrix(self$covariatesModel, covariates)
      fieldIndex <- inla.spde.make.index("spatial", n.spde=self$getSPDEObject()$n.spde)
      A <- inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc=coordinates)
      
      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept=1))) else list(fieldIndex)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      self$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
    
    addValidationStack = function(sp, covariates, offset, tag="val") {
      self$addObservationStack(sp=sp, response=NA, covariates=covariates, offset=offset, tag=tag)
    },
    
    addPredictionStack = function(sp, tag="pred") {
      if (missing(sp))
        stop("Required argument 'sp' must be given.")
      if (!inherits(sp, "SpatialPoints"))
        stop("Argument 'sp' must be of class 'SpatialPoints'.")
      
      dataList <- list(response=NA)
      if (!is.null(self$getLinkFunction())) dataList$link <- self$getLinkFunction()
      
      coordinates <- self$scaleCoordinates(sp::coordinates(sp))
      fieldIndex <- inla.spde.make.index("spatial", n.spde=self$getSPDEObject()$n.spde)
      effects <- if (self$hasIntercept()) list(c(fieldIndex, coordinates, list(intercept=1))) else list(c(fieldIndex, coordinates))
      AList <- list(1)
      
      self$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
        
    estimate = function(verbose=F) {
      if (is.null(self$getFullStack()))
        stop("Data stack must be specified first.")
      
      dataStack <- inla.stack.data(self$getFullStack(), spde=self$getSPDEObject())
      self$result <- try(inla(self$getLinearModel(), family=self$getLikelihood(), data=dataStack, E=dataStack$E,
                                 control.predictor=list(A=inla.stack.A(self$getFullStack()), link=1, compute=TRUE),
                                 control.compute=list(waic=TRUE, config=TRUE),
                                 control.inla=list(reordering="metis"),
                                 verbose=verbose))
      if (inherits(self$result, "try-error") || self$result$ok == FALSE)
        stop("Estimation failed. Use verbose=TRUE to find the possible cause.")
      
      return(invisible(self))
    },

    getIndex = function(tag="obs") {
      if (is.null(self$getFullStack()))
        stop("No index has been specified.")
      
      index <- inla.stack.index(self$getFullStack(), tag)$data
      if (is.null(index))
        stop(paste("No index found with tag", obs))
      return(index)
    },
    
    getOffset = function(tag="obs") {
      index <- self$getIndex(tag)
      offset <- inla.stack.LHS(self$getFullStack())$E[index]
      if (is.null(offset)) offset <- 1
      offset
    },

    getObserved = function(tag="obs") {
      index <- self$getIndex(tag)
      inla.stack.LHS(self$getFullStack())$response[index]
    },

    getFittedResponse = function(tag="obs") {
      index <- self$getIndex(tag)
      offset <- self$getOffset(tag)
      data <- list()
      data$responseMean <- self$getResult()$summary.fitted.values$mean[index] * offset
      data$responseSd <- self$getResult()$summary.fitted.values$sd[index] * offset
      return(as.data.frame(data))
    },
    
    getFittedLinearPredictor = function(tag="obs") {
      index <- self$getIndex(tag)
      data <- list()
      data$etaMean <- self$getResult()$summary.linear.predictor$mean[index]
      data$etaSd <- self$getResult()$summary.linear.predictor$sd[index]
      return(as.data.frame(data))      
    },
    
    getFittedSpatialEffect = function() {
      data <- list()
      data$spatialMean <- self$getResult()$summary.random$spatial$mean
      data$spatialSd <- self$getResult()$summary.random$spatial$sd
      # location, year
      return(as.data.frame(data))
    },
    
    getFittedFixedEffects = function() {
      return(self$getResult()$summary.fixed)
    },
    
    getFittedHyperparameters = function() {
      return(self$getResult()$summary.hyperpar)
    },
    
    getWAIC = function() {
      return(self$getResult()$waic)
    }
  )
)
