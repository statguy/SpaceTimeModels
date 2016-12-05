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
    interceptPrecision = 0.0,
    
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
      obsStack <- inla.stack(data = data, A = A, effects = effects, tag = tag)
      self$fullStack <- if (is.null(self$fullStack)) INLA::inla.stack(obsStack)
      else {
        if (tag %in% names(self$fullStack$data$index))
          stop("Stack with tag '", tag, "' already exists.")
        INLA::inla.stack(self$fullStack, obsStack)
      }
      return(invisible(self))
    },
    
    getSpatialMesh = function() return(self$spaceMesh),
    getSPDE = function() return(self$spde),
    getFullStack = function() return(self$fullStack),
    
    clearStack = function() {
      self$fullStack <- NULL
      return(invisible(self))
    },
    
    setSpatialMesh = function(mesh) {
      if (missing(mesh))
        stop("Argument 'mesh' must be specified.")
      if (!inherits(mesh, "Mesh"))
        stop("Argument 'mesh' must be of class 'SpaceTimeModels::Mesh' or descedant.")
      if (is.null(mesh$getINLAMesh()))
        stop("Mesh has not been initialized.")
      self$spaceMesh <- mesh
      return(invisible(self))
    },
    
    setSPDE = function(spde) {
      if (missing(spde))
        stop("Argument 'spde' must be specified.")
      if (!inherits(spde, "inla.spde2"))
        stop("Argument 'spde' must be of class 'INLA::inla.spde'.")
      self$spde <- spde
      return(invisible(self))
    },
    
    setCovariatesModel = function(covariatesModel, covariates) {
      if (missing(covariatesModel)) covariatesModel <- ~ 1
      else { 
        if (!inherits(covariatesModel, "formula"))
          stop("Argument 'covariatesModel' must be class of 'formula' or descedant.")
      }
      
      x <- if (missing(covariates)) terms(covariatesModel)
      else {
        if (!inherits(covariates, "data.frame"))
          stop("Argument 'covariates' must be class of 'data.frame' or descedant.")
        SpaceTimeModels::assertCompleteCovariates(covariatesModel, covariates)
        terms(covariatesModel, data = covariates)
      }
      
      if (attr(x, "response") != 0)
        stop("The covariates model formula must be right-sided.")
      
      if (!is.null(self$covariatesModel)) {
        warning("Covariates model has been respecified. To enable reuse of the model object, clearStack() method must be called and the data stack needs to be reconstructed.")
      }
      
      self$covariatesModel <- covariatesModel
      if (length(SpaceTimeModels::getCovariateNames(self$covariatesModel)) > 0 && missing(covariates))
        stop("Covariates specified in the model but argument 'covariates' missing.")
      covariateNames <- colnames(SpaceTimeModels::getINLAModelMatrix(covariatesModel, covariates))
      intercept <- if (attr(x, "intercept")[1] == 0) NULL else "intercept"
      randomEffect <- self$getRandomEffectTerm()
      self$linearModel <- reformulate(termlabels = c(intercept, covariateNames, randomEffect), response = "response", intercept = FALSE)
      
      return(invisible(self))
    },
    
    setSpatialPrior = function(rho, rho.init = 0.5, sigma = 1, sigma.init = 0.5) {
      mesh <- self$getSpatialMesh()
      if (is.null(mesh))
        stop("Mesh must be defined first.")
      rho <- ifelse(missing(rho), mesh$getSize() / 2, rho)
      spde <- SpaceTimeModels::local.inla.spde2.matern.new(mesh = mesh$getINLAMesh(), prior.pc.rho = c(rho, rho.init), prior.pc.sig = c(sigma, sigma.init))
      self$setSPDE(spde)
      return(invisible(self))
    },

    setSpatialPriorDefault = function() {
      mesh <- self$getSpatialMesh()
      if (is.null(mesh))
        stop("Mesh must be defined first.")
      spde <- INLA::inla.spde2.matern(mesh = mesh$getINLAMesh())
      self$setSPDE(spde)
      return(invisible(self))
    },
    
    setInterceptPrecision = function(prec = 0.0) {
      self$interceptPrecision <- prec
      return(invisible(self))
    },
    
    getInterceptPrecision = function() return(self$interceptPrecision),
    
    addObservationStack = function(sp, response, covariates, offset, tag = "obs") {
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
        stop("Argument 'sp' must be of class 'sp::SpatialPoints'.")
      if (missing(response))
        stop("Required argument 'response' must be given.")
      
      dataList <- list(response=response)
      if (!missing(offset)) dataList$E <- offset / self$getOffsetScale()
      if (!is.null(self$getLinkFunction())) dataList$link <- self$getLinkFunction()
      
      coordinates <- self$scaleCoordinates(sp::coordinates(sp))
      if (!missing(covariates)) SpaceTimeModels::assertCompleteCovariates(self$covariatesModel, covariates)
      if (length(SpaceTimeModels::getCovariateNames(self$covariatesModel)) > 0 && missing(covariates))
        stop("Covariates specified in the model but argument 'covariates' missing.")
      modelMatrix <- SpaceTimeModels::getINLAModelMatrix(self$covariatesModel, covariates)
      fieldIndex <- INLA::inla.spde.make.index("spatial", n.spde = self$getSPDE()$n.spde)
      A <- INLA::inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc = coordinates)
      
      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept = 1))) else list(fieldIndex)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      } else list(A)
      
      self$addStack(data = dataList, A = AList, effects = effects, tag = tag)
      
      return(invisible(self))
    },
    
    addValidationStack = function(sp, covariates, offset, tag = "val") {
      self$addObservationStack(sp = sp, response = NA, covariates = covariates, offset = offset, tag = tag)
    },
    
    addPredictionStack = function(sp, tag = "pred") {
      if (missing(sp))
        stop("Required argument 'sp' must be given.")
      if (!inherits(sp, "SpatialPoints"))
        stop("Argument 'sp' must be of class 'sp::SpatialPoints'.")
      
      dataList <- list(response = NA)
      if (!is.null(self$getLinkFunction())) dataList$link <- self$getLinkFunction()
      
      coordinates <- self$scaleCoordinates(sp::coordinates(sp))
      fieldIndex <- INLA::inla.spde.make.index("spatial", n.spde = self$getSPDE()$n.spde)

      effects <- if (self$hasIntercept()) list(c(fieldIndex, list(intercept = 1))) else list(c(fieldIndex))
      mesh <- self$getSpatialMesh()
      nodes <- coordinates(mesh$getScaledMeshNodes())
      AList <- list(INLA::inla.spde.make.A(mesh$getINLAMesh(), loc = nodes))
      
      self$addStack(data = dataList, A = AList, effects = effects, tag = tag)
      
      return(invisible(self))
    },
    
    estimate = function(waic = TRUE, dic = FALSE, cpo = FALSE, verbose = FALSE) {
      if (is.null(self$getFullStack()))
        stop("Data stack must be specified first.")
      
      dataStack <- inla.stack.data(self$getFullStack(), spde = self$getSPDE())
      self$result <- try(
        INLA::inla(self$getLinearModel(), family = self$getLikelihood(), data = dataStack, E = dataStack$E,
                   control.predictor = list(A = INLA::inla.stack.A(self$getFullStack()), link = 1, compute = TRUE),
                   control.fixed = list(prec.intercept = self$getInterceptPrecision()),
                   control.compute = list(waic = waic, dic = dic, cpo = cpo, config = TRUE),
                   control.inla = list(reordering = "metis"),
                   verbose = verbose)
        )
      
      if (inherits(self$result, "try-error") || self$result$ok == FALSE)
        stop("Estimation failed. Use verbose=TRUE to find the possible cause.")
      
      return(invisible(self))
    },

    getIndex = function(tag = "obs") {
      if (is.null(self$getFullStack()))
        stop("No index has been specified.")
      
      index <- INLA::inla.stack.index(self$getFullStack(), tag)$data
      if (is.null(index))
        stop(paste("No index found with tag", obs))
      return(index)
    },
    
    getOffset = function(tag = "obs") {
      index <- self$getIndex(tag)
      offset <- INLA::inla.stack.LHS(self$getFullStack())$E[index]
      if (is.null(offset) || all(is.na(offset))) offset <- 1
      offset
    },

    getObserved = function(tag = "obs") {
      index <- self$getIndex(tag)
      INLA::inla.stack.LHS(self$getFullStack())$response[index]
    },

    getFittedResponse = function(variable = "mean", withOffset = FALSE, tag = "obs") {
      index <- self$getIndex(tag)
      offset <- ifelse(withOffset, self$getOffset(tag), 1)
      return(self$getResult()$summary.fitted.values[index, variable] * offset)
    },
    
    getFittedLinearPredictor = function(variable = "mean", tag = "obs") {
      index <- self$getIndex(tag)
      return(self$getResult()$summary.linear.predictor[index, variable])
    },
    
    getFittedSpatialEffect = function(variable = "mean") {
      return(self$getResult()$summary.random$spatial[variable])
    },
    
    getFittedFixedEffects = function() {
      return(self$getResult()$summary.fixed)
    },
    
    getFittedHyperparameters = function() {
      return(self$getResult()$summary.hyperpar)
    },
    
    getDIC = function() {
      return(self$getResult()$dic)
    },
    
    getWAIC = function() {
      return(self$getResult()$waic)
    },
    
    getCPO = function() {
      return(self$getResult()$cpo)
    },
    
    getSPDEResult = function() {
      if (is.null(self$result) || is.null(self$spde))
        stop("The model has not been estimated.")
      return(INLA::inla.spde2.result(self$result, "spatial", self$spde))
    },
    
    summarySpatialParameters = function() {
      spdeResult <- self$getSPDEResult()
      range <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.range.nominal[[1]], coordinatesScale = self$getSpatialMesh()$getScale())
      variance <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.variance.nominal[[1]])
      kappa <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.kappa[[1]], coordinatesScale = 1 / self$getSpatialMesh()$getScale())
      tau <- SpaceTimeModels::summaryINLAParameter(spdeResult$marginals.tau[[1]])
      x <- rbind(kappa = kappa, tau = tau, range = range, variance = variance)
      colnames(x) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
      x
    }
  )
)
