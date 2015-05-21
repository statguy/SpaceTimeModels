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
  inherit = SpaceTime::SpaceModel,
  private = list(
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
      if (is.null(private$linearModel))
        stop("Linear model must be defined first.")
      return("intercept" %in% attr(terms(private$linearModel), "term.labels"))
    },
    
    addStack = function(data, A, effects, tag) {
      obsStack <- inla.stack(data=data, A=A, effects=effects, tag=tag)
      private$fullStack <- if (is.null(private$fullStack)) inla.stack(obsStack)
      else {
        if (names(private$fullStack$data$index) == tag)
          warning("Stack with tag '", tag, "' already exists. Overwriting...")
        # TODO: this doesn't exactly overwrite?
        inla.stack(private$fullStack, obsStack)
      }
      return(invisible(self))
    }
  ),
  public = list(
    getSpatialMesh = function() return(private$spaceMesh),
    getSPDEObject = function() return(private$spde),
    getFullStack = function() return(private$fullStack),
    
    setSpatialMesh = function(mesh) {
      if (!inherits(mesh, "Mesh"))
        stop("Argument 'mesh' must be of class 'Mesh'.")
      private$spaceMesh <- mesh
      return(invisible(self))
    },
    
    setCovariatesModel = function(covariatesModel, covariates) {
      if (missing(covariatesModel)) covariatesModel <- ~ 1
      x <- if (missing(covariates)) terms(covariatesModel)
      else terms(covariatesModel, data=covariates)
      
      if (attr(x, "response") != 0)
        stop("The covariates model formula must be right-sided.")
      
      private$covariatesModel <- covariatesModel
      covariates <- colnames(getINLAModelMatrix(covariatesModel, covariates))
      intercept <- if (attr(x, "intercept")[1] == 0) NULL else "intercept"
      randomEffect <- private$getRandomEffectTerm()
      private$linearModel <- reformulate(termlabels=c(intercept, covariates, randomEffect), response="response", intercept=FALSE)
      
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
      private$spde <- inla.spde2.matern(mesh=self$getSpatialMesh()$getINLAMesh(),
                                        B.tau=cbind(log(tau0), -1, +1),
                                        B.kappa=cbind(log(kappa0), 0, -1),
                                        theta.prior.mean=c(0, 0),
                                        theta.prior.prec=c(0.1, 1))
      
      return(invisible(self))        
    },
    
    addObservationStack = function(coordinates, response=NA, covariates, offset, offsetScale, tag="obs") {
      # TODO: allow defining link function
      
      if (missing(coordinates))
        stop("Required argument 'coordinates' missing.")
      if (is.null(self$getSpatialMesh()))
        stop("Mesh must be defined first.")
      if (is.null(private$spde))
        stop("Spatial prior must be defined first.")
      if (is.null(private$covariatesModel))
        stop("Covariates model must be defined first.")
      
      dataList <- list(response=response)
      if (!missing(offset)) {
        private$offsetScale <- if (missing(offsetScale)) SpaceTime::findScale(offset[1])
        else offsetScale
        dataList$E <- offset / private$offsetScale
      }
      
      coordinates <- private$scaleCoordinates(coordinates)
      modelMatrix <- SpaceTime::getINLAModelMatrix(private$covariatesModel, covariates)
      fieldIndex <- inla.spde.make.index("spatial", n.spde=self$getSPDEObject()$n.spde)
      A <- inla.spde.make.A(self$getSpatialMesh()$getINLAMesh(), loc=coordinates)
      
      effects <- if (private$hasIntercept()) list(c(fieldIndex, list(intercept=1))) else list(fieldIndex)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      private$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
    
    addPredictionStack = function(coordinates, response=NA, covariates, offset, offsetScale, tag="pred") {
      # TODO: finish this function
      
      effects <- if (private$hasIntercept()) list(c(fieldIndex, coordinates, list(intercept=1))) else list(c(index, coordinates))
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      private$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
        
    estimate = function(verbose=F) {
      if (is.null(self$getFullStack()))
        stop("Data stack must be defined first.")
      
      dataStack <- inla.stack.data(self$getFullStack(), spde=self$getSPDEObject())
      private$result <- try(inla(self$getLinearModel(), family=self$getLikelihood(), data=dataStack, E=dataStack$E,
                                 control.predictor=list(A=inla.stack.A(self$getFullStack()), compute=TRUE), # link=dataStack$link
                                 control.compute=list(waic=TRUE, config=TRUE),
                                 verbose=verbose))
      if (inherits(private$result, "try-error") || private$result$ok == FALSE)
        stop("Estimation failed. Use verbose=TRUE to find the possible cause.")
      
      return(invisible(self))
    },
    
    save = function(fileName) {
      stop("Unimplemented abstract method 'save'.")
    }
  )
)
