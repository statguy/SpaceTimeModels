#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass SpaceModel
#' @export SpaceModel
#' @keywords internal
SpaceModel <- R6::R6Class(
  "SpaceModel",
  lock_objects = FALSE,
  private = list(
    offsetScale = 1,
    covariatesModel = NULL,
    linearModel = NULL,
    likelihood = "gaussian",
    result = NULL,
    
    getRandomEffectTerm = function() {
      stop("Unimplemented abstract method 'getRandomEffect'.")
    }
  ),
  public = list(
    initialize = function(offsetScale=1, ...) {
      private$offsetScale <- offsetScale
    },
    
    getDistanceUnit = function() return(distanceUnit),
    getOffsetScale = function() return(private$offsetScale),
    getLikelihood = function() return(private$likelihood),
    getLinearModel = function() return(private$linearModel),
    getResult = function() return(private$result),
            
    setCovariatesModel = function(covariatesModel, covariates) {
      stop("Unimplemented abstract method 'setCovariatesModel'.")
    },
    
    setSmoothingModel = function() {
      return(self$setCovariatesModel(~ 1))
    },
    
    setLikelihood = function(likelihood) {
      if (missing(likelihood))
        stop("Required argument 'likelihood' missing.")
      private$likelihood <- likelihood
      return(invisible(self))
    },

    estimate = function(verbose=T) {
      stop("Unimplemented abstract method 'estimate'.")
    },
    
    save = function(fileName) {
      save(self, file=fileName)
      return(invisible(self))
    },
    
    load = function(fileName) {
      load(fileName, env=self)
      return(invisible(self))
      #tempEnv <- new.env()
      #load(fileName, env=tempEnv)
      #return(invisible(tempEnv$self))
    },

    summary = function() {
      print(summary(self$getResult()))
      return(invisible(self))
    },

    getObserved = function(tag="obs") {
      stop("Unimplemented abstract method 'getFittedObserved'.")
    },
    
    getFittedResponse = function(tag="obs") {
      stop("Unimplemented abstract method 'getFittedResponse'.")
    },
    
    getFittedLinearPredictor = function(tag="obs") {
      stop("Unimplemented abstract method 'getFittedLinearPredictor'.")
    },
    
    getFittedSpatialEffect = function() {
      data <- list()
      data$spatialMean <- self$getResult()$summary.random$spatial$mean
      data$spatialSd <- self$getResult()$summary.random$spatial$sd
      return(as.data.frame(data))
    }
  )
)
