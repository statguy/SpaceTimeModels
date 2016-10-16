#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass SpaceModel
#' @export SpaceModel
#' @keywords internal
SpaceModel <- R6::R6Class(
  "SpaceModel",
  lock_objects = FALSE,
  public = list(
    offsetScale = 1,
    covariatesModel = NULL,
    linearModel = NULL,
    likelihood = "gaussian",
    link = NULL,
    result = NULL,
    
    getRandomEffectTerm = function() {
      stop("Unimplemented abstract method 'getRandomEffect'.")
    },
    
    initialize = function(offsetScale=1, ...) {
      self$offsetScale <- offsetScale
    },
    
    getDistanceUnit = function() return(distanceUnit),
    getOffsetScale = function() return(self$offsetScale),
    getLikelihood = function() return(self$likelihood),
    getLinkFunction = function() return(self$link),
    getLinearModel = function() return(self$linearModel),
    getResult = function() return(self$result),
            
    setCovariatesModel = function(covariatesModel, covariates) {
      stop("Unimplemented abstract method 'setCovariatesModel'.")
    },
    
    setSmoothingModel = function() {
      return(self$setCovariatesModel(~ 1))
    },
    
    setLikelihood = function(likelihood) {
      if (missing(likelihood))
        stop("Required argument 'likelihood' missing.")
      self$likelihood <- likelihood
      return(invisible(self))
    },

    setLinkFunction = function(link) {
      if (missing(link))
        stop("Required argument 'link' missing.")
      self$link <- link
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
