#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass SpaceModel
#' @export SpaceModel
#' @keywords internal
SpaceModel <- R6::R6Class(
  "SpaceModel",
  private = list(
    distanceUnit = "m",
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
    initialize = function(distanceUnit="m") {
      private$distanceUnit <- distanceUnit
    },
    
    getDistanceUnit = function() return(distanceUnit),
    getOffsetScale = function() return(private$offsetScale),
    #getOffset = function() return(private$offset * self$getOffsetScale()),
    getLikelihood = function() return(private$likelihood),
    getLinearModel = function() return(private$linearModel),
    getResult = function() return(private$result),
            
    setCovariatesModel = function(covariatesModel, covariates) {
      stop("Unimplemented abstract method 'setCovariatesModel'.")
#       
#       if (missing(covariatesModel)) covariatesModel <- ~ 1
#       x <- if (missing(covariates)) terms(covariatesModel)
#       else terms(covariatesModel, data=covariates)
#       
#       if (attr(x, "response") != 0)
#         stop("The covariates model formula must be right-sided.")
#       
#       private$covariatesModel <- covariatesModel
#       covariates <- colnames(getINLAModelMatrix(covariatesModel, covariates))
#       intercept <- if (attr(x, "intercept")[1] == 0) FALSE else TRUE
#       randomEffect <- private$getRandomEffectTerm()
#       private$linearModel <- reformulate(termlabels=c(covariates, randomEffect), response="response", intercept=intercept)
#       
#       return(invisible(self))
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
      stop("Unimplemented abstract method 'save'.")
    },
    
    load = function(fileName) {
      load(fileName, env=private)
      return(invisible(self))
    }
  )
)
