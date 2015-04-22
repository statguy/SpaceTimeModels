#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass SpaceTimeModel
#' @export SpaceTimeModel
#' @keywords internal
SpaceTimeModel <- R6::R6Class(
  "SpaceTimeModel",
  private = list(
    distanceUnit = "m",
    timeUnit = "a"
  ),
  public = list(
    initialize = function(distanceUnit="m", timeUnit="a") {
      private$distanceUnit <- distanceUnit
      private$timeUnit <- timeUnit
    },
    
    getDistanceUnit = function() return(distanceUnit),
    getTimeUnit = function() return(timeUnit),
    
    estimate = function(verbose=T) {
      stop("Unimplemented method 'estimate'.")
    },
    
    save = function(fileName) {
      stop("Unimplemented method 'save'.")
    },
    
    load = function(fileName) {
      load(fileName, env=private)
      return(invisible(self))
    }
  )
)
