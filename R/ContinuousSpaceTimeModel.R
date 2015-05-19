#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass ContinousSpaceTimeModel
#' @export ContinuousSpaceTimeModel
#' @keywords internal
ContinuousSpaceTimeModel <- R6::R6Class(
  "SpaceTimeModel",
  inherit = SpaceTime::ContinuousSpaceModel,
  private = list(
    timeUnit = "a"
  ),
  public = list(
    initialize = function(distanceUnit="m", timeUnit="a") {
      super$initialize(distanceUnit=distanceUnit)
      private$timeUnit <- timeUnit
    },
    
    getTimeUnit = function() return(timeUnit)    
  )
)
