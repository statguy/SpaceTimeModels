#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass ContinousSpaceTimeModel
#' @export ContinuousSpaceTimeModel
#' @keywords internal
ContinuousSpaceTimeModel <- R6::R6Class(
  "ContinuousSpaceTimeModel",
  lock_objects = FALSE,
  inherit = SpaceTimeModels::ContinuousSpaceModel,
  public = list(
    temporalModel = NULL,
    temporalPrior = NULL,
    
    setTemporalPrior = function(model, prior) {
      stop("Unimplemented abstract method 'setTemporalPrior'.")
    }
  )
)
