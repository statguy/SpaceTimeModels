#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass ContinousSpaceTimeModel
#' @export ContinuousSpaceTimeModel
#' @keywords internal
ContinuousSpaceTimeModel <- R6::R6Class(
  "SpaceTimeModel",
  lock_objects = FALSE,
  inherit = SpaceTimeModels::ContinuousSpaceModel,
  public = list(
  )
)
