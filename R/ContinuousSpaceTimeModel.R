#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass ContinousSpaceTimeModel
#' @export ContinuousSpaceTimeModel
#' @keywords internal
ContinuousSpaceTimeModel <- R6::R6Class(
  "SpaceTimeModel",
  inherit = SpaceTimeModels::ContinuousSpaceModel,
  private = list(
  ),
  public = list(
  )
)
