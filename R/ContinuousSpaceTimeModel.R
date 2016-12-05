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
    temporalParams = NULL,
    
    setTemporalPrior = function(model, prior, ...) {
      if (missing(model) && missing(prior))
        stop("Required arguments 'model' and/or 'prior' missing.")
      if (!missing(model)) self$temporalModel <- model
      if (!missing(prior)) self$temporalPrior <- prior
      if (!missing(...)) self$temporalParams <- list(...)
      return(invisible(self))
    }
  )
)
