#' @title 2D mesh
#' @description Constructs 2D mesh.
#' @references Lindgren, F. & Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software, 63(19).
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass SpatialMesh
#' @export SpatialMesh
SpatialMesh <- R6::R6Class(
  "SpatialMesh",
  lock_objects = FALSE,
  inherit = SpaceTimeModels::Mesh,
  public = list(
    construct = function(cutoff = NULL, maxEdge = NULL, offset = NULL, minAngle = NULL, locDomain = NULL) {
      if (missing(cutoff))
        stop("Required argument 'cutoff' missing.")
      if (missing(maxEdge))
        stop("Required argument 'maxEdge' missing.")
      if (!is.null(locDomain) && !inherits(locDomain, "SpatialPoints"))
        stop("Argument 'locDomain' must be of class SpatialPoints.")
      
      meshCoordinates <- self$getMeshKnots()
      locDomain <- if (!is.null(locDomain)) SpaceTimeModels::nullScale(sp::coordinates(locDomain), self$getScale()) else NULL
      self$mesh <- INLA::inla.mesh.2d(loc = meshCoordinates,
                                      loc.domain = locDomain,
                                      cutoff = SpaceTimeModels::nullScale(cutoff, self$getScale()),
                                      max.edge = SpaceTimeModels::nullScale(maxEdge, self$getScale()),
                                      offset = SpaceTimeModels::nullScale(offset, self$getScale()),
                                      min.angle = minAngle)
    },

    initialize = function(..., cutoff = NULL, maxEdge = NULL, offset = NULL, minAngle = NULL, locDomain = NULL) {
      super$initialize(...)
      self$construct(cutoff = cutoff, maxEdge = maxEdge, offset = offset, minAngle = minAngle, locDomain = locDomain)
    }
  )
)
