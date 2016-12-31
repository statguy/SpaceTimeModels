#' @title Non-convex hull mesh
#' @description Constructs non-convex hull mesh.
#' @references Lindgren, F. & Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software, 63(19).
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass NonConvexHullMesh
#' @export NonConvexHullMesh
NonConvexHullMesh <- R6::R6Class(
  "NonConvexHullMesh",
  lock_objects = FALSE,
  inherit = SpaceTimeModels::Mesh,
  public = list(
    construct = function(cutoff = NULL, maxEdge = NULL, offset = NULL, minAngle = NULL, innerConvex = -0.15, outerConvex) {
      if (missing(cutoff))
        stop("Required argument 'cutoff' missing.")
      if (missing(maxEdge))
        stop("Required argument 'maxEdge' missing.")
      
      meshCoordinates <- self$getMeshKnots()
      boundary1 <- INLA::inla.nonconvex.hull(points = meshCoordinates, convex = innerConvex)
      boundary <- if (!missing(outerConvex)) {
        boundary2 <- INLA::inla.nonconvex.hull(points = meshCoordinates, convex = outerConvex)
        list(boundary1, boundary2)
      }
      else boundary1
      self$mesh <- INLA::inla.mesh.2d(boundary = boundary,
                                      cutoff = SpaceTimeModels::nullScale(cutoff, self$getScale()),
                                      max.edge = SpaceTimeModels::nullScale(maxEdge, self$getScale()),
                                      offset  = SpaceTimeModels::nullScale(offset, self$getScale()),
                                      min.angle = minAngle)
    },
    
    initialize = function(..., cutoff = NULL, maxEdge = NULL, offset = NULL, minAngle = NULL, innerConvex, outerConvex) {
      super$initialize(...)
      self$construct(cutoff = cutoff, maxEdge = maxEdge, offset = offset, minAngle = minAngle, innerConvex = innerConvex, outerConvex = outerConvex)
    }
  )
)
