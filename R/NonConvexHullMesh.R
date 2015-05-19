#' @title Non-convex hull mesh
#' @description Constructs non-convex hull mesh.
#' @references Lindgren, F. & Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software, 63(19).
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass NonConvexHullMesh
#' @export NonConvexHullMesh
#' @keywords internal
NonConvexHullMesh <- R6::R6Class(
  "NonConvexHullMesh",
  inherit = SpaceTime::Mesh,
  private = list(
  ),
  public = list(
    construct = function(cutoff=NULL, maxEdge=NULL, offset=NULL, minAngle=NULL, convex) {
      if (missing(cutoff))
        stop("Required argument 'cutoff' missing.")
      if (missing(maxEdge))
        stop("Required argument 'maxEdge' missing.")
          
      meshCoordinates <- private$getMeshKnots()
      boundary <- inla.nonconvex.hull(points=meshCoordinates, convex=convex)
      private$mesh <- inla.mesh.2d(boundary=boundary,
                                   cutoff=SpaceTime::nullScale(cutoff, self$getScale()),
                                   max.edge=SpaceTime::nullScale(maxEdge, self$getScale()),
                                   offset=SpaceTime::nullScale(offset, self$getScale()),
                                   min.angle=minAngle)
      
      return(invisible(self))
    }
  )
)
