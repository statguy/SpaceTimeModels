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
  "TwoDMesh",
  inherit = SpaceTime::Mesh,
  private = list(
  ),
  public = list(
    construct = function(cutoff=NULL, maxEdge=NULL, offset=NULL, minAngle=NULL, locDomain=NULL) {
      if (missing(cutoff))
        stop("Required argument 'cutoff' missing.")
      if (missing(maxEdge))
        stop("Required argument 'maxEdge' missing.")
      
      meshCoordinates <- private$getMeshKnots()
      locDomain <- SpaceTime::nullScale(locDomain, self$getScale())
      private$mesh <- inla.mesh.2d(loc=meshCoordinates,
                                   loc.domain=locDomain,
                                   cutoff=SpaceTime::nullScale(cutoff, self$getScale()),
                                   max.edge=SpaceTime::nullScale(maxEdge, self$getScale()),
                                   offset=SpaceTime::nullScale(offset, self$getScale()),
                                   min.angle=minAngle)
      
      return(invisible(self))
    }
  )
)
