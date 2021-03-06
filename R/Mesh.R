#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass Mesh
#' @export Mesh
#' @keywords internal
Mesh <- R6::R6Class(
  "Mesh",
  lock_objects = FALSE,
  public = list(
    knotsScale = NULL,
    knots = NULL,
    mesh = NULL,
    crs = NULL,
    
    #getMeshKnots = function() return(unique(self$knots) / self$getScale()),
    getMeshKnots = function() return(unique(self$knots)),
    
    initialize = function(knots, knotsScale = 1) {
      if (missing(knots))
        stop("Required argument 'knots' missing.")
      if (!inherits(knots, c("SpatialPoints", "ST")))
        stop("Argument 'knots' must be of class 'SpatialPoints' or 'ST'.")
      self$knotsScale <- knotsScale
      self$knots <- sp::coordinates(knots) / knotsScale
      self$crs <- sp::CRS(sp::proj4string(knots))
    },
    
    getScale = function() return(self$knotsScale),
    
    getKnots = function() return(sp::SpatialPoints(self$knots * self$getScale(), proj4string = self$getCRS())),
    #getScaledKnots = function() return(sp::SpatialPoints(self$knots / self$getScale(), proj4string = self$getCRS())),
    getScaledKnots = function() return(sp::SpatialPoints(self$knots, proj4string = self$getCRS())),
    getMeshNodes = function() return(sp::SpatialPoints(self$getINLAMesh()$loc[,1:2] * self$getScale())),
    getScaledMeshNodes = function() return(sp::SpatialPoints(self$getINLAMesh()$loc[,1:2])),
  
    getINLAMesh = function() return(self$mesh),
    getCRS = function() return(self$crs),
    getNumNodes = function() return(self$getINLAMesh()$n),
    
    getRange = function() {
      nodes <- coordinates(self$getScaledMeshNodes())
      return(c(diff(base::range(nodes[,1])), diff(base::range(nodes[,2]))))
    },
    
    getSize = function() return(min(self$getRange())),
    
    plot = function() {
      if (is.null(self$mesh))
        stop("Mesh must be constructed first.")
      plot(self$getINLAMesh())
      points(self$getMeshKnots(), pch = '*', col = 'red')
      return(invisible(self))
    }
  )
)
