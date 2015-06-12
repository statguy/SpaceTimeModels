#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass Mesh
#' @export Mesh
#' @keywords internal
Mesh <- R6::R6Class(
  "Mesh",
  private = list(
    knotsScale = NULL,
    knots = NULL,
    mesh = NULL,
    crs = NULL,
    
    getMeshKnots = function() return(unique(private$knots) / self$getScale())
  ),
  public = list(
    initialize = function(knots, knotsScale=1) {
      if (missing(knots))
        stop("Required argument 'knots' missing.")
      if (!inherits(knots, "SpatialPoints"))
        stop("Argument 'knots' must be of class 'SpatialPoints'.")
      private$knotsScale <- knotsScale
      #private$knots <- as.matrix(knots)
      private$knots <- sp::coordinates(obs)
      private$crs <- sp::CRS(sp::proj4string(obs))
    },
    
    getScale = function() return(private$knotsScale),
    getKnots = function() return(sp::SpatialPoints(private$knots, proj4string=self$getCRS())),
    getScaledKnots = function() return(sp::SpatialPoints(private$knots / self$getScale(), proj4string=self$getCRS())),
    getINLAMesh = function() return(private$mesh),
    getCRS = function() return(private$crs),
    
    construct = function(...) {
      stop("Unimplemented abstract method 'construct'.")
    },
    
    plot = function() {
      if (is.null(private$mesh))
        stop("Mesh must be constructed first.")
      plot(self$getINLAMesh())
      points(private$getMeshKnots(), pch='*', col='red')
      return(invisible(self))
    }
  )
)
