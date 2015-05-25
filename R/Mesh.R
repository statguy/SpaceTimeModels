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
    
    getMeshKnots = function() return(unique(self$getKnots()) / self$getScale())
  ),
  public = list(
    initialize = function(knots, knotsScale=1) {
      if (missing(knots))
        stop("Required argument 'knots' missing.")
      private$knotsScale <- knotsScale
      private$knots <- as.matrix(knots)
    },
    
    getScale = function() return(private$knotsScale),
    getKnots = function() return(private$knots),
    getScaledKnots = function() return(self$getKnots() / self$getScale()),
    getINLAMesh = function() return(private$mesh),
    
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
