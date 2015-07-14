#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass Kernel
#' @export Kernel
#' @keywords internal
Kernel <- R6::R6Class(
  "Kernel",
  private = list(
    scale = NULL,
    resolution = NULL,
    scaledRadius = NULL,
    kernel = NULL,
    
    constructKernelInternal = function(radius, scale) {
      stop("Unimplemented abstract method 'constructKernel'.")
    }
  ),
  public = list(
    constructKernel = function(resolution, scale) {
      private$resolution <- if (inherits(resolution, "RasterLayer")) raster::res(resolution)[1] else resolution
      private$scale = scale
      rasterScale <- private$scale / private$resolution
      private$scaledRadius <- round(private$scale / private$resolution)
      private$kernel <- private$constructKernelInternal(private$scaledRadius, rasterScale)
      return(invisible(self))
    },
      
    getScale = function() return(private$scale),
    getResolution = function() return(private$resolution),
    getScaledRadius = function() return(private$scaledRadius),
    asMatrix = function() return(private$kernel)
  )
)

#' @title Identity kernel for smoothing
#' @description Identity kernel for smoothing.
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass IdentityKernel
#' @export IdentityKernel
IdentityKernel <- R6::R6Class(
  "IdentityKernel",
  inherit = SpaceTimeModels::Kernel,
  private = list(
    constructKernelInternal = function(radius, scale) {
      kernel <- matrix(0, ncol=2*radius+1, nrow=2*radius+1)
      kernel[radius+1, radius+1] <- 1
      return(kernel)
    }
  )
)

#' @title Exponential kernel for smoothing
#' @description Exponential kernel for smoothing.
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass ExponentialKernel
#' @export ExponentialKernel
ExponentialKernel <- R6::R6Class(
  "ExponentialKernel",
  inherit = SpaceTimeModels::Kernel,
  private = list(
    constructKernelInternal = function(radius, scale) {
      x <- matrix(-radius:radius, ncol=2*radius+1, nrow=2*radius+1)
      xy <- sqrt(x^2+t(x)^2)
      xy[xy > radius] <- Inf
      kernel <- exp(-xy / scale)
      return(kernel)
    }
  )
)

#' @title Gaussian kernel for smoothing
#' @description Gaussian kernel for smoothing.
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass GaussianKernel
#' @export GaussianKernel
GaussianKernel <- R6::R6Class(
  "GaussianKernel",
  inherit = SpaceTimeModels::Kernel,
  private = list(
    constructKernelInternal = function(radius, scale) {
      x <- matrix(-radius:radius, ncol=2*radius+1, nrow=2*radius+1)
      xy <- x^2+t(x)^2
      xy[xy > radius] <- Inf
      kernel <- exp(-xy / scale)
      return(kernel)
    }
  )
)
