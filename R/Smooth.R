#' @import raster
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @export smoothDiscreteSubset
#' @keywords internal
smoothDiscreteSubset <- function(r, x, y, kernel, smoothValues, edgeValues) {
  col <- colFromX(r, x)
  row <- rowFromY(r, y)  
  if (r[row, col] %in% edgeValues || is.na(r[row, col])) {
    warning("The point is outside the effective area. The smoothing cannot be proceeded.")
    return(data.frame(x=x, y=y, scale=kernel$getScale(), value=NA))
  }
  
  k <- kernel$asMatrix()
  # Cut the kernel if partially outside the effective area
  kernelRadius1 <- kernel$getScaledRadius() + 1
  startRow <- max(0, row-kernelRadius1) + 1
  startCol <- max(0, col-kernelRadius1) + 1
  nrows <- min(dim(r)[1]+1, row+kernelRadius1) - startRow
  ncols <- min(dim(r)[2]+1, col+kernelRadius1) - startCol
  xmin <- startRow + (row-kernelRadius1) - 1
  ymin <- startCol + (col-kernelRadius1) - 1
  if (!(dim(k)[1] == nrows & dim(k)[2] == ncols)) {
    #message("Cut kernel ", dim(k)[1], " X ", dim(k)[2], " to ", nrows, " X ", ncols)
    k <- k[1:nrows+xmin, 1:ncols+ymin]
  }
  
  # Get edges and process values around the specified point that matches the size of the kernel
  edgeRaster <- raster::getValuesBlock(r, startRow, nrows, startCol, ncols)
  # Set kernel zero at edges for edge correction
  k[edgeRaster %in% edgeValues | is.na(edgeRaster)] <- 0
  # Rescale kernel to constraint smooth values between 0...1
  k <- k / sum(k)
  # Get binary mask of the values to be smoothed
  smoothMaskRaster <- edgeRaster %in% smoothValues
  # Find convolution
  smoothValue <- sum(k * processRaster, na.rm=T)
  x <- data.frame(x=x, y=y, scale=kernel$getScale(), value=smoothValue)
  return(x)
}

#' @import raster
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @export smoothDiscreteSubsets
smoothDiscreteSubsets <- function(r, coords, kernel, scales, processValues, edgeValues, wide=T, .parallel=T) {
  library(raster)
  library(plyr)
  library(sp)
  library(tidyr)
  
  if (missing(r))
    stop("Argument 'r' missing.")
  if (!inherits(r, "RasterLayer"))
    stop("Argument 'r' must be of type 'RasterLayer'")
  if (missing(coords))
    stop("Argument 'coords' missing.")
  
  if (!inherits(r, "RasterLayer"))
    stop("Argument 'r' must of class 'raster'.")
  if (!inherits(coords, "matrix") && !inherits(coords, "data.frame") && !inherits(coords, "SpatialPoints"))
    stop("Argument 'coords' must of class 'matrix', 'data.frame' or 'SpatialPoints'.")
  
  if (res(r)[1] != res(r)[2])
    stop("Rasters of unequal resolution unsupported.")
  if (nlayers(r) > 1)
    stop("Multiband rasters unsupported. Please supply bands separately.")
  
  if (inherits(coords, "SpatialPoints")) coords <- coordinates(coords)
  
  if (!inherits(kernel, "Kernel"))
    stop("Argument 'kernel' must be of class 'Kernel'.")
  
  library(plyr)
  smoothPixels <- plyr::ldply(rev(sort(scales)), function(scale, coords) {
    n.coords <- nrow(coords)
    kernel$constructKernel(resolution=r, scale=scale)
    message("Kernel size = ", dim(kernel$asMatrix())[1], " X ", dim(kernel$asMatrix())[2])
    smoothPixels <- plyr::ldply(1:n.coords, function(i, coords, n.coords, scale, kernel) {
      message("Smoothing scale = ", scale, ", for coord = ", i, "/", n.coords, " (", coords[i,1], ",", coords[i,2], ")")
      x <- smoothDiscreteSubset(r=r, x=coords[i,1], y=coords[i,2], kernel=kernel, processValues=processValues, edgeValues=edgeValues)
      return(x)
    }, coords=coords, n.coords=n.coords, scale=scale, kernel=kernel, .parallel=.parallel) # Inner loop parallel strategy slower for small kernels, but faster for big kernels
    return(smoothPixels)
  }, coords=coords)
  
  if (wide) smoothPixels <- tidyr::spread(smoothPixels, scale, value)
  return(smoothPixels)
}

#' @import raster
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @export smoothContinuousSubset
#' @keywords internal
smoothContinuousSubset <- function(r, x, y, kernel, edgeValues=c()) {
  col <- colFromX(r, x)
  row <- rowFromY(r, y)  
  if (r[row, col] %in% edgeValues || is.na(r[row, col])) {
    warning("The point is outside the effective area. The smoothing cannot be proceeded.")
    return(data.frame(x=x, y=y, scale=kernel$getScale(), value=NA))
  }
  
  k <- kernel$asMatrix()
  # Cut the kernel if partially outside the effective area
  kernelRadius1 <- kernel$getScaledRadius() + 1
  startRow <- max(0, row-kernelRadius1) + 1
  startCol <- max(0, col-kernelRadius1) + 1
  nrows <- min(dim(r)[1]+1, row+kernelRadius1) - startRow
  ncols <- min(dim(r)[2]+1, col+kernelRadius1) - startCol
  xmin <- startRow - (row-kernelRadius1) - 1
  ymin <- startCol - (col-kernelRadius1) - 1
  if (!(dim(k)[1] == nrows & dim(k)[2] == ncols)) {
    message("Cut kernel ", dim(k)[1], " X ", dim(k)[2], " to ", nrows, " X ", ncols)
    k <- k[1:nrows+xmin, 1:ncols+ymin]
  }
  
  # Get edges and process values around the specified point that matches the size of the kernel
  processRaster <- raster::getValuesBlock(r, startRow, nrows, startCol, ncols)
  # Set kernel zero at edges for edge correction
  k[processRaster %in% edgeValues | is.na(processRaster)] <- 0
  # Rescale kernel to constraint smooth values between 0...1
  k <- k / sum(k)
  # Find convolution
  smoothValue <- sum(k * processRaster, na.rm=T)
  x <- data.frame(x=x, y=y, scale=kernel$getScale(), value=smoothValue)
  return(x)
}

#' @import raster
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @export smoothContinuousSubsets
smoothContinuousSubsets <- function(r, coords, kernel, scales, edgeValues=c(), wide=T, .parallel=T) {
  library(raster)
  library(plyr)
  library(sp)
  library(tidyr)
  
  if (missing(r))
    stop("Argument 'r' missing.")
  if (!inherits(r, "RasterLayer"))
    stop("Argument 'r' must be of type 'RasterLayer'")
  if (missing(coords))
    stop("Argument 'coords' missing.")
  
  if (!inherits(r, "RasterLayer"))
    stop("Argument 'r' must of class 'raster'.")
  if (!inherits(coords, "matrix") && !inherits(coords, "data.frame") && !inherits(coords, "SpatialPoints"))
    stop("Argument 'coords' must of class 'matrix', 'data.frame' or 'SpatialPoints'.")
  
  if (res(r)[1] != res(r)[2])
    stop("Rasters of unequal resolution unsupported.")
  if (nlayers(r) > 1)
    stop("Multiband rasters unsupported. Please supply bands separately.")
  
  if (inherits(coords, "SpatialPoints")) coords <- coordinates(coords)
  
  if (!inherits(kernel, "Kernel"))
    stop("Argument 'kernel' must be of class 'Kernel'.")
  
  library(plyr)
  smoothPixels <- plyr::ldply(rev(sort(scales)), function(scale, coords) {
    n.coords <- nrow(coords)
    kernel$constructKernel(resolution=r, scale=scale)
    message("Kernel size = ", dim(kernel$asMatrix())[1], " X ", dim(kernel$asMatrix())[2])
    smoothPixels <- plyr::ldply(1:n.coords, function(i, coords, n.coords, scale, kernel) {
      message("Smoothing scale = ", scale, ", for coord = ", i, "/", n.coords, " (", coords[i,1], ",", coords[i,2], ")")
      x <- smoothContinuousSubset(r=r, x=coords[i,1], y=coords[i,2], kernel=kernel, edgeValues=edgeValues)
      return(x)
    }, coords=coords, n.coords=n.coords, scale=scale, kernel=kernel, .parallel=.parallel)
    return(smoothPixels)
  }, coords=coords)
  
  if (wide) smoothPixels <- tidyr::spread(smoothPixels, scale, value)
  return(smoothPixels)
}
