#' @title Space-time raster
#' @description Class to hold space time rasters.
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo \email{jvj@@iki.fi}
#' @exportClass SpaceTimeRaster
#' @export SpaceTimeRaster
SpaceTimeRaster <- R6::R6Class(
  "SpaceTimeRaster",
  lock_objects = FALSE,
  public = list(
    template = raster::raster(),
    layers = raster::stack(),

    initialize = function(x, height = 180, width = 360, crs) {
      if (!missing(x)) {
        self$template <- if (inherits(x, "Extent")) raster::raster(x, nrows = height, ncols = width, crs = crs)
        else if (inherits(x, "Spatial")) raster::raster(raster::extent(x), nrows = height, ncols = width, crs = x@proj4string)
        else if (inherits(x, "RasterLayer")) x
        else stop("Parameter 'extent' must be of class 'raster::Raster', 'raster::Extent', 'sp::Spatial' or descedant.")
      }
    },
    
    setLayers = function(layers) {
      self$layers <- layers
      invisible(self)
    },
    
    addLayer = function(layer) {
      self$layers <- raster::addLayer(self$layers, layer)
      invisible(self)
    },
    
    getLayer = function(index) {
      if (index < 1 || index > raster::nlayers(self$layers))
        stop("Parameter 'index' is out of range.")
      return(self$layers[[index]])
    },
    
    getLayers = function() {
      return(self$layers)
    },
    
    project = function(mesh, predictions, timeLabels) {
      if (missing(mesh)) stop("Required parameter 'mesh' missing.")
      if (missing(predictions)) stop("Required parameter 'predictions' missing.")
      
      inlaMesh <- if (inherits(mesh, "Mesh")) mesh$getINLAMesh()
      else if (inherits(mesh, "inla.mesh")) mesh
      else stop("Parameter 'mesh' must be of type 'SpaceTimeModels::Mesh' or 'INLA::inla.mesh'.")
      projector <- INLA::inla.mesh.projector(inlaMesh,
                                             dims = c(ncol(self$template), nrow(self$template)),
                                             xlim = c(xmin(self$template), xmax(self$template)),
                                             ylim = c(ymin(self$template), ymax(self$template)))
      
      if (!inherits(predictions, "matrix")) stop("Parameter 'predictions' must be of type 'matrix'.")
      
      timeLabels <- if (missing(timeLabels)) paste0("t", 1:ncol(predictions)) else timeLabels
      
      for (i in 1:ncol(predictions)) {
        projection <- INLA::inla.mesh.project(projector, predictions[,i])
        raster::values(self$template) <- t(projection[,ncol(projection):1])
        names(self$template) <- timeLabels[i]
        self$addLayer(self$template)
      }
      
      invisible(self)
    },
    
    getColorBreaks = function(n = 6) {
      vmin <- min(raster::minValue(self$layers))
      vmax <- max(raster::maxValue(self$layers))
      return(seq(vmin, vmax, length.out = n))
    }
  )
)
