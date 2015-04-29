#' @title Discrete time, continuous space model
#' @description Building and estimating discrete time, continuous space models.
#' @references Lindgren, F. & Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software, 63(19).
#' @usage NULL
#' @format NULL
#' @import R6
#' @author Jussi Jousimo <\email{jvj@@iki.fi}>
#' @exportClass DiscreteTimeContinuousSpaceModel
#' @export DiscreteTimeContinuousSpaceModel
DiscreteTimeContinuousSpaceModel <- R6::R6Class(
  "DiscreteTimeContinuousSpaceModel",
  inherit = SpaceTime::SpaceTimeModel,
  private = list(
    coordinatesScale = 1,
    offsetScale = 1,
    coordinates = NULL,
    mesh = NULL,
    spde = NULL,
    covariatesModel = NULL,
    linearModel = NULL,   
    fullStack = NULL,
    likelihood = "gaussian",
    result = NULL,
    
    hasIntercept = function() {
      if (is.null(private$linearModel))
        stop("Linear model must be defined first.")
      return("intercept" %in% attr(terms(private$linearModel), "term.labels"))
    },
    
    addStack = function(data, A, effects, tag) {
      obsStack <- inla.stack(data=data, A=A, effects=effects, tag=tag)
      private$fullStack <- if (is.null(private$fullStack)) inla.stack(obsStack)
      else {
        if (names(private$fullStack$data$index) == tag)
          warning("Stack with tag '", tag, "' already exists. Overwriting...")
        inla.stack(private$fullStack, obsStack)
      }
      return(invisible(self))
    }
  ),
  public = list(
    getCoordinatesScale = function() return(private$coordinatesScale),
    getOffsetScale = function() return(private$offsetScale),
    getCoordinates = function() return(private$coordinates),
    getScaledCoordinates = function() return(self$getCoordinates() / self$getCoordinatesScale()),
    getTime = function() return(private$time),
    getOffset = function() return(private$offset * self$getOffsetScale()),
    getResponse = function() return(private$response),
    getCovariates = function() return(private$covariates),
    
    getMesh = function() return(private$mesh),
    getSPDEObject = function() return(private$spde),
    getModelMatrix = function() return(private$modelMatrix),
    getLinearModel = function() return(private$linearModel),
    getFullStack = function() return(private$fullStack),
    getResult = function() return(private$result),
 
    constructMesh = function(coordinates, coordinatesScale, cutoff=NULL, maxEdge=NULL, offset=NULL, minAngle=NULL, locDomain=NULL, convex) {
      if (missing(coordinates))
        stop("Required argument 'coordinates' missing.")
      if (missing(cutoff))
        stop("Required argument 'cutoff' missing.")
      if (missing(maxEdge))
        stop("Required argument 'maxEdge' missing.")
      
      private$coordinatesScale <- if (missing(coordinatesScale)) SpaceTime::findScale(coordinates[1,1]) else coordinatesScale
      private$coordinates <- as.matrix(coordinates)
      meshCoordinates <- unique(self$getCoordinates()) / self$getCoordinatesScale()
      
      if (missing(convex)) {
        locDomain <- SpaceTime::nullScale(locDomain, self$getCoordinatesScale())
        private$mesh <- inla.mesh.2d(loc=meshCoordinates, loc.domain=locDomain, cutoff=SpaceTime::nullScale(cutoff, self$getCoordinatesScale()), max.edge=SpaceTime::nullScale(maxEdge, self$getCoordinatesScale()), offset=SpaceTime::nullScale(offset, self$getCoordinatesScale()), min.angle=minAngle)
      }
      else {
        boundary <- inla.nonconvex.hull(points=meshCoordinates, convex=convex)
        private$mesh <- inla.mesh.2d(boundary=boundary, cutoff=SpaceTime::nullScale(cutoff, self$getCoordinatesScale()), max.edge=SpaceTime::nullScale(maxEdge, self$getCoordinatesScale()), offset=SpaceTime::nullScale(offset, self$getCoordinatesScale()), min.angle=minAngle)
      }
      
      return(invisible(self))
    },
    
    plotMesh = function() {
      if (is.null(private$mesh))
        stop("Mesh must be defined first.")
      
      plot(self$getMesh())
      points(self$getScaledCoordinates(), pch='*', col='red')
      return(invisible(self))
    },
    
    setSpatialPrior = function(range, rangeFactor=5) {
      if (is.null(private$mesh))
        stop("Mesh must be defined first.")
      
      sigma0 <- 1
      range0 <- if (missing(range)) {
        size <- min(c(diff(base::range(private$mesh$loc[,1])), diff(base::range(private$mesh$loc[,2]))))
        size / rangeFactor
      }
      else range0 <- range
      kappa0 <- sqrt(8) / range0
      tau0 <- 1 / (sqrt(4 * pi) * kappa0 * sigma0)
      private$spde <- inla.spde2.matern(mesh=private$mesh,
                                        B.tau=cbind(log(tau0), -1, +1),
                                        B.kappa=cbind(log(kappa0), 0, -1),
                                        theta.prior.mean=c(0, 0),
                                        theta.prior.prec=c(0.1, 1))
      
      return(invisible(self))        
    },
        
    setCovariatesModel = function(covariatesModel, covariates) {
      if (missing(covariatesModel))
        covariatesModel <- ~ 1
      x <- if (missing(covariates)) terms(covariatesModel)
      else terms(covariatesModel, data=covariates)
      
      if (attr(x, "response") != 0)
        stop("The covariates model formula must be right-sided.")
      
      private$covariatesModel <- covariatesModel
      covariates <- colnames(getINLAModelMatrix(covariatesModel, covariates))
      intercept <- if (attr(x, "intercept")[1] == 0) NULL else "intercept"
      randomEffect <- "f(spatial, model=spde, group=spatial.group, control.group=list(model=\"ar1\"))"
      private$linearModel <- reformulate(termlabels=c(intercept, covariates, randomEffect), response="response", intercept=FALSE)
            
      return(invisible(self))
    },
    
    setSmoothingModel = function() {
      return(self$setCovariatesModel(~ 1))
    },
    
    addObservationStack = function(coordinates, time, response=NA, covariates, offset, offsetScale, tag="obs") {
      # TODO: allow defining link function
      
      if (missing(coordinates))
        stop("Required argument 'coordinates' missing.")
      if (missing(time))
        stop("Required argument 'time' missing.")
      
      if (is.null(private$mesh))
        stop("Mesh must be defined first.")
      if (is.null(private$spde))
        stop("Spatial prior must be defined first.")
      if (is.null(private$covariatesModel))
        stop("Covariates model must be defined first.")
      
      dataList <- list(response=response)
      if (!missing(offset)) {
        offsetScale <- if (missing(offsetScale)) SpaceTime::findScale(offset[1])
        else offsetScale
        dataList$E <- offset / private$offsetScale
      }
      
      coordinates <- as.matrix(coordinates) / self$getCoordinatesScale()
      modelMatrix <- getINLAModelMatrix(private$covariatesModel, covariates)
      nTime <- length(unique(time))
      timeIndex <- time - min(time) + 1
      fieldIndex <- inla.spde.make.index("spatial", n.spde=private$spde$n.spde, n.group=nTime)
      A <- inla.spde.make.A(private$mesh, loc=coordinates, group=timeIndex, n.group=nTime)
      
      effects <- if (private$hasIntercept()) list(c(fieldIndex, list(intercept=1))) else list(index)
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
            
      private$addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
        
    addPredictionStack = function(coordinates, time, response=NA, covariates, offset, offsetScale, tag="pred") {
      # TODO: finish this function
      
      effects <- if (private$hasIntercept()) list(c(fieldIndex, coordinates, list(intercept=1))) else list(c(index, coordinates))
      AList <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      addStack(data=dataList, A=AList, effects=effects, tag=tag)
      
      return(invisible(self))
    },
    
    setLikelihood = function(likelihood) {
      if (missing(likelihood))
        stop("Required argument 'likelihood' missing.")
      private$likelihood <- likelihood
      return(invisible(self))
    },
    
    estimate = function(verbose=F) {
      if (is.null(private$fullStack))
        stop("Data stack must be defined first.")
      
      dataStack <- inla.stack.data(private$fullStack, spde=private$spde)
      private$result <- try(inla(private$linearModel, family=private$likelihood, data=dataStack, E=dataStack$E,
                             control.predictor=list(A=inla.stack.A(private$fullStack), link=dataStack$link, compute=TRUE),
                             control.compute=list(waic=TRUE, config=TRUE),
                             verbose=verbose))
      if (inherits(private$result, "try-error") || private$result$ok == FALSE)
        stop("Estimation failed. Use verbose=TRUE to find the possible cause.")
      
      return(invisible(self))      
    },
    
    save = function(fileName) {
      save(distanceUnit, timeUnit, coordinatesScale, offsetScale, coordinates, mesh, spde, covariatesModel, linearModel, fullStack, likelihood, result, envir=private, file=fileName)
      return(invisible(self))
    },
    
    getSPDEResult = function() {
      if (is.null(private$result) || is.null(private$spde))
        stop("The model has not been estimated.")
      return(inla.spde2.result(private$result, "spatial", private$spde))
    },
    
    summarySpatialParameters = function() {
      spdeResult <- self$getSPDEResult()
      range <- SpaceTime::summaryINLAParameter(spdeResult$marginals.range.nominal[[1]], coordinatesScale=self$getCoordinatesScale())
      variance <- SpaceTime::summaryINLAParameter(spdeResult$marginals.variance.nominal[[1]])
      kappa <- SpaceTime::summaryINLAParameter(spdeResult$marginals.kappa[[1]], coordinatesScale=1/self$getCoordinatesScale())
      tau <- SpaceTime::summaryINLAParameter(spdeResult$marginals.tau[[1]])
      x <- rbind(kappa=kappa, tau=tau, range=range, variance=variance)
      colnames(x) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
      print(x)
      return(invisible(self))
    },
    
    summary = function() {
      if (is.null(private$result))
        stop("The model has not been estimated.")
      print(summary(private$result))
      return(invisible(self))
    },
    
    getFitted = function() {
      indexObserved <- inla.stack.index(private$fullStack, "obs")$data
      fitted <- private$result$summary.fitted.values$mean[indexObserved] * private$offsetScale
      return(fitted)      
    },
    
    summaryTemporalVariation = function() {
      observed <- private$response / private$offset * private$offsetScale
      fitted <- self$getFitted()
      x <- data.frame(time=private$time, observed=observed, fitted=fitted)
      x <- ddply(x, .(time), function(x) data.frame(observed=sum(x$observed), fitted=sum(x$fitted)))
      print(x)
      return(invisible(self))
    },
    
    plotTemporalVariation = function() {
      observed <- private$response / private$offset * private$offsetScale
      fitted <- self$getFitted()
      x <- data.frame(time=private$time, observed=observed, fitted=fitted)
      x <- ddply(x, .(time), function(x) data.frame(observed=sum(x$observed), fitted=sum(x$fitted)))
      x <- melt(x, id.vars="time", measure.vars=c("observed", "fitted"))
      p <- ggplot(x, aes(time, value, colour=variable)) + geom_line()
      print(p)
      return(invisible(self))
    },
    
    plotSpatialVariation = function(timeIndex) {
      
    }
  )
)
