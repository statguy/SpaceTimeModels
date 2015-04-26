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
    coordsScale = 1,
    offsetScale = 1,
    
    coords = NULL,
    time = NULL,
    response = NULL,
    offset = 1,
    
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
    }
  ),
  public = list(
    getCoordinatesScale = function() return(private$coordsScale),
    getOffsetScale = function() return(private$offsetScale),
    
    getCoordinates = function() return(private$coords * private$coordsScale),
    getTime = function() return(private$time),
    getOffset = function() return(private$offset * private$offsetScale),
    getResponse = function() return(private$response),
    getCovariates = function() return(private$covariates),
    
    getMesh = function() return(private$mesh),
    getSPDEObject = function() return(private$spde),
    getModelMatrix = function() return(private$modelMatrix),
    getLinearModel = function() return(private$linearModel),
    getFullStack = function() return(private$fullStack),
    getResult = function() return(private$result),
 
    constructMesh = function(coords, coordsScale, cutoff=NULL, maxEdge=NULL, offset=NULL, minAngle=NULL, locDomain=NULL, convex) {
      if (missing(coords))
        stop("Required argument 'coords' missing.")
      if (missing(cutoff))
        stop("Required argument 'cutoff' missing.")
      if (missing(maxEdge))
        stop("Required argument 'maxEdge' missing.")
      
      private$coordsScale <- if (missing(coordsScale)) SpaceTime::findScale(coords[1,1])
      else coordsScale
      private$coords <- as.matrix(coords) / private$coordsScale
      meshCoords <- unique(private$coords)
      
      if (missing(convex)) {
        locDomain <- SpaceTime::nullScale(locDomain, private$coordsScale)
        private$mesh <- inla.mesh.2d(loc=private$coords, loc.domain=locDomain, cutoff=SpaceTime::nullScale(cutoff, private$coordsScale), max.edge=SpaceTime::nullScale(maxEdge, private$coordsScale), offset=SpaceTime::nullScale(offset, private$coordsScale), min.angle=minAngle)
      }
      else {
        boundary <- inla.nonconvex.hull(points=meshCoords, convex=convex)
        private$mesh <- inla.mesh.2d(boundary=boundary, cutoff=SpaceTime::nullScale(cutoff, private$coordsScale), max.edge=SpaceTime::nullScale(maxEdge, private$coordsScale), offset=SpaceTime::nullScale(offset, private$coordsScale), min.angle=minAngle)
      }
      
      return(invisible(self))
    },
    
    plotMesh = function() {
      if (is.null(private$mesh))
        stop("Mesh must be defined first.")
      
      plot(private$mesh)
      points(private$coords, pch='.', col='red')
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
    
    addObservationStack = function(time, response, covariates, offset, coordsScale, offsetScale, tag="obs") {
      # TODO: allow building stack for predictions
      # TODO: allow defining link function
      
      if (missing(time))
        stop("Required argument 'time' missing.")
      if (missing(response))
        stop("Required argument 'response' missing.")
      
      if (is.null(private$coords))
        stop("Coordinates must be defined first.")
      if (is.null(private$mesh))
        stop("Mesh must be defined first.")
      if (is.null(private$spde))
        stop("Spatial prior must be defined first.")
      if (is.null(private$covariatesModel))
        stop("Covariates model must be defined first.")
      
      if (!missing(offset)) {
        if (length(offset) == 1) offset <- rep(offset, length(response))
        offsetScale <- if (missing(offsetScale)) SpaceTime::findScale(offset[1])
        else offsetScale
        private$offset <- offset / private$offsetScale
      }
      
      private$time <- time
      private$response <- response
      modelMatrix <- getINLAModelMatrix(private$covariatesModel, covariates)
      
      nTime <- length(unique(private$time))
      timeIndex <- private$time - min(private$time) + 1
      index <- inla.spde.make.index("spatial", n.spde=private$spde$n.spde, n.group=nTime)
      A <- inla.spde.make.A(private$mesh, loc=private$coords, group=timeIndex, n.group=nTime)
      
      effects <- if (private$hasIntercept()) list(c(index, list(intercept=1)))
      else list(index)
      Alist <- if (!is.null(modelMatrix)) {
        effects[[2]] <- modelMatrix
        list(A, 1)
      }
      else list(A)
      
      obsStack <- inla.stack(data=list(response=private$response, E=private$offset, link=1),
                             A=Alist, effects=effects, tag=tag)
      private$fullStack <- if (is.null(private$fullStack)) inla.stack(obsStack)
      else {
        if (names(private$fullStack$data$index) == tag)
          stop("Stack with tag '", tag, "' already exists.")
        inla.stack(private$fullStack, obsStack)
      }
      
      return(invisible(self))
    },
    
    setLikelihoodModel = function(likelihood) {
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
      save(distanceUnit, timeUnit, coordsScale, offsetScale, coords, time, response, offset, mesh, spde, covariatesModel, linearModel, fullStack, likelihood, result, envir=private, file=fileName)
      return(invisible(self))
    },
    
    getSPDEResult = function() {
      if (is.null(private$result) || is.null(private$spde))
        stop("The model has not been estimated.")
      return(inla.spde2.result(private$result, "spatial", private$spde))
    },
    
    summarySpatialParameters = function() {
      spdeResult <- self$getSPDEResult()
      range <- SpaceTime::summaryINLAParameter(spdeResult$marginals.range.nominal[[1]], coordsScale=private$coordsScale)
      variance <- SpaceTime::summaryINLAParameter(spdeResult$marginals.variance.nominal[[1]])
      kappa <- SpaceTime::summaryINLAParameter(spdeResult$marginals.kappa[[1]], coordsScale=1/private$coordsScale)
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
