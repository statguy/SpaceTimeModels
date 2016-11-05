# Test the ContinuousSpaceDiscreteTimeModel class with Cameletti et al. (2012) data and model

library(SpaceTimeModels)

# Download data
Piemonte_data <- read.csv("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_data_byday.csv", header=TRUE, sep=",")
coordinates <- read.csv("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/coordinates.csv", header=TRUE, sep=",")
borders <- read.table("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_borders.csv", header=TRUE, sep=",")
Piemonte_data_validation <-read.table("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_data_byday_validation.csv", header=TRUE, sep=",")
coordinates_validation <-read.table("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/coordinates_validation.csv", header=TRUE, sep=",")
rownames(coordinates) = coordinates[,"Station.ID"]
rownames(coordinates_validation) = coordinates_validation[,"Station.ID"]

# Prepare data
n_stations <- length(coordinates$Station.ID)
n_stations_val <- length(coordinates_validation$Station.ID)
n_data <- length(Piemonte_data$Station.ID)
n_days <- as.integer(n_data/n_stations)
Piemonte_data$logPM10 <- log(Piemonte_data$PM10)
Piemonte_data$time <- rep(1:n_days, each=n_stations)
Piemonte_data_validation$logPM10 <- log(Piemonte_data_validation$PM10)
Piemonte_data_validation$time <- rep(1:n_days,each=n_stations_val)

mean_covariates <- apply(Piemonte_data[,3:10], 2, mean)
sd_covariates <- apply(Piemonte_data[,3:10], 2, sd)
Piemonte_data[,3:10] <- scale(Piemonte_data[,3:10], mean_covariates, sd_covariates)
Piemonte_data_validation[,3:10] <- scale(Piemonte_data_validation[,3:10], mean_covariates, sd_covariates)

# Put observations and validation points to a spatio-temporal data frame
obs <- spacetime::STIDF(sp::SpatialPoints(coordinates[Piemonte_data$Station.ID, c("UTMX","UTMY")]), as.Date(Piemonte_data$Date, "%d/%m/%y"), Piemonte_data)
val <- spacetime::STIDF(sp::SpatialPoints(coordinates_validation[Piemonte_data_validation$Station.ID, c("UTMX","UTMY")]), as.Date(Piemonte_data_validation$Date, "%d/%m/%y"), Piemonte_data_validation)

if (F) {
  # Take a subset of the data for quick testing
  obs <- obs[,obs@time["/2005-10-05"]]
  val <- val[,val@time["/2005-10-05"]]
  nrow(obs)
  nrow(val)
}

# Build estimation mesh
mesh <- SpaceTimeModels::SpatialMesh$new(knots=obs, locDomain=sp::SpatialPoints(borders), offset=c(10, 140), maxEdge=c(50, 1000), minAngle=c(26, 21), cutoff=0)
mesh$plot()

# Build model
formula <- ~ A + UTMX + UTMY + WS + TEMP + HMIX + PREC + EMI
model <- SpaceTimeModels::ContinuousSpaceDiscreteTimeModel $ 
  new() $
  setSpatialMesh(mesh) $
  setSpatialPrior() $
#  setSmoothingModel() $
  setCovariatesModel(formula, obs@data) $
  setLikelihood("gaussian") $
  setLinkFunction(gaussian()$link) $
  addObservationStack(sp = obs, response = obs@data$logPM10, covariates=obs@data) $
  addValidationStack(sp = val, covariates = val@data) $
  addPredictionStack(sp = obs)

# Print the linear model specification
model$getLinearModel()

# Estimate the model
model$estimate(verbose = T)

# Get summary of the estimated parameters
model$summary()
model$summarySpatialParameters()

# Print observed and fitted values in time
model$summaryTemporalVariation(timeIndex = time(obs))

# Plot temporal variation
model$plotTemporalVariation(timeIndex = time(obs))

# Quick plot estimates on a map
model$plotSpatialVariation(timeIndex = 1)
model$plotSpatialVariation(timeIndex = 2)

# Plot estimates on a map more neatly
rasters <- model$getSpatialVariationRaster(template = raster::extend(raster::extent(as.matrix(borders)), 10), width = 200, height = 200)
rasterVis::gplot(rasters$getLayer(1)) + ggplot2::geom_raster(aes(fill = value)) +
  ggplot2::geom_path(data = borders, aes(UTM_X, UTM_Y)) +
  ggplot2::geom_point(data = data.frame(obs@sp), aes(UTMX, UTMY)) +
  ggplot2::scale_fill_gradientn(colours = terrain.colors(40), guide = guide_legend(title = expression(PM[10]))) +
  ggplot2::coord_equal() +
  SpaceTimeModels::theme_raster() +
  ggplot2::theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 10))

# Validate the model
validation0 <- list(p=rep(NA, length(obs$logPM10)))
etaMean <- model$getFittedLinearPredictor()
etaSd <- model$getFittedLinearPredictor(variable = "sd")
validation0$res <- obs$logPM10 - etaMean
validation0$res.std <- validation0$res / sqrt(etaSd^2 + 1 / model$getFittedHyperparameters()[1, "mean"])
validation0$p <- pnorm(validation0$res.std)

validation <- list()
etaMean <- model$getFittedLinearPredictor(tag = "val")
etaSd <- model$getFittedLinearPredictor(variable = "sd", tag = "val")
validation$res <- val$logPM10 - etaMean
validation$res.std <- validation$res / sqrt(etaSd^2 + 1 / model$getFittedHyperparameters()[1, "mean"])
validation$p <- pnorm(validation$res.std)

validation$rmse <- sqrt(mean(validation$res^2, na.rm=TRUE))
validation$cor <- cor(val$logPM10, etaMean, use="pairwise.complete.obs", method="pearson")
validation$cover <- mean((validation$p > 0.025) & (validation$p < 0.975), na.rm = TRUE)

validation0
validation
