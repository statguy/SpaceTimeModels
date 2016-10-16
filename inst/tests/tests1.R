# Test the ContinuousSpaceDiscreteTimeModel class with Cameletti et al. (2012) data and model

library(SpaceTimeModels)

Piemonte_data <- read.csv("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_data_byday.csv", header=TRUE, sep=",")
coordinates <- read.csv("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/coordinates.csv", header=TRUE, sep=",")
borders <- read.table("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_borders.csv", header=TRUE, sep=",")
Piemonte_data_validation <-read.table("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_data_byday_validation.csv", header=TRUE, sep=",")
coordinates_validation <-read.table("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/coordinates_validation.csv", header=TRUE, sep=",")
rownames(coordinates) = coordinates[,"Station.ID"]
rownames(coordinates_validation) = coordinates_validation[,"Station.ID"]

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

obs <- spacetime::STIDF(sp::SpatialPoints(coordinates[Piemonte_data$Station.ID, c("UTMX","UTMY")]), as.Date(Piemonte_data$Date, "%d/%m/%y"), Piemonte_data)
val <- spacetime::STIDF(sp::SpatialPoints(coordinates_validation[Piemonte_data_validation$Station.ID, c("UTMX","UTMY")]), as.Date(Piemonte_data_validation$Date, "%d/%m/%y"), Piemonte_data_validation)

# Take a subset of the data for quick testing
if (F) {
  obs <- obs[,obs@time["/2005-10-05"]]
  val <- val[,val@time["/2005-10-05"]]
  nrow(obs)
  nrow(val)
}

mesh <- SpaceTimeModels::SpatialMesh$new(knots=obs, locDomain=sp::SpatialPoints(borders), offset=c(10, 140), maxEdge=c(50, 1000), minAngle=c(26, 21), cutoff=0)
mesh$plot()

formula <- ~ A + UTMX + UTMY + WS + TEMP + HMIX + PREC + EMI

model <- SpaceTimeModels::ContinuousSpaceDiscreteTimeModel$new()
model$setSpatialMesh(mesh)
model$setSpatialPrior()
#model$setSmoothingModel()
model$setCovariatesModel(formula, obs@data)
model$getLinearModel()
#model$setLikelihood("gaussian")
#model$setLinkFunction(gaussian()$link)
model$addObservationStack(sp=obs, response=obs@data$logPM10, covariates=obs@data)
model$addValidationStack(sp=val, covariates=val@data)
model$addPredictionStack(sp=obs)
model$estimate(verbose=T)
model$summary()
model$summarySpatialParameters()
model$summaryTemporalVariation(timeIndex=time(obs))
model$plotTemporalVariation(timeIndex=time(obs))
model$plotSpatialVariation(timeIndex=1)
model$plotSpatialVariation(timeIndex=2)

validation0 <- list(p=rep(NA, length(obs$logPM10)))
etaMean <- model$getFittedLinearPredictor()
etaSd <- model$getFittedLinearPredictor(variable = "sd")
validation0$res <- obs$logPM10 - etaMean
validation0$res.std <- validation0$res / sqrt(etaSd^2 + 1/model$getFittedHyperparameters()[1, "mean"])
validation0$p <- pnorm(validation0$res.std)

validation <- list()
etaMean <- model$getFittedLinearPredictor(tag = "val")
etaSd <- model$getFittedLinearPredictor(variable = "sd", tag = "val")
validation$res <- val$logPM10 - etaMean
validation$res.std <- validation$res / sqrt(etaSd^2 + 1/model$getFittedHyperparameters()[1, "mean"])
validation$p <- pnorm(validation$res.std)

validation$rmse <- sqrt(mean(validation$res^2, na.rm=TRUE))
validation$cor <- cor(val$logPM10, etaMean, use="pairwise.complete.obs", method="pearson")
validation$cover <- mean((validation$p > 0.025) & (validation$p < 0.975), na.rm = TRUE)

validation0
validation
