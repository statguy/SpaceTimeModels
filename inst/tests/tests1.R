# Test the ContinuousSpaceDiscreteTimeModel class with Cameletti et al. (2012) data and model

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

mesh <- SpaceTime::SpatialMesh$new(knots=cbind(coordinates$UTMX, coordinates$UTMY))
mesh$construct(locDomain=borders, offset=c(10, 140), maxEdge=c(50, 1000), minAngle=c(26, 21), cutoff=0)
mesh$plot()

model <- SpaceTime::ContinuousSpaceDiscreteTimeModel$new(timeUnit="day")
model$setSpatialMesh(mesh)
model$setSpatialPrior()
model$setSmoothingModel()
model$getLinearModel()
model$addObservationStack(coordinates=coordinates_validation[Piemonte_data$Station.ID, c("UTMX","UTMY")], time=Piemonte_data$time, response=Piemonte_data$logPM10, tag="obs")
model$addObservationStack(coordinates=coordinates_validation[Piemonte_data_validation$Station.ID, c("UTMX","UTMY")], time=Piemonte_data_validation$time, tag="val")
model$estimate(verbose=T)
model$summary()
model$summarySpatial()
