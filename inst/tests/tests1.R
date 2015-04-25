# Test the DiscreteTimeContinuousSpaceModel class with Cameletti et al. (2012) data and model

Piemonte_data <- read.csv("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_data_byday.csv", header=TRUE, sep=",")
coordinates <- read.csv("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/coordinates.csv", header=TRUE, sep=",")
borders <- read.table("http://www.math.ntnu.no/inla/r-inla.org/case-studies/Cameletti2012/Piemonte_borders.csv", header=TRUE, sep=",")

n_stations <- length(coordinates$Station.ID)
n_data <- length(Piemonte_data$Station.ID)
n_days <- as.integer(n_data/n_stations)
Piemonte_data$logPM10 <- log(Piemonte_data$PM10)
Piemonte_data$time <- rep(1:n_days, each=n_stations)

model <- DiscreteTimeContinuousSpaceModel$new(timeUnit="day")
    model$constructMesh(coords=coordinates[Piemonte_data$Station.ID, c("UTMX","UTMY")], locDomain=borders, offset=c(10, 140), maxEdge=c(50, 1000), minAngle=c(26, 21), cutoff=0)
model$plotMesh()
model$setSpatialPrior()
model$setSmoothingModel()
model$getLinearModel()
model$addObservationStack(time=Piemonte_data$time, response=Piemonte_data$logPM10)
model$estimate(verbose=T)
model$summary()
model$summarySpatial()
