# SpaceTime - An interface for fitting parametric spatial and spatio-temporal models with R and INLA

TODO: 1D meshes, discrete spatial models, priors

## Introduction

The SpaceTime [R](http://www.r-project.org/) package provides a simplified interface for parametrizing
a few standard spatial and spatio-temporal models with [R-INLA](http://www.r-inla.org/). The models enable
parametric (ie. the parameters are estimated from the data) smoothing over space and time and determining
effect of covariates on response.

### Autocorrelation

Things close in space and time often resemble each other. We can observe this e.g. in nature where
large scale processes such as climate and geomorphic processes shape the biological processes. For example,
coniferous trees are found more in areas of cool climate and deciduous in warm climate.
Measurements obtained close to each other tend to explain each other better than distant ones.
This degree of similarity is known as autocorrelation, i.e. the correlation within the process itself.

From the statistical modelling point of view, dependency in data can cause estimates in a model with 
an independence assumption to be biased, e.g. providing too small p-values.
Several spatial and spatio-temporal models have been developed to take autocorrelation properly
into account.

Since there is often no repeated measurements available that would allow estimating the autocorrelation
structure directly, several model assumptions are made. These often include observations depending via
some neighbourhood structure (e.g. distance), stationarity (variance is constant across space/time) and
isotropy (independence of direction in space).

### Discrete and continuous data

There are two types of data for both time and space. Discrete data have
a clear separation such as years which can be indexed from $t=1, 2, \dots, T$ or regions $s=1, 2, \dots, S$
which are separated by boundaries. Continous data is given within a range such as duration
from $t=1$ to $t=2$, where there are infinitely many time points within this range. Similarly,
continuous space has infinitely many locations within the study area.
The package provides model classes for combinations of discrete and continuous data for
space and time.

### Models

The following model classes are currently implemented in the package:

* `DiscreteSpaceModel` for discrete spatial data
* `ContinuousSpaceModel` for continuous spatial data
* `ContinuousSpaceDiscreteTimeModel` for continuous spatial and discrete temporal data
* `ContinuousSpaceContinuousTimeModel` for continuous spatial and continuous temporal data

### Continuous models

Specifying continuous models require specifying an estimation mesh. The estimates are provided at
the mesh nodes and estimates at the observation or prediction locations are interpolated from
the node estimates. A higher number of nodes improves the estimation results, but increases
computational time.

The following mesh classes are currently implemented in the package:

* `TemporalMesh` for temporal data
* `SpatialMesh` for spatial data
* `NonConvexHullMesh` for spatial data for creating a non-convex hull around the observations

The continuous spatial models assume that dependencies between the mesh nodes are specified as
a function of distance with unknown scale and variance parameters to be estimated from the data.
The function is of Matérn class, described in more detail in
[the reference](http://www.math.ntnu.no/inla/r-inla.org/papers/jss/lindgren.pdf).

### Discrete models

Discrete models are specified with a neighbourhood structure. For time, it is often assumed
that the current time point t depends on the previous time point t-1 in some fashion.
For example, the autoregressive model is of the form $y_t = \phi y_{t-1} + \epsilon$,
where $\phi$ is the degree of dependency between subsequent observations and $\epsilon$
is the zero-centered Gaussian error term.

Regions are typically considered neighbours if they share the same border or are within
certain radius from the centre points.

Please refer to the R-INLA [documentation](http://www.r-inla.org/models/latent-models) for more
details of the autoregressive and the Besag models.

## Installation

Test version of R-INLA is required to be installed first, see [http://www.r-inla.org/download](here) for
the installation instructions. The SpaceTime package is installed with the `devtools` package using
the command `devtools::install_github("statguy/SpaceTime")`. Additional packages are installed
automatically from CRAN if needed. The package will be ready to use with the command
```
library(SpaceTime)
```

## Usage

The assembly line for constructing the models is the following:

1. Create model object.
2. Specify mesh for continuous spatial models / Specify the neighborhood structure for discrete spatial models.
3. Specify priors for space and time components. If omitted, default priors are used.
4. Specify covariates or intercept only (smoothing-only model).
5. Add 
  + Add observation locations, time points, observed responses and covariates
  + Add validation locations, time points and covariates
  + Add prediction locations, time points and covariates
6. Specify likelihood. If omitted, Gaussian likelihood is used as default.

Once the model is specified, the unknown parameters are ready to be estimated. Once the model is estimated,
the results can be extract from the model object.

### Model object

Objects are created with the `new()` method (constructor), for example, the model object
```
model <- SpaceTime::ContinuousSpaceDiscreteTimeModel$new()
```

### Coordinates

Due to numerical accuracy, spatial coordinates may need to be scaled down. For example, if the coordinates
are such that $x=6100000$ and $y=50000$, $1000000$ should be added to $y$ to "match" with $x$. Furthermore,
the coordinates should be scaled, e.g., by dividing by $1000000$ and thus obtaining $x=6.1$ and $y=1.05$,
which INLA handles better.

### Spatial mesh creation

Spatial mesh object is obtained, for example, with
```
mesh <- SpaceTime::NonConvexHullMesh$new(knots = cbind(x, y))
```
Coordinates in the variables `x` and `y` are provided for the constructor via the argument `knots`
that stores the coordinates in the mesh object for future use.

The mesh is specified with the `construct` method and depending on the type of the mesh, several
parameters are needed to be supplied for the method that specify topology of the mesh. The parameters
affect, for example, the number of nodes in the mesh.
The mesh is found by a triangulation over the study area and the process should be completed rather
quickly.
However, the triangulation gets stuck sometimes and the INLA triangulation process or R has to be
force-terminated (killed).
A set of suitable mesh parameters is often found through iteration. It is adviced to run the models
with s small number of mesh nodes first, especially with the spatio-temporal models, which may take
very long time to estimate with large meshes.

The mesh parameters are listed in the help (which is obtained with R command `?x`
where `x` is the class name, e.g. `?NonConvexHullMesh`) and explained in more detail in
[the reference](http://www.math.ntnu.no/inla/r-inla.org/papers/jss/lindgren.pdf), which the reader
is adviced to go through.

The command `mesh$plot()` plots the mesh with the observation locations.
Once created, the mesh is supplied for the model object with
````
model$setSpatialMesh(mesh)
```

### Priors

TODO

### Likelihood

Likelihood is specified with
```
model$setLikelihood("x")
```
where `x` is the likelihood. Please refer to the [R-INLA documentation](http://www.r-inla.org/models/likelihoods)
for the supported likelihoods.
When using count data likelihoods such as binomial and Poisson, the offset term is supplied
for the `addObservationStack` method via the `offset` argument. Otherwise the offset is assumed
to equal to 1. See the section Data stacks.

### Covariates and intercept-only model

Covariates are specified with the `setCovariatesModel` method, which takes arguments

* `covariatesModel` as for right-sided equation of the covariates to be included in the model.
* `covariates` as for the covariates data.

For example
```
model$setCovariatesModel(covariatesModel = ~ a + b, covariates = covariates)$
```

Intercept-only model can be specified with the `setSmoothingModel()` method,
which includes only intercept and spatial or spatio-temporal random effect in the model.
Such models provide only noise filtered (smoothed) observations.

To view the full specification of the linear part of the model, use the `getLinearModel()` method.
This will also print the random effect term as supplied to R-INLA.
The part for the covariates differs from the specified for the categorial variables (factors)
as they are expanded to consist of dummy variables.

### Data stacks

Observation, validation and prediction data is supplied in chunks that are tagged and stacked.
Such design allows indexing the data so that it can be later referenced by the tags.
The chunks consist of coordinates, time indices (for spatio-temporal models), covariates (optional)
and the observations for the observation chunk. The observations consist of responses and optionally
offsets for count data. For example
```
model$addObservationStack(coordinates=cbind(x.obs, y.obs), time=t.obs, response=response, offset=offset.obs, covariates=cov.obs, tag="observed")
model$addValidationStack(coordinates=cbind(x.val, y.val), time=t.val, covariates=cov.val, tag="validation")
model$addPredictionStack(coordinates=cbind(x.pred, y.pred), time=t.pred, covariates=cov.pred, tag="predicted")
```

### Estimation

Once the model is specified, estimation is started with the `estimate()` method.
Argument `verbose=TRUE` is recommended to be supplied to follow progress of the estimation.
Note that spatio-temporal models may take considerable amount of time and memory to be estimated.

### Extracting results

The following methods provide basic summaries of the estimated models

* `summary` for overall summary of the model.
* `summarySpatialParameters` for summary of the spatial parameters.

To access the R-INLA result object directly, `getResult()` method is provided.

## Examples

* [Cameletti et al. (2002)](http://www.r-inla.org/examples/case-studies/cameletti-et-al) reproduced with the SpaceTime package,
[link](https://github.com/statguy/SpaceTime/blob/master/inst/tests/tests1.R).

## References and supporting material

* [Lindgren, F., Rue, H. (2015). Bayesian Spatial Modelling with R-INLA. Journal of Statistical Software.](http://www.math.ntnu.no/inla/r-inla.org/papers/jss/lindgren.pdf)
* [Gelman, A., Hwang, J., Vehtari, A. (2013). Understanding predictive information criteria for Bayesian models. Statistics and Computing.](http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf)

## Contact

Jussi Jousimo, <jvj@iki.fi>
