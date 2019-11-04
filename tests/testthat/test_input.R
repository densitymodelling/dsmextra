library(testthat)
library(tidyverse)
library(janitor)
library(dsmextra)

# Context of the set of tests

testthat::context("test inputs")

# load the Gulf of Mexico dolphin data
data(spermwhales)

# Extract the data
segs <- spermwhales$segs
predgrid <- spermwhales$predgrid

# Define relevant coordinate system
my_crs <- sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Define covariates of interest
my_cov <- c("Depth", "DistToCAS", "SST", "EKE", "NPP")

testthat::test_that("Missing columns",{

  # Missing covariates

  testthat::expect_error(object = compute_extrapolation(segments = segs,
                        covariate.names = my_cov,
                        prediction.grid = predgrid[, 1:3],
                        coordinate.system = my_crs),
                        regexp = "Missing covariates in the prediction grid")

  testthat::expect_error(object = compute_extrapolation(segments = segs[,1:4],
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid,
                                                        coordinate.system = my_crs),
                         regexp = "Missing covariates in the segment data")

  # Missing coordinates

  testthat::expect_error(object = compute_extrapolation(segments = segs,
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid[, !names(predgrid)%in%c("x")] ,
                                                        coordinate.system = my_crs),
                         regexp = "pred.grid must contain x and y coordinates")

})

testthat::test_that("Wrong inputs",{

  # Segment or prediction datasets provided as a spatial objects ----------------------------------

  segs.pts <- sp::SpatialPointsDataFrame(coords = segs[, c("x", "y")], data = segs, proj4string = my_crs)
  pred.pts <- sp::SpatialPointsDataFrame(coords = predgrid[, c("x", "y")], data = predgrid, proj4string = my_crs)

  testthat::expect_error(object = compute_extrapolation(segments = segs.pts,
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid,
                                                        coordinate.system = my_crs),
                         regexp = "no method for coercing this S4 class to a vector")

  # Covariate names don't match ----------------------------------

  segs.wrongnames <- janitor::clean_names(segs.wrongnames)
  pred.wrongnames <- janitor::clean_names(predgrid)

  testthat::expect_error(object = compute_extrapolation(segments = segs.wrongnames,
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid,
                                                        coordinate.system = my_crs),
                         regexp = "Missing/unrecognised covariates in the segment data")

  testthat::expect_error(object = compute_extrapolation(segments = segs,
                                                        covariate.names = my_cov,
                                                        prediction.grid = pred.wrongnames,
                                                        coordinate.system = my_crs),
                         regexp = "Missing/unrecognised covariates in the prediction grid")

  # Unrecognised coordinate system ----------------------------------

  testthat::expect_error(object = compute_extrapolation(segments = segs,
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid,
                                                        coordinate.system = "batman"),
                         regexp = "Unrecognised coordinate system")

  # Missing coordinate system ----------------------------------

  testthat::expect_error(object = compute_extrapolation(segments = segs,
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid),
                         regexp = "argument \"coordinate.system\" is missing, with no default")

  # Invalid inputs ----------------------------------

  # compare_covariates

  testthat::expect_error(object = compare_covariates(extrapolation.type = "bananas",
                                                     segments = segs,
                                                     covariate.names = my_cov,
                                                     n.covariates = 3,
                                                     prediction.grid = predgrid,
                                                     coordinate.system = my_crs),
                         regexp = "Unknown extrapolation type")

  testthat::expect_error(object = compare_covariates(extrapolation.type = "both",
                                                     segments = segs,
                                                     covariate.names = my_cov,
                                                     n.covariates = 15,
                                                     prediction.grid = predgrid,
                                                     coordinate.system = my_crs),
                         regexp = "n.covariates exceeds the number of covariates available")

  # compute_nearby

  testthat::expect_error(object = compute_nearby(segments = segs,
                                                 prediction.grid = predgrid,
                                                 coordinate.system = my_crs,
                                                 covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
                                                 nearby = 0),
                         regexp = "nearby must be strictly positive")

  testthat::expect_error(object = compute_nearby(segments = segs,
                                                 prediction.grid = predgrid,
                                                 coordinate.system = my_crs,
                                                 covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
                                                 nearby = 1,
                                                 max.size = "big data"),
                         regexp = "Non-numeric input to argument: max.size")

  testthat::expect_error(object = compute_nearby(segments = segs,
                                                 prediction.grid = predgrid,
                                                 coordinate.system = my_crs,
                                                 covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
                                                 nearby = 1, no.partitions = 1e48),
                         regexp = "Number of partitions too large")

  # map_extrapolation

  extrapolation.calc <- compute_extrapolation(segments = segs,
                                              print.summary = FALSE,
                                              covariate.names = my_cov,
                                              prediction.grid = predgrid,
                                              coordinate.system = my_crs)

  testthat::expect_error(object = map_extrapolation(map.type = "Humpback whale",
                                                    extrapolation.values = extrapolation.calc,
                                                    covariate.names = my_cov,
                                                    prediction.grid = predgrid,
                                                    coordinate.system = my_crs),
                         regexp = "Unknown map type")


  testthat::expect_error(object = map_extrapolation(map.type = NULL,
                                                    extrapolation.values = extrapolation.calc,
                                                    covariate.names = my_cov,
                                                    prediction.grid = predgrid,
                                                    coordinate.system = my_crs),
                         regexp = "Argument 'maptype' must be specified")

  testthat::expect_error(object = map_extrapolation(map.type = "nearby",
                                                    extrapolation.values = extrapolation.calc,
                                                    covariate.names = my_cov,
                                                    prediction.grid = predgrid,
                                                    coordinate.system = my_crs),
                         regexp = "Argument 'gower.values' cannot be NULL when maptype is set to 'nearby'")

})
