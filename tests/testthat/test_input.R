library(testthat)
library(tidyverse)
library(janitor)
library(dsmextra)

par.tol <- 1e-5

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

  # Segment or prediction datasets provided as a spatial objects

  segs.pts <- sp::SpatialPointsDataFrame(coords = segs[, c("x", "y")], data = segs, proj4string = my_crs)
  pred.pts <- sp::SpatialPointsDataFrame(coords = predgrid[, c("x", "y")], data = predgrid, proj4string = my_crs)

  testthat::expect_error(object = compute_extrapolation(segments = segs.pts,
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid,
                                                        coordinate.system = my_crs),
                         regexp = "no method for coercing this S4 class to a vector")

  # Covariate names don't match

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

  # Unrecognised coordinate system

  testthat::expect_error(object = compute_extrapolation(segments = segs,
                                                        covariate.names = my_cov,
                                                        prediction.grid = predgrid,
                                                        coordinate.system = "batman"),
                         regexp = "Unrecognised coordinate system")



})
