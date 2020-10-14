library(testthat)
library(tidyverse)
library(dsmextra)

# Context of the set of tests

testthat::context("test outputs")

# load the Gulf of Mexico dolphin data
data(spermwhales)

# Extract the data
segs <- spermwhales$segs[1:50,]
predgrid <- spermwhales$predgrid[1:500,]

# Define relevant coordinate system
my_crs <- sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Define covariates of interest
my_cov <- c("Depth", "DistToCAS", "SST", "EKE", "NPP")

suppressWarnings(spermw.extra <- compute_extrapolation(samples = segs,
                                      covariate.names = my_cov,
                                      prediction.grid = predgrid,
                                      coordinate.system = my_crs))

testthat::test_that("Outputs are correct",{

  # compute_extrapolation returns a list containing data, rasters, and a summary if save.summary = TRUE

  testthat::expect_s3_class(spermw.extra, "extrapolation_results")
  testthat::expect_equal(length(spermw.extra), 8)
  testthat::expect_equal(names(spermw.extra), c("type", "data", "rasters", "summary", "covariate.names",
                                                "samples", "prediction.grid", "coordinate.system"))

  # Map_extrapolation returns an html / leaflet object

  suppressWarnings(map1 <- map_extrapolation(map.type = "extrapolation",
                  extrapolation.object = spermw.extra))

  testthat::expect_equal(object = class(map1), expected = c("leaflet", "htmlwidget"))

})
