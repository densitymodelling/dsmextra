library(testthat)
library(raster)
library(dsmextra)


par.tol <- 1e-5

# Context of the set of tests

testthat::context("Mid-Atlantic sperm whales")

# load the Gulf of Mexico dolphin data
data(spermwhales)

segs <- spermwhales$segs[1:50,]
predgrid <- spermwhales$predgrid[1:500,]

# Define relevant coordinate system
my_crs <- sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Define covariates of interest
my_cov <- c("Depth", "DistToCAS", "SST", "EKE", "NPP")


testthat::test_that("Do we get the same results?",{

  suppressWarnings(spermw.extra <- compute_extrapolation(samples = segs,
                        covariate.names = my_cov,
                        prediction.grid = predgrid,
                        coordinate.system = my_crs))

  testthat::expect_equal(spermw.extra$summary$extrapolation$univariate.n, 105, tolerance = par.tol)

  suppressWarnings(spermw.near <- compute_nearby(samples = segs,
                                covariate.names = my_cov,
                                prediction.grid = predgrid,
                                coordinate.system = my_crs,
                                nearby = 1))

  testthat::expect_equal(range(na.omit(raster::getValues(spermw.near))), c(0,54))

})
