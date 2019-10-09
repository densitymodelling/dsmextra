#' Full extrapolation assessment
#'
#' Performs a complete evaluation of both univariate (Type I) and combinatorial (Type II) extrapolation in density surface models of line transect data, by calling relevant functions from \code{dsmextra}. As such, arguments \code{extrapolation_analysis} mirror those of the individual functions from which they are taken:
#'  \tabular{ll}{
#'   \code{summarise/summary.arguments} \tab Arguments from \code{\link{summarise_extrapolation}} \cr
#'   \code{compare.arguments} \tab Arguments from \code{\link{compare_covariates}} \cr
#'   \code{nearby.arguments} \tab Arguments from \code{\link{compute_nearby}} \cr
#'   \code{map.arguments} \tab Arguments from \code{\link{map_extrapolation}} \cr
#'  }
#'
#' @inheritParams compute_extrapolation
#' @param summarise.extrapolation Logical. If TRUE, run \code{\link{summarise_extrapolation}}.
#' @param summary.print.precision Integer. Number of significant figures to be used when printing the extrapolation summary. Default value of 2.
#'
#' @param compare.covariates Logical. If TRUE, run \code{\link{compare_covariates}}.
#' @param compare.extrapolation.type Character string indicating the type of extrapolation to be assessed. One of \code{univariate}, \code{combinatorial}, or \code{both} (default).
#' @param compare.n.covariates Integer. Maximum number of covariates. The function will compare all combinations of 1 to \code{n.covariates} covariates.
#' @param compare.create.plots Logical, defaults to \code{FALSE}. Whether to produce summary plots.
#' @param compare.display.percent Logical, defaults to \code{TRUE}. Scales the y-axis of the summary plots as a percentage of the total number of grid cells in \code{prediction.grid}.
#'
#' @param nearby.compute Logical. If TRUE, run \code{\link{compute_nearby}}.
#' @param nearby.nearby Scalar indicating which reference data points are considered to be 'nearby' (i.e. withing ‘nearby’ mean geometric Gower's distances of) prediction points. Defaults to 1.
#' @param nearby.max.size Minimum size threshold for partitioning computations. Calculated as \code{\link[base]{prod}(\link[base]{nrow}(segments),\link[base]{nrow}(prediction.grid))}. Has a default value of \code{1e7}.
#' @param nearby.no.partitions Integer. Number of desired partitions of the data (default of 10).
#'
#' @param map.generate Logical. If TRUE, run \code{\link{map_extrapolation}}.
#' @param map.sightings Species observations (optional). Can be supplied as a \code{matrix} of coordinates, a \code{data.frame}, a \code{\link[sp]{SpatialPoints}} object or a \code{\link[sp]{SpatialPointsDataFrame}} object. Circle markers will be proportional to group size if the data contain a column labelled \code{size}.
#' @param map.tracks Survey tracks (optional). Can be supplied as a \code{matrix} of coordinates, a \code{data.frame}, a \code{\link[sp]{SpatialLines}} object or a \code{\link[sp]{SpatialLinesDataFrame}} object. A \code{TransectID} field is required for matrix or data.frame inputs.
#'
#' @export
#' @author Phil J. Bouchet
#' @references Bouchet PJ, Miller DL, Roberts JJ, Mannocci L, Harris CM and Thomas L (2019). From here and now to there and then: Practical recommendations for extrapolating cetacean density surface models to novel conditions. CREEM Technical Report 2019-01, 59 p. \href{https://research-repository.st-andrews.ac.uk/handle/10023/18509}{https://research-repository.st-andrews.ac.uk/handle/10023/18509}
#'
#' Mesgaran MB, Cousens RD, Webber BL (2014). Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity & Distributions, 20: 1147-1159, DOI: \href{https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209}{10.1111/ddi.12209}
#'
#' @examples
#' library(dsmextra)
#'
#' # Load the Mid-Atlantic sperm whale data (see ?spermwhales)
#' data(spermwhales)
#'
#' # Extract the data
#' segs <- spermwhales$segs
#' predgrid <- spermwhales$predgrid
#'
#' # Define relevant coordinate system
#' my_crs <- sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0
#'  +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
#'
#' # Define covariates of interest
#' my_cov <- c("Depth", "DistToCAS", "SST", "EKE", "NPP")
#'
#' spermw.analysis <- extrapolation_analysis(segments = segs,
#'                                          covariate.names = my_cov,
#'                                          prediction.grid = predgrid,
#'                                          coordinate.system = my_crs,
#'                                          summarise.extrapolation = TRUE,
#'                                          summary.print.precision = 2,
#'                                          compare.covariates = TRUE,
#'                                          compare.extrapolation.type = "both",
#'                                          compare.n.covariates = NULL,
#'                                          compare.create.plots = TRUE,
#'                                          compare.display.percent = TRUE,
#'                                          nearby.compute = TRUE,
#'                                          nearby.nearby = 1,
#'                                          nearby.max.size = 1e7,
#'                                          nearby.no.partitions = 10,
#'                                          map.generate = TRUE)
extrapolation_analysis <- function(segments,
                                  covariate.names,
                                  prediction.grid,
                                  coordinate.system,
                                  summarise.extrapolation = TRUE,
                                  summary.print.precision = 2,
                                  compare.covariates = FALSE,
                                  compare.extrapolation.type = "both",
                                  compare.n.covariates = NULL,
                                  compare.create.plots = FALSE,
                                  compare.display.percent = TRUE,
                                  nearby.compute = TRUE,
                                  nearby.nearby = 1,
                                  nearby.max.size = 1e7,
                                  nearby.no.partitions = 10,
                                  map.generate = TRUE,
                                  map.sightings = NULL,
                                  map.tracks = NULL){

  #---------------------------------------------
  # Initiate list (in which results will be stored)
  #---------------------------------------------

  resl <- list()

  #---------------------------------------------
  # Assess extrapolation
  #---------------------------------------------

  message("=== Assessing extrapolation ===\n")

  if(summarise.extrapolation){

    ex1 <- compute_extrapolation(segments = segments,
                                 covariate.names = covariate.names,
                                 prediction.grid = prediction.grid,
                                 coordinate.system = coordinate.system,
                                 print.summary = TRUE,
                                 print.precision = summary.print.precision,
                                 save.summary = TRUE)

  }else{

    ex1 <- compute_extrapolation(segments = segments,
                                 covariate.names = covariate.names,
                                 prediction.grid = prediction.grid,
                                 coordinate.system = coordinate.system,
                                 print.summary = FALSE,
                                 print.precision = summary.print.precision,
                                 save.summary = FALSE)

  }

  resl$extrapolation <- ex1 # Store results

  #---------------------------------------------
  # Compare covariate combinations
  #---------------------------------------------

  message("\n=== Testing covariate combinations ===\n")

  if(compare.covariates){

    compare_covariates(extrapolation.type = compare.extrapolation.type,
                       covariate.names = covariate.names,
                       n.covariates = compare.n.covariates,
                       segments = segments,
                       prediction.grid = prediction.grid,
                       coordinate.system = coordinate.system,
                       create.plots = compare.create.plots,
                       display.percent = compare.display.percent)

  }

  #---------------------------------------------
  # Neighbourhood metrics
  #---------------------------------------------

  message("\n=== Calculating neighbourhood metrics ===\n")

  if(nearby.compute){

    ex2 <- compute_nearby(segments = segments,
                          prediction.grid = prediction.grid,
                          coordinate.system = coordinate.system,
                          covariate.names = covariate.names,
                          nearby = nearby.nearby,
                          max.size = nearby.max.size,
                          no.partitions = nearby.no.partitions)

    resl$nearby <- ex2 # Store results

  }

  #---------------------------------------------
  # Generate maps
  #---------------------------------------------

  message("\n=== Generating maps ===\n")

  if(map.generate){

    m1 <- map_extrapolation(map.type = "extrapolation",
                            extrapolation.values = ex1,
                            covariate.names = covariate.names,
                            prediction.grid = prediction.grid,
                            coordinate.system = coordinate.system,
                            sightings = map.sightings,
                            tracks = map.tracks)

    print(m1)

    m2 <- map_extrapolation(map.type = "mic",
                            extrapolation.values = ex1,
                            covariate.names = covariate.names,
                            prediction.grid = prediction.grid,
                            coordinate.system = coordinate.system,
                            sightings = map.sightings,
                            tracks = map.tracks)

    print(m2)

    resl$maps <- list(extrapolation = m1, mic = m2)

    if(nearby.compute){

      m3 <- map_extrapolation(map.type = "nearby",
                              gower.values = ex2,
                              covariate.names = covariate.names,
                              prediction.grid = prediction.grid,
                              coordinate.system = coordinate.system,
                              sightings = map.sightings,
                              tracks = map.tracks)}

    resl$maps <- list(extrapolation = m1, mic = m2, nearby = m3) # Store results

    print(m3)

  }

  message("\n=== Analysis complete! ===\n")
  return(resl)

}
