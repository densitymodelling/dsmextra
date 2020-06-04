#' Extrapolation detection in density surface models
#'
#' The \code{dsmextra} package provides functions for detecting, diagnosing and visualising extrapolation in multivariate environmental space, with applications to density surface models (DSMs) created in the \code{\link[dsm]{dsm}} package.
#'
#' Further information on distance sampling methods and example code is available at \href{http://distancesampling.org/R/}{http://distancesampling.org/R/}.
#'
#' For help with distance sampling, there is a Google Group \href{https://groups.google.com/forum/#!forum/distance-sampling}{https://groups.google.com/forum/#!forum/distance-sampling}.
#'
#' @author Phil J. Bouchet, Laura Mannocci, David Miller
#' @docType package
#' @name dsmextra
#' @importFrom dplyr %>%
#' @importFrom purrr %||%
NULL

# Quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Sperm whale sightings in the Mid-Atlantic
#'
#' Data from a combination of several NOAA visual surveys conducted for cetaceans in the U.S. Mid-Atlantic. Survey tracklines were split into segments for analysis in the \code{\link[dsm]{dsm}} R package. Spatial coordinates (\code{x}, \code{y}) and values of both static and dynamic environmental covariates are availble for each segment mid-point and for the centroids of a 10 km grid overlaid on the study area. Covariates include: depth, sea surface temperature (\code{SST}), net primary production (\code{NPP}), distance to nearest canyon or seamount (\code{DistToCAS}), and eddy kinetic energy (\code{EKE}).
#'
#' An example extrapolation assessment is provided at \url{https://densitymodelling.github.io/dsmextra/articles/dsmextra-vignette.html}.
#'
#' @references NOAA Northeast Fisheries Science Center (2004). A survey for abundance and distribution of cetaceans in the U.S. Mid-Atlantic with an emphasis on pilot whales. Report of a survey conducted aboard NOAA Ship Gordon Gunter Cruise GU-04-03 (028).
#'
#' NOAA Northeast Fisheries Science Center (2004). Report on the 2004 Mid-Atlantic Marine Mammal Shipboard Abundance Survey aboard the R/V Endeavor, Cruise No. EN 04-395/396.
#'
#' Palka, DL 2012. Cetacean abundance estimates in US northwestern Atlantic Ocean waters from summer 2011 line transect survey. Northeast Fisheries Science Center Reference Document 12-29, 37 p.
#'
#' Palka, DL 2006 Summer abundance estimates of cetaceans in US North Atlantic Navy operating areas. Northeast Fisheries Science Center Reference Document 06-03, 41 p.
#'
#' @format A list of two
#' \describe{
#'   \item{segs}{segment-level data, as used in the \code{\link[dsm]{dsm}} package}
#'   \item{predgrid}{prediction data.frame}
#' }
#'
#' @name spermwhales
#' @docType data
#' @source Data provided by Debi Palka (NOAA North East Fisheries Science Center) and Lance Garrison (NOAA South East Fisheries Science Center). Initial data processing by Jason Roberts (Marine Geospatial Ecology Lab, Duke University). A subset of the data are held and described on OBIS-SEAMAP at: \url{http://seamap.env.duke.edu/dataset/396}.
#' @keywords datasets
NULL
