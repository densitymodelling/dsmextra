#' Quantify extrapolation in multivariate environmental space
#'
#' Assesses univariate (Type I) and combinatorial (Type II) extrapolation in spatial ecological models such as density surface models of line transect data. Models are built in a reference (calibration) system and projected into one or more target (prediction) system(s). The function is based on original code from the \href{https://cran.r-project.org/web/packages/ecospat/index.html}{ecospat} package (Broennimann et al. 2016). Although the required inputs mirror those of the \href{https://cran.r-project.org/web/packages/dsm/index.html}{dsm} package (Miller et al. 2015), the function is not restricted to line/strip-transect data and can be applied to other survey types and predictive modelling scenarios. See the 'Examples' section and Bouchet et al. (2019) for more information.
#'
#' The function calculates values of the ExDet (EXtrapolation DETection) metric as originally proposed by Mesgaran et al. (2014). ExDet takes on strictly negative values during univariate extrapolation (i.e. when predictions are made outside the range of individual covariates), is strictly >1 during combinatorial extrapolation (i.e. when predictions are made within the range of individual covariates, but for combinations of environmental conditions not encountered in the sample), and lies within the range 0-1 when predictions are made in conditions analogous to those found in the reference system. The function also determines which covariates make the largest contribution to each type of extrapolation; this is the most influential covariate (MIC). See Mesgaran et al. (2014) for details.
#'
#' Note that \code{compute_extrapolation} returns results in both numerical and raster format. The latter is used to support mapping functions and requires the locations in \code{prediction.grid} to be evenly spaced. If this is not the case, \code{dsmextra} will attempt to automatically generate a raster with a resolution given by the \code{resolution} argument (and expressed in the units of \code{coordinate.system}). An error may be returned if no \code{resolution} is specified.
#'
#' The \code{data} list captures ExDet values at prediction locations (i.e. cells in \code{prediction.grid}) and is organised into multiple \code{data.frame} objects, as follows:
#'
#' \tabular{ll}{
#'   \code{all} \tab All prediction locations\cr
#'   \code{univariate} \tab Prediction locations subject to univariate extrapolation (only)\cr
#'   \code{combinatorial} \tab Prediction locations subject to combinatorial extrapolation (only)\cr
#'   \code{analogue} \tab Prediction locations where conditions are analogous to sampled conditions (only)\cr
#'  }
#'
#'  Each \code{data.frame} contains four columns:
#'
#'  \tabular{ll}{
#'   \code{ExDet} \tab ExDet values\cr
#'   \code{mic_univariate} \tab Integer identifying the univariate MIC\cr
#'   \code{mic_combinatorial} \tab Integer identifying the combinatorial MIC\cr
#'   \code{mic} \tab Integer identifying the MIC\cr
#'  }
#'
#' The \code{rasters} list comprises two elements, named \code{ExDet} and \code{mic}. Each contains individual rasters mapping ExDet and MIC values, respectively.
#'
#' @import purrr
#'
#' @param samples Sample (reference) dataset used for model building and calibration. This corresponds to the \code{segment.data} used when building density surface models in \code{dsm}. It must contain one column for each of the covariates in \code{covariate.names}.
#' @param covariate.names Character string. Names of the covariates of interest.
#' @param prediction.grid Prediction data.frame. This contains both geographic coordinates (\code{x}, \code{y}) and covariate values associated with the target locations for which predictions are desired. Typically, these locations are taken as the centroids of the grid cells in a spatial prediction grid/raster. See \code{\link[dsm]{predict.dsm}}.
#' @param coordinate.system Projected coordinate system relevant to the study location. Can be either a character string or an object of class \code{\link[sp]{CRS}}.
# @param print.summary Logical, defaults to \code{TRUE}. Outputs a summary of the results to the R console.
# @param print.precision Integer. Number of significant figures to be used when printing the summary. Default value of 2.
# @param save.summary Logical, defaults to \code{FALSE}. Adds summary statistics to the output list.
#' @param resolution Resolution of the output raster (in units relevant to \code{coordinate.system}). Only required if \code{prediction.grid} is irregular, and thus needs to be rasterised. Defaults to \code{NULL}.
#' @param verbose Logical. Show or hide possible warnings and messages.
#'
#' @return A list object containing extrapolation values in both \code{data.frame} and \code{\link[raster]{raster}} format. Also included are a summary object of class \code{extrapolation_results_summary} and a copy of function inputs (i.e, \code{coordinate.system}, \code{covariate.names}, and \code{prediction.grid}).
#'
#' @author Phil J. Bouchet
#'
#' @references Bouchet PJ, Miller DL, Roberts JJ, Mannocci L, Harris CM and Thomas L (2019). From here and now to there and then: Practical recommendations for extrapolating cetacean density surface models to novel conditions. CREEM Technical Report 2019-01, 59 p. \href{https://research-repository.st-andrews.ac.uk/handle/10023/18509}{https://research-repository.st-andrews.ac.uk/handle/10023/18509}
#'
#' Broennimann O, Di Cola V, Guisan A (2016). ecospat: Spatial Ecology Miscellaneous Methods. R package version 2.1.1. \href{https://CRAN.R-project.org/package=ecospat}{https://CRAN.R-project.org/package=ecospat}
#'
#' Mesgaran MB, Cousens RD, Webber BL (2014). Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity & Distributions, 20: 1147-115. DOI: \href{https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209}{10.1111/ddi.12209}
#'
#' Miller DL, Rexstad E, Burt L, Bravington MV, Hedley S (2015). dsm: Density Surface Modelling of Distance Sampling Data. R package version 2.2.9. \href{https://CRAN.R-project.org/package=dsm}{https://CRAN.R-project.org/package=dsm}
#' @export
#' @examples
#' library(dsmextra)
#'
#' #  --- EXAMPLE 1: Line transect surveys of sperm whales ---
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
#' # Assess extrapolation in the multivariate space defined by five covariates
#' spermw.extrapolation <- compute_extrapolation(samples = segs,
#'       covariate.names = my_cov,
#'       prediction.grid = predgrid,
#'       coordinate.system = my_crs)
#'
#' #  --- EXAMPLE 2: Crowdsourced records of Acacia cyclops (see Mesgaran et al. 2014)  ---
#'
#' library(dsmextra)
#' library(raster)
#' library(sp)
#' library(magrittr)
#'
#' # Native and introduced range
#' data(acacia)
#'
#' # Download BioClim data
#' r <- raster::getData(name = "worldclim", var = "bio", res = 10)
#'
#' # Define variables of interest
#' bioclim.variables <- c("bio1", "bio5", "bio6", "bio12", "bio13", "bio14")
#'
#' # Reference system (South Australia)
#' ref <- raster::crop(x = r[[bioclim.variables]], y = acacia$south_australia) %>%
#' raster::mask(x = ., mask = acacia$south_australia) %>%
#' raster::as.data.frame(x = ., xy = TRUE, na.rm = TRUE)
#'
#' # Target system (South Africa)
#' target <- raster::crop(x = r[[bioclim.variables]], y = acacia$south_africa) %>%
#' raster::mask(x = ., mask = acacia$south_africa) %>%
#' raster::as.data.frame(., xy = TRUE, na.rm = TRUE)
#'
#' # Assess extrapolation
#' bioclim.ex <- compute_extrapolation(samples = ref,
#'                                     covariate.names = bioclim.variables,
#'                                     prediction.grid = target,
#'                                     coordinate.system = sp::proj4string(r))
#'
#' # Make a map
#' map_extrapolation(map.type = "extrapolation", extrapolation.object = bioclim.ex)

compute_extrapolation <- function(samples,
                                  segments,
                                  covariate.names,
                                  prediction.grid,
                                  coordinate.system,
                                  resolution = NULL,
                                  verbose = TRUE){

  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------

  calls <- names(sapply(match.call(), deparse))[-1]

  if(any("segments" %in% calls)) {
    if(verbose) warning("The 'segments' argument is deprecated, please use 'samples' instead.")
    samples <- segments
  }

  if(!"data.frame"%in%class(prediction.grid)) stop("pred.grid must be of class data.frame")
  if(!"data.frame"%in%class(samples)) stop("samples must be of class data.frame")

  if(!"x"%in%names(prediction.grid) | !"y"%in%names(prediction.grid)) stop("pred.grid must contain x and y coordinates")

  if(!all(covariate.names%in%names(samples))) stop("Missing/unrecognised covariates in the sample data")
  if(!all(covariate.names%in%names(prediction.grid))) stop("Missing/unrecognised covariates in the prediction grid")

  coordinate.system <- check_crs(coordinate.system = coordinate.system)

  samples <- na.omit(samples) # Cannot have NA values
  prediction.grid <- na.omit(prediction.grid)

  #---------------------------------------------
  # Check if prediction grid is regular
  #---------------------------------------------

  check.grid <- prediction.grid %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(z = 1)

  grid.regular <- try(raster::rasterFromXYZ(check.grid), silent = TRUE)

  # grid.regular <- try(raster::rasterFromXYZ(check.grid,
  #                     res = ifelse(is.null(resolution), c(NA,NA), resolution)),
  #                     silent = TRUE)

  #---------------------------------------------
  # If grid is irregular, rasterise prediction.grid based on specified resolution
  #---------------------------------------------

  if (class(grid.regular) == "try-error") {
    if (is.null(resolution)) stop("Prediction grid cells are not regularly spaced.\nA target raster resolution must be specified. See package documentation for details.")

    if (verbose) warning("Prediction grid cells are not regularly spaced.\nData will be rasterised and covariate values averaged. See package documentation for details.")

    RasteriseGrid <- TRUE
  } else if (class(grid.regular) == "RasterLayer" & !is.null(resolution)) {
    if (verbose) warning("New resolution specified.\nData will be rasterised and covariate values averaged. See package documentation for details.")

    RasteriseGrid <- TRUE

  } else {
    RasteriseGrid <- FALSE
  }

  if(RasteriseGrid){

    check.grid$z <- NULL
    sp::coordinates(check.grid) <- ~x+y
    sp::proj4string(check.grid) <- coordinate.system

    # Create empty raster with desired resolution

    ras <- raster::raster(raster::extent(check.grid), res = resolution)
    raster::crs(ras) <- coordinate.system

    # Create individual rasters for each covariate

    ras.list <- purrr::map(.x = covariate.names,
                           .f = ~raster::rasterize(as.data.frame(check.grid), ras,
                                                   prediction.grid[,.x], fun = mean_ras)) %>%
      purrr::set_names(., covariate.names)

    # Combine all rasters

    ras.list <- raster::stack(ras.list)

    # Update prediction grid

    prediction.grid <- raster::as.data.frame(ras.list, xy = TRUE, na.rm = TRUE)


  } # End if

  if(verbose) message("Computing ...")

  #---------------------------------------------
  # Define reference and target systems
  #---------------------------------------------

  reference <- samples[, covariate.names]
  target <- prediction.grid[, covariate.names]

  #---------------------------------------------
  # Run the exdet tool from Mesgaran et al. (2014)
  #---------------------------------------------

  mesgaran <- ExDet(ref = reference,
                    tg = target,
                    xp = covariate.names)

  #---------------------------------------------
  # Add coordinates
  #---------------------------------------------

  mesgaran <- prediction.grid %>%
    dplyr::select(x,y) %>%
    cbind(., mesgaran)

  #---------------------------------------------
  # Return a list with univariate, combinatorial, and analog conditions as separate elements
  #---------------------------------------------

  reslist <- list(data=NULL, rasters=NULL)

  reslist$data$all <- mesgaran

  reslist$data$univariate <- mesgaran %>%
    dplyr::filter(., ExDet < 0)

  reslist$data$combinatorial <- mesgaran %>%
    dplyr::filter(., ExDet > 1)

  reslist$data$analogue <- mesgaran %>%
    dplyr::filter(., ExDet >= 0 & ExDet <= 1)

  #---------------------------------------------
  # Create rasters from extrapolation/MIC values
  #---------------------------------------------

  reslist$rasters$ExDet <- reslist$data %>%
    purrr::map(., ~ dplyr::select(., x, y, ExDet) %>%
                  safe_raster(.))%>%
    purrr::map(., "result")

  reslist$rasters$mic <- reslist$data %>%
    purrr::map(., ~ dplyr::select(., x, y, mic) %>%
                  safe_raster(.)) %>%
    purrr::map(., "result")

  #---------------------------------------------
  # Check that rasters have been produced for each extrapolation type
  #---------------------------------------------

  null.check <- purrr::map_lgl(.x = reslist$rasters$ExDet, .f = ~is.null(.x))

  ms <- names(null.check[null.check])
  ms <- purrr::map_dbl(.x = reslist$data[ms], .f = ~nrow(.x))
  # ms <- names(ms[ms>0])

  if(length(ms)>0){

    for(i in 1:length(ms)){

      if(ms[i]>0){

        # Extract data

        ds <- reslist$data[names(ms[i])]

        # Build raster

        predr <- prediction.grid %>%
          dplyr::select(x,y) %>%
          dplyr::mutate("ID" = 1) %>%
          raster::rasterFromXYZ(xyz = ., crs = coordinate.system)

        rs <- ps <- purrr::map(.x = ds,
                               .f= ~raster::rasterize(x = .x[,c("x", "y")], y = predr))

        # Reassign values

        for(i in 1:length(rs)){
          rs[[i]][!is.na(rs[[i]])] <- ds[[i]]$ExDet
          ps[[i]][!is.na(ps[[i]])] <- ds[[i]]$mic
        }

        reslist$rasters$ExDet <- append(reslist$rasters$ExDet, rs) %>%
          purrr::discard(is.null)

        reslist$rasters$mic <- append(reslist$rasters$mic, ps) %>%
          purrr::discard(is.null)


      } # End if ms[i] > 0

    } # End for loop length(ms)
  } # End if(length(ms)>0)

  #---------------------------------------------
  # Project rasters
  #---------------------------------------------

  for(r in 1:length(reslist$rasters$ExDet)){
    if(!is.null(reslist$rasters$ExDet[[r]]))raster::projection(reslist$rasters$ExDet[[r]]) <- coordinate.system}

  for(r in 1:length(reslist$rasters$mic)){
    if(!is.null(reslist$rasters$mic[[r]]))raster::projection(reslist$rasters$mic[[r]]) <- coordinate.system}

#  #---------------------------------------------
#  # Print/save summary
#  #---------------------------------------------

  sumres <- summarise_extrapolation(extrapolation.object = reslist,
                                    covariate.names = covariate.names,
                                    extrapolation = TRUE,
                                    mic = TRUE)

  class(sumres) <- c("extrapolation_results_summary", class(sumres))
  reslist <- append(x = reslist, values = list(summary = sumres))

  # Add function inputs to obviate need to specify them in map()
  reslist <- append(x = reslist, values = list(
                                  covariate.names = covariate.names,
                                  samples = samples,
                                  prediction.grid = prediction.grid,
                                  coordinate.system = coordinate.system))

  reslist <- append(list(type = c("extrapolation", "mic")), reslist)

  # Keep it classy
  class(reslist) <- c("extrapolation_results", class(reslist))

  if(verbose) message("Done!")
  return(reslist)

}
