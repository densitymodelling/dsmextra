#' Quantify extrapolation in multivariate environmental space
#'
#' Assesses univariate (Type I) and combinatorial (Type II) extrapolation in density surface models of line transect data. Models are built in a reference (calibration) system and projected into one or more target (prediction) systems. The function is based on original code from the \code{\href{https://cran.r-project.org/web/packages/ecospat/index.html}{ecospat}} package, and tailored to outputs returned by the \code{\href{https://cran.r-project.org/web/packages/dsm/index.html}{dsm}} package, although it can be applied to other predictive modelling scenarios. See Bouchet et al. (2019) and Mesgaran et al. (2014) for details.
#'
#' The function calculates values of the ExDet (EXtrapolation DETection) metric as originally proposed by Mesgaran et al. (2014). ExDet is strictly negative during univariate extrapolation (i.e. when predictions are made outside the range of individual covariates), strictly positive and >1 during combinatorial extrapolation (i.e. when predictions are made within the range of individual covariates, but for combinations of envionmental conditions not encountered in the sample), and lies within the range 0-1 when predictions are made in conditions analogous to those found in the reference system. The function also determines which covariate makes the largest contribution to each type of extrapolation; this is the most influential covariate (MIC). See Mesgaran et al. (2014) for details.
#'
#' The `data` list captures ExDet values in prediction locations (i.e. cells in prediction.grid) and is organised into multiple data.frame objects, as follows:
#'
#' \tabular{ll}{
#'   \code{all} \tab All prediction locations\cr
#'   \code{univariate} \tab Prediction locations subject to univariate extrapolation (only)\cr
#'   \code{combinatorial} \tab Prediction locations subject to combinatorial extrapolation (only)\cr
#'   \code{analogue} \tab Prediction locations where conditions are analogous to sampled conditions (only)\cr
#'  }
#'
#'  Each data.frame contains four columns:
#'
#'  \tabular{ll}{
#'   \code{ExDet} \tab ExDet values\cr
#'   \code{mic_univariate} \tab Integer identifying the univariate MIC\cr
#'   \code{mic_combinatorial} \tab Integer identifying the combinatorial MIC\cr
#'   \code{mic} \tab Integer identifying the MIC()\cr
#'  }
#'
#' The `rasters` list consists of two elements, named `ExDet` and `mic`. Each contains individual rasters mapping ExDet values and MICs, respectively.
#'
#' @param segments Segment data.frame (i.e. surveyed transects divivided into segments for analysis). This is the reference dataset used for model building and calibration. This must contain one column for each of the covariates in `covariate.names`.
#' @param covariate.names Character string. Names of the covariates of interest.
#' @param prediction.grid Prediction data.frame. This contains both the geographic coordinates (`x`, `y`) and the covariate values associated with the target locations for which predictions are desired. Typically, these locations are taken as the centroids of the grid cells in a spatial prediction grid/raster.
#' @param coordinate.system Projected coordinate system relevant to the study location. Can be either a character string or an object of class \code{\link[sp]{CRS}}.
#' @param print.summary Logical, defaults to TRUE. Outputs a summary of the results to the console
#' @param print.precision Integer. Number of significant figures to be used when printing the summary. Default value of 2.
#' @param save.summary Logical. Adds summary statistics to the output list. Defaults to FALSE.
#' @param resolution Resolution of the output raster (in units relevant to `coordinate.system`). Only required if prediction.grid is irregular, and thus needs to be rasterised. Defaults to NULL.
#' @return A list object with two elements, namely values of the extrapolation metrics and their spatial representation in raster form
#' @author Phil J. Bouchet
#' @references Bouchet PJ, Miller DL, Roberts JJ, Mannocci L, Harris CM and Thomas L (2019). From here and now to there and then: Practical recommendations for extrapolating cetacean density surface models to novel conditions. CREEM Technical Report 2019-01, 59 p. \href{https://research-repository.st-andrews.ac.uk/handle/10023/18509}{https://research-repository.st-andrews.ac.uk/handle/10023/18509}
#'
#' Mesgaran MB, Cousens RD, Webber BL (2014). Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity & Distributions, 20: 1147-1159, DOI: \href{https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209}{10.1111/ddi.12209}
#' @export
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
#' # Assess extrapolation in the multivariate space defined by five covariates
#'spermw.extrapolation <- compute_extrapolation(segments = segs,
#'       covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
#'       prediction.grid = predgrid,
#'       coordinate.system = my_crs,
#'       print.summary = TRUE,
#'       save.summary = TRUE,
#'       print.precision = 2)
compute_extrapolation <- function(segments,
                                  covariate.names,
                                  prediction.grid,
                                  coordinate.system,
                                  print.summary = TRUE,
                                  print.precision = 2,
                                  save.summary = FALSE,
                                  resolution = NULL){

  #'---------------------------------------------
  # Prints out progress bar if the function is
  # called within a call to compare_covariates
  #'---------------------------------------------

  if(exists("call.compare")){
    if(exists("pb")){pb$tick()$print()}
  }

  #'---------------------------------------------
  # Perform function checks
  #'---------------------------------------------

  if(!class(prediction.grid)=="data.frame") stop("pred.grid must be of class data.frame")

  if(!"x"%in%names(prediction.grid) | !"y"%in%names(prediction.grid)) stop("pred.grid must contain x and y coordinates")

  if(!all(covariate.names%in%names(segments))) stop("Missing covariates in the segment data")
  if(!all(covariate.names%in%names(prediction.grid))) stop("Missing covariates in the prediction grid")

  if(!class(coordinate.system)=="CRS"){

    coord.err <- tryCatch(expr = sp::CRS(coordinate.system),
                          error = function(e) return(NA))

    if(is.na(coord.err)){stop('Unrecognised coordinate system')
    }else{coordinate.system <- sp::CRS(coordinate.system)}
  }


  #'---------------------------------------------
  # Check if prediction grid is regular
  #'---------------------------------------------

  check.grid <- prediction.grid %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(z = 1)

  grid.regular <- try(raster::rasterFromXYZ(check.grid), silent = TRUE)

  #'---------------------------------------------
  # If grid is irregular, rasterise prediction.grid based on specified resolution
  #'---------------------------------------------

  if(class(grid.regular)=="try-error"){

    if(is.null(resolution)) stop('Prediction grid cells are not regularly spaced.\nA target raster resolution must be specified.')

    warning('Prediction grid cells are not regularly spaced.\nData will be rasterised and covariate values averaged.')

    check.grid$z <- NULL
    coordinates(check.grid) <- ~x+y
    sp::proj4string(check.grid) <- coordinate.system

    # Create empty raster with desired resolution

    ras <- raster::raster(extent(check.grid), res = resolution)
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

    warning('New prediction grid (pred.grid) saved to global environment.')
    assign(x = 'pred.grid', prediction.grid, envir = .GlobalEnv)


  } # End if class(grid.regular)

  message("Computing ...")

  #'---------------------------------------------
  # Define reference and target systems
  #'---------------------------------------------

  reference <- segments[, covariate.names]
  target <- prediction.grid[, covariate.names]

  #'---------------------------------------------
  # Run the exdet tool from Mesgaran et al. (2014)
  #'---------------------------------------------

  mesgaran <- ExDet(ref = reference,
                    tg = target,
                    xp = covariate.names)

  #'---------------------------------------------
  # Add coordinates
  #'---------------------------------------------

  mesgaran <- prediction.grid %>%
    dplyr::select(x,y) %>%
    cbind(., mesgaran)

  #'---------------------------------------------
  # Return a list with univariate, combinatorial, and analog conditions as separate elements
  #'---------------------------------------------

  reslist <- list(data=NULL, rasters=NULL)

  reslist$data$all <- mesgaran

  reslist$data$univariate <- mesgaran %>%
    dplyr::filter(., ExDet < 0)

  reslist$data$combinatorial <- mesgaran %>%
    dplyr::filter(., ExDet > 1)

  reslist$data$analogue <- mesgaran %>%
    dplyr::filter(., ExDet >= 0 & ExDet <= 1)

  #'---------------------------------------------
  # Create rasters from extrapolation/MIC values
  #'---------------------------------------------

  reslist$rasters$ExDet <- reslist$data %>%
    purrr:::map(., ~ dplyr::select(., x, y, ExDet) %>%
                  safe_raster(.))%>%
    purrr::map(., "result")

  reslist$rasters$mic <- reslist$data %>%
    purrr:::map(., ~ dplyr::select(., x, y, mic) %>%
                  safe_raster(.)) %>%
    purrr::map(., "result")

  #'---------------------------------------------
  # Check that rasters have been produced for each extrapolation type
  #'---------------------------------------------

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
  #'---------------------------------------------
  # Project rasters
  #'---------------------------------------------

  for(r in 1:length(reslist$rasters$ExDet)){
    if(!is.null(reslist$rasters$ExDet[[r]]))raster::projection(reslist$rasters$ExDet[[r]]) <- coordinate.system}

  for(r in 1:length(reslist$rasters$mic)){
    if(!is.null(reslist$rasters$mic[[r]]))raster::projection(reslist$rasters$mic[[r]]) <- coordinate.system}

  message("Done!")

  #'---------------------------------------------
  # Print/save summary
  #'---------------------------------------------

  if(print.summary){

    if(save.summary){

      sumres <- summarise_extrapolation(extrapolation.object = reslist,
                                        covariate.names = covariate.names,
                                        extrapolation = TRUE,
                                        mic = TRUE,
                                        print.precision = print.precision)

      reslist <- append(x = reslist, values = list(summary = sumres))

    }else{

      summarise_extrapolation(extrapolation.object = reslist,
                              covariate.names = covariate.names,
                              extrapolation = TRUE,
                              mic = TRUE,
                              print.precision = print.precision)
    }

  }else{

    if(save.summary){

      sink("/dev/null")
      sumres <- summarise_extrapolation(extrapolation.object = reslist,
                                        covariate.names = covariate.names,
                                        extrapolation = TRUE,
                                        mic = TRUE,
                                        print.precision = print.precision)
      sink()
      reslist <- append(x = reslist, values = list(summary = sumres))

    }else{

    }
  }

  return(reslist)

}
