#' Quantify the percentage of data nearby
#'
#' Calculates the proportion of sample (reference) points lying in the vicinity of prediction (target) points in multivariate environmental space. See the 'Details' section and King & Zeng (2007) for a technical explanation.
#'
#' While extrapolation is often seen as a binary concept (i.e. it either does or does not take place), it is reasonable to expect that predictions made at target points situated just outside the sampled environmental space may be more reliable than those made at points far outside it. The ExDet tool available through \code{\link{compute_extrapolation}} inherently quantifies this notion of \emph{distance} from the envelope of the reference data.
#'
#' However, the multivariate distribution of reference data points is often far from homogeneous. It is possible, therefore, for target points representing analogue conditions to fall within sparsely sampled regions of the reference space; or conversely, for two target points reflecting an equal degree of extrapolation to have very different amounts of reference data within their vicinity.
#'
#' The notion of neighbourhood (or percentage/proportion of data nearby, \%N) captures this idea, and provides an additional measure of the reliability of extrapolations in multivariate environmental space (Virgili et al. 2017; Mannocci et al. 2018). In practice, \%N for any target point can be defined as the proportion of reference data within a radius of one geometric mean Gower’s distance (G^{2}), calculated between all pairs of reference points (King and Zeng 2007). The Gower’s distance between two points \emph{i} and \emph{j} defined along the axes of \emph{K} covariates is calculated as the average absolute distance between the values of these two points in each dimension, divided by the range of the data, such that:
#'
#' \eqn{G_{ij}^2=\frac{1}{K}\sum_{k=1}^{K}\frac{\left|x_{ik}-x_{jk}\right|}{\textrm{max}(X_k)-\textrm{min}(X_k)}}
#'
#' The \code{\link{compute_nearby}} function is adapted from the code given in Mannocci et al. (2018) and allows the calculation of Gower’s distances as a basis for defining the neighbourhood.
#'
#' In addition, the \code{\link[WhatIf]{whatif}} function from the \href{https://CRAN.R-project.org/package=WhatIf}{Whatif} package (Gandrud et al. 2017), which \code{compute_nearby} calls internally, may not run on very large datasets. Running calculations on partitions of the data may circumvent this problem and lead to speed gains. Two arguments can be used to do this:
#' \tabular{ll}{
#'   \code{max.size} \tab Threshold above which partitioning will be triggered \cr
#'   \code{no.partitions} \tab Number of required partitions \cr
#'  }
#' In practice, a run of \code{compute_nearby} begins with a quick assessment of the dimensions of the input data, i.e. the reference and target data.frames. If the product of their dimensions (i.e. number of samples multiplied by number of prediction grid cells) exceeds the value set for \code{max.size}, then \code{no.partitions} subsets of the data will be created and the computations run on each using \code{\link[purrr]{map}} functions from the \href{https://cran.r-project.org/web/packages/purrr/index.html}{purrr} package (Henry and Wickham 2019). This means that a smaller \code{max.size} will trigger partitioning on correspondingly smaller datasets. By default, \code{max.size} is set to \code{1e7}. This value was chosen arbitrarily, and should be sufficiently large as to obviate the need for partitioning on most datasets.
#'
#' @inheritParams compute_extrapolation
#' @param nearby Scalar indicating which reference data points are considered to be 'nearby' (i.e. within ‘nearby’ mean geometric Gower's distances of) prediction points. Defaults to 1, as per Mannocci et al. (2018) and Virgili et al. (2017).
#' @param max.size Minimum size threshold for partitioning computations. Calculated as \code{\link[base]{prod}(\link[base]{nrow}(samples),\link[base]{nrow}(prediction.grid))}. Has a default value of \code{1e7}. See the 'Details' section.
#' @param no.partitions Integer. Number of desired partitions of the data (default of 10). See the 'Details' section.
#' @param resolution Resolution of the output raster (in units relevant to \code{coordinate.system}). Only required if \code{prediction.grid} is irregular, and thus needs to be rasterised. Defaults to NULL.
#' @return A raster object mapping the proportion of reference data nearby each point in \code{prediction.grid}.
#' @author Phil J. Bouchet, Laura Mannocci.
#' @references Bouchet PJ, Miller DL, Roberts JJ, Mannocci L, Harris CM and Thomas L (2019). From here and now to there and then: Practical recommendations for extrapolating cetacean density surface models to novel conditions. CREEM Technical Report 2019-01, 59 p. \href{https://research-repository.st-andrews.ac.uk/handle/10023/18509}{https://research-repository.st-andrews.ac.uk/handle/10023/18509}
#'
#' Gandrud C, King G, Stoll H, Zeng L (2017). WhatIf: Evaluate Counterfactuals. R package version 1.5-9. \href{https://CRAN.R-project.org/package=WhatIf}{https://CRAN.R-project.org/package=WhatIf}.
#'
#' Henry L, Wickham H (2019). purrr: Functional Programming Tools. R package version 0.3.2. \href{https://CRAN.R-project.org/package=purrr}{https://CRAN.R-project.org/package=purrr}.
#'
#' King G, Zeng L (2007). When can history be our guide? The pitfalls of counterfactual inference. International Studies Quarterly 51, 183–210. DOI: \href{https://www.jstor.org/stable/pdf/4621707.pdf?seq=1#page_scan_tab_contents}{10.1111/j.1468-2478.2007.00445.x}
#'
#' Mannocci L, Roberts JJ, Halpin PN, Authier M, Boisseau O, Bradai MN, Canãdas A, Chicote C, David L, Di-Méglio N, Fortuna CM, Frantzis A, Gazo M, Genov T, Hammond PS, Holcer D, Kaschner K, Kerem D, Lauriano G, Lewis T, Notarbartolo Di Sciara G, Panigada S, Raga JA, Scheinin A, Ridoux V, Vella A, Vella J (2018). Assessing cetacean surveys throughout the mediterranean sea: A gap analysis in environmental space. Scientific Reports 8, art3126. DOI: \href{https://www.nature.com/articles/s41598-018-19842-9}{10.5061/dryad.4pd33}.
#'
#' Virgili A, Racine M, Authier M, Monestiez P, Ridoux V (2017). Comparison of habitat models for scarcely detected species. Ecological Modelling 346, 88–98. DOI: \href{https://www.sciencedirect.com/science/article/pii/S0304380016308146}{10.1016/j.ecolmodel.2016.12.013}.
#'
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
#' # Assess the percentage of data nearby
#' spermw.nearby <- compute_nearby(samples = segs,
#'                                prediction.grid = predgrid,
#'                                coordinate.system = my_crs,
#'                                covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
#'                                nearby = 1)
compute_nearby <- function (samples,
                            covariate.names,
                            prediction.grid,
                            coordinate.system,
                            nearby,
                            max.size = 1e7,
                            no.partitions = 10,
                            resolution = NULL) {

  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------

  if(nearby<=0) stop("nearby must be strictly positive")
  if(!is.numeric(nearby)) stop("Non-numeric input to argument: nearby")
  if(max.size<=0) stop("max.size must be strictly positive")
  if(!is.numeric(max.size)) stop("Non-numeric input to argument: max.size")
  if(no.partitions>nrow(prediction.grid)) stop("Number of partitions too large")

  coordinate.system <- check_crs(coordinate.system = coordinate.system)

  #---------------------------------------------
  # Check if prediction grid is regular
  #---------------------------------------------

  check.grid <- prediction.grid %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(z = 1)

  grid.regular <- try(raster::rasterFromXYZ(check.grid), silent = TRUE)

  #---------------------------------------------
  # If grid is irregular, rasterise prediction.grid based on specified resolution
  #---------------------------------------------

  if(class(grid.regular)=="try-error"){

    if(is.null(resolution)) stop('Prediction grid cells are not regularly spaced.\nA target raster resolution must be specified.')

    warning('Prediction grid cells are not regularly spaced.\nData will be rasterised and covariate values averaged.')

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

    # warning('New prediction grid (pred.grid) saved to global environment.')
    # assign(x = 'pred.grid', prediction.grid, envir = .GlobalEnv)


  } # End if class(grid.regular)

  #---------------------------------------------
  # Check size of input datasets
  #---------------------------------------------

  big.data <- ifelse(prod(nrow(samples), nrow(prediction.grid)) > max.size, TRUE, FALSE)

  #---------------------------------------------
  # Compute counterfactuals
  #---------------------------------------------

  if(big.data){

    counterfact <- whatif.opt(formula = NULL,
                              data = make_X(calibration_data = samples,
                                            test_data = samples,
                                            covariate.names),
                              cfact = make_X(calibration_data = samples,
                                             test_data = prediction.grid,
                                             covariate.names),
                              nearby = nearby,
                              no.partitions = no.partitions)

  }else{

    counterfact <- whatif(formula = NULL,
                           data = make_X(calibration_data = samples,
                                         test_data = samples,
                                         covariate.names),
                           cfact = make_X(calibration_data = samples,
                                          test_data = prediction.grid,
                                          covariate.names),
                           nearby = nearby,
                           choice = "distance")

  }

  #---------------------------------------------
  # Convert to raster
  #---------------------------------------------

  rgow <- cbind(prediction.grid[, c("x", "y")],
                100*counterfact$sum.stat) # in percent
  names(rgow)[3] <- 'perc_nearby'

  rgow <- raster::rasterFromXYZ(xyz = rgow,
                                crs = coordinate.system)

  message('Done!')
  return(rgow)
}
