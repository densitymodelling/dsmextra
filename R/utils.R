#---------------------------------------------
# Function to tally the number/% of cells subject to extrapolation
#---------------------------------------------

n_and_p <- function(x){

  exl <- list(univariate.n = length(x[x < 0]),
              univariate.p = 100 * length(x[x < 0])/length(x),
              combinatorial.n = length(x[x > 1]),
              combinatorial.p = 100 * length(x[x > 1])/length(x),
              analogue.n = length(x[x >= 0 & x <=1]),
              analogue.p = 100 * length(x[x >= 0 & x <=1])/length(x))
  return(exl)
}

#---------------------------------------------
# Take mean of raster values - used in rasterize
#---------------------------------------------

mean_ras <- function(x, ...){mean(x, na.rm = TRUE)}

#---------------------------------------------
# Function to concatenate variables names where needed
#---------------------------------------------

collapse.vars <- function(l){

  if(length(l)==1){
    l <- unlist(l)
    if(length(l)>1) l <- paste(l, collapse = ", ")}

  if(length(l)>1){
    l <- purrr::map(.x = l, .f = ~data.frame(t(as.character(.x)))) %>%
      plyr::rbind.fill(.)

    for(i in 1:ncol(l)) l[,i] <- as.character(l[,i])
    l[is.na(l)] <- ""

    l <- apply(X = l,
               MARGIN = 1,
               function(x) {paste(x, collapse = " ") %>%
                   trimws(., which = "both") %>%
                   gsub(pattern = " ", replacement = ", ", x = .)})
  }
  return(l)
}

#---------------------------------------------
# Function to flip legend in leaflet
#---------------------------------------------
# Taken from: https://github.com/rstudio/leaflet/issues/256

addLegend_decreasing <- function (map, position = c("topright", "bottomright", "bottomleft", "topleft"), pal, r.values, na.label = "NA", bins = 7, colors, opacity = 0.5, labels = NULL, labFormat = labelFormat(), title = NULL, className = "info legend", layerId = NULL, group = NULL, data = getMapData(map), decreasing = FALSE) {

  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors))
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(r.values, "formula"))
      title <- deparse(r.values[[2]])
    r.values <- leaflet::evalFormula(r.values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && grDevices::col2rgb(na.color, alpha = TRUE)[[4]] ==
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins))
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1)
        pretty(r.values, bins)
      else bins

      if (length(bins) > 2)
        if (!all(abs(diff(bins, differences = 2)) <=
                 sqrt(.Machine$double.eps)))
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(r.values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")

    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }

    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- stats::quantile(r.values, probs = p, na.rm = TRUE)
      mids <- stats::quantile(r.values, probs = (p[-1] + p[-n])/2, na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(stats::na.omit(r.values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(r.values)))
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels))
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)),
                 na_color = na.color, na_label = na.label, opacity = opacity,
                 position = position, type = type, title = title, extra = extra,
                 layerId = layerId, className = className, group = group)
  invokeMethod(map, data, "addLegend", legend)
}

#---------------------------------------------
# Safe version of the rasterFromXYZ
#---------------------------------------------

# Returns a list with two elements: result and error (if one occurred)

safe_raster <- purrr::safely(function(x) suppressWarnings(raster::rasterFromXYZ(xyz = x)), otherwise = NULL)

safe_pts <- purrr::safely(function(x) suppressWarnings(SpatialPointsDataFrame(coords = x[,1:2], data = data.frame(x[,3]))))


#---------------------------------------------
# Functions used in Gower's distance calculations
#---------------------------------------------
# From Mannocci et al. 2018

# Function to standardise covariate values

rescale_cov <- function (ynew, y) { return ((ynew - mean(y, na.rm = TRUE)) / (stats::sd(y, na.rm = TRUE))) }

# Standardising prediction data simplifies computation A LOT!

make_X <- function (calibration_data,
                    test_data,
                    var_name){
  # Changed from: rescale_cov(ynew = test_data[, k], y = calibration_data[, k]) -- June 1st, 2021
  X <- sapply(var_name, function (k) {rescale_cov(ynew = test_data[[k]], y = calibration_data[[k]])})
  X <- as.data.frame (X)
  names (X) <- var_name
  return (X)}

#---------------------------------------------
# Function to project rasters to lat/lon for plotting
#---------------------------------------------

proj_rasters <- function(ll, coordinate.system){

  suppressWarnings(crs.webmerc <- sp::CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

  llr <- ll # Copy list

  # Univariate extrapolation is negative by definition
  # When only a small number of cells are subject to UE, the resampling
  # may result in the loss of some of them.
  # By recording the indices of UE cells, we can perform a simplistic
  # correction to make sure they show up on the map.

  analogue.xy <- raster::as.data.frame(llr$analogue, xy = TRUE) %>% stats::na.omit(.)
  analogue.xy <- sp::SpatialPointsDataFrame(coords = analogue.xy[, c("x", "y")],
                                            data = analogue.xy,
                                            proj4string = coordinate.system)
  analogue.xy <- sp::spTransform(analogue.xy, CRSobj = crs.webmerc)

  univariate.ind <- raster::Which(llr$univariate < 0, cells = TRUE)
  univariate.values <- llr$univariate[univariate.ind]

  univariate.xy <- raster::as.data.frame(llr$univariate, xy = TRUE) %>% stats::na.omit(.)
  univariate.xy <- sp::SpatialPointsDataFrame(coords = univariate.xy[, c("x", "y")],
                                              data = univariate.xy,
                                              proj4string = coordinate.system)
  univariate.xy <- sp::spTransform(univariate.xy, CRSobj = crs.webmerc)

  llr$all <- NULL
  llr <- purrr::discard(.x = llr, is.null)

  suppressWarnings(
  llr <- purrr::map(.x = llr, # Same extent as the full raster, allows correct alignment
                    .f = ~raster::projectRaster(from = .x,
                                                to = ll$all,
                                                method = 'ngb')) %>%
    purrr::map(.x = ., # CRS used by leaflet
               .f = ~raster::projectRaster(from = .,
                                           crs = crs.webmerc,
                                           method = 'ngb')))

  llr.univariate.ind <- raster::cellFromXY(object = llr$univariate,
                                           xy = sp::coordinates(univariate.xy))

  llr$univariate[llr.univariate.ind[which(is.na(llr$univariate[llr.univariate.ind]))]] <-    univariate.values[which(is.na(llr$univariate[llr.univariate.ind]))]

  r1 <- raster::as.data.frame(llr$univariate, xy = TRUE)
  r2 <- raster::as.data.frame(llr$analogue, xy = TRUE)
  names(r1) <- names(r2) <- c("x", "y", "ExDet")

  duplicate.cells <- rbind(r1, r2) %>%
    stats::na.omit(.) %>%
    dplyr::select(., x, y) %>%
    .[duplicated(.),]

  llr.analogue.ind <- raster::cellFromXY(object = llr$analogue,
                                         xy = duplicate.cells)


  llr$analogue[llr.analogue.ind] <- NA


  return(llr)}

#---------------------------------------------
# Function to check correct CRS input
#---------------------------------------------

check_crs <- function(coordinate.system){

if(!class(coordinate.system)=="CRS"){

  coord.err <- tryCatch(expr = sp::CRS(coordinate.system),
                        error = function(e) return(NA))

  if(is.na(coord.err)){stop('Unrecognised coordinate system')
  }else{supressWarnings(coordinate.system <- sp::CRS(coordinate.system))}
}

  return(coordinate.system)

}

