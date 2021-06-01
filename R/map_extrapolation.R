#' Interactive maps of extrapolation
#'
#' Produces interactive html maps of extrapolation values in the prediction area. The function relies on the \code{\link[leaflet]{leaflet}} package (Cheng et al. 2018), and thus requires an internet connection (i.e. will not work offline).
#'
#' @importFrom raster as.data.frame as.factor
#' @import leaflet
#'
#' @param map.type Character string. Type of map to be returned. Either \code{extrapolation} for an extrapolation map, \code{mic} for a map of the most influential covariates, or \code{nearby} for a map of the percentage of data nearby.
#' @param extrapolation.object List object as returned by \link{compute_extrapolation} or \link{compute_nearby}.
#' @param base.layer Base layer used for mapping. The default is \code{ocean}, which uses the ESRI.OceanBasemap. Use \code{world} for ESRI.WorldImagery and \code{gray} for ESRI.WorldGrayCanvas. Available map tiles can be viewed at \href{https://leaflet-extras.github.io/leaflet-providers/preview/}{https://leaflet-extras.github.io/leaflet-providers/preview/}.
#' @param sightings Species observations (optional). Can be supplied as a \code{matrix} of coordinates, a \code{data.frame}, a \code{\link[sp]{SpatialPoints}} object or a \code{\link[sp]{SpatialPointsDataFrame}} object. Circle markers will be proportional to group size if the data contain a column labelled \code{size}.
#' @param tracks Survey tracks (optional). Can be supplied as a \code{matrix} of coordinates, a \code{data.frame}, a \code{\link[sp]{SpatialLines}} object or a \code{\link[sp]{SpatialLinesDataFrame}} object. A \code{TransectID} field is required for matrix or data.frame inputs.
#' @param verbose Logical. Show or hide possible warnings and messages.
#'
#' @return An interactive html map.
#'
#' @author Phil J. Bouchet
#'
#' @references Bouchet PJ, Miller DL, Roberts JJ, Mannocci L, Harris CM and Thomas L (2019). From here and now to there and then: Practical recommendations for extrapolating cetacean density surface models to novel conditions. CREEM Technical Report 2019-01, 59 p. \href{https://research-repository.st-andrews.ac.uk/handle/10023/18509}{https://research-repository.st-andrews.ac.uk/handle/10023/18509}
#'
#' Cheng J, Karambelkar B, Xie Y (2018). leaflet: Create Interactive Web Maps with the JavaScript 'Leaflet' Library. R package version 2.0.2. \href{https://CRAN.R-project.org/package=leaflet}{https://CRAN.R-project.org/package=leaflet}
#'
#' @export
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
#'  # Define covariates of interest
#'  my_cov <- c("Depth", "DistToCAS", "SST", "EKE", "NPP")
#'
#' # Assess extrapolation in the multivariate space defined by five covariates
#' spermw.extrapolation <- compute_extrapolation(samples = segs,
#'       covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
#'       prediction.grid = predgrid,
#'       coordinate.system = my_crs)
#'
#' # Assess the percentage of data nearby
#' spermw.nearby <- compute_nearby(samples = segs,
#'                                prediction.grid = predgrid,
#'                                coordinate.system = my_crs,
#'                                covariate.names = my_cov,
#'                                nearby = 1)
#'
#' # Generate maps
#' map_extrapolation(map.type = "extrapolation", extrapolation.object = spermw.extrapolation)
#' map_extrapolation(map.type = "mic", extrapolation.object = spermw.extrapolation)
#' map_extrapolation(map.type = "nearby", extrapolation.object = spermw.nearby)

map_extrapolation <- function(map.type = NULL,
                              extrapolation.object = NULL,
                              base.layer = "ocean",
                              sightings = NULL,
                              tracks = NULL,
                              verbose = TRUE){

  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------

  if(is.null(map.type)) stop("Argument 'maptype' must be specified.")
  if(!map.type%in%c("extrapolation", "nearby", "mic")) stop("Unknown map type requested.")
  if(is.null(extrapolation.object) & map.type == "extrapolation") stop("Argument 'extrapolation.values' cannot be NULL when maptype is set to 'extrapolation'.")
  if(is.null(extrapolation.object) & map.type == "mic") stop("Argument 'extrapolation.values' cannot be NULL when maptype is set to 'mic'.")
  if(!map.type %in% extrapolation.object$type) stop("Map type undefined for the input extrapolation.object")

  #---------------------------------------------
  # Extract data
  #---------------------------------------------

  covariate.names <- extrapolation.object$covariate.names
  prediction.grid <- extrapolation.object$prediction.grid
  coordinate.system <- extrapolation.object$coordinate.system
  if(map.type == "nearby") gower.values <- extrapolation.object$raster
  if(map.type %in% c("extrapolation", "mic")) extrapolation.values <- extrapolation.object


  if(!is.null(sightings)){

    if(!class(sightings)%in%c("data.frame", "matrix", "SpatialPoints", "SpatialPointsDataFrame")) stop("Sightings must be of class matrix, data.frame, SpatialPoints or SpatialPointsDataFrame")

    if(is.data.frame(sightings) | is.matrix(sightings)){

      if(sum(c("x","y")%in%names(sightings))<2) {
        if(verbose) message("Missing x,y coordinates: sightings not shown")
        sightings = NULL
        latlon.sightings = NULL}

      if(!is.null(sightings)) latlon.sightings <- sum(all(range(sightings$x)>=-180 & range(sightings$x)<=180))

    }else{coords.sightings <- sp::coordinates(sightings)
    latlon.sightings <- sum(all(as.vector(apply(coords.sightings,2, range))>=-180 & as.vector(apply(coords.sightings,2, range))<=180))
    }
  }

  if(!is.null(tracks)){

    if(!class(tracks)%in%c("data.frame", "matrix", "SpatialLine", "SpatialLines", "SpatialLinesDataFrame")) stop("Tracks must be of class matrix, data.frame, SpatialLines or SpatialLinesDataFrame")

    if(is.data.frame(tracks) | is.matrix(tracks)){

      if(sum(c("x","y")%in%names(tracks))<2) {
        if(verbose) message("Missing x,y coordinates: tracks not shown")
        tracks = NULL
        latlon.tracks = NULL}

      if(!is.null(tracks)) latlon.tracks <- sum(all(range(tracks$x)>=-180 & range(tracks$x)<=180))

    }else{

      coords.tracks <- sp::coordinates(tracks)

      if(is.list(coords.tracks)){
        coords.tracks <- purrr::flatten(coords.tracks)
        coords.tracks <- do.call("rbind", coords.tracks)}

      latlon.tracks <- sum(all(as.vector(apply(coords.tracks, 2, range))>=-180 & as.vector(apply(coords.tracks, 2, range))<=180))
    }
  }

  baselyr <- switch(base.layer,
                    "ocean" = "Esri.OceanBasemap",
                    "world" = "Esri.WorldImagery",
                    "gray" = "Esri.WorldGrayCanvas")

  coordinate.system <- check_crs(coordinate.system = coordinate.system)

  #---------------------------------------------
  # Define coordinate systems
  #---------------------------------------------

  suppressWarnings(crs.ll <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  # EPSG:3857 (Web Mercator) is needed by leaflet for plotting

  suppressWarnings(crs.webmerc <- sp::CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

  #---------------------------------------------
  # Check which type of extrapolation occurred
  #---------------------------------------------

  if(map.type %in% c("mic", "extrapolation")){
  types <- purrr::map_dbl(.x = extrapolation.values$data[2:4],
                          .f = ~nrow(.))
  types <- c("Univariate", "Combinatorial", "Analogue")[types>0]}

  #---------------------------------------------
  # Defines survey extent (in lat/lon coords)
  #---------------------------------------------

  survey_extent <- methods::as(raster::extent(range(prediction.grid$x),
                              range(prediction.grid$y)),
                      "SpatialPolygons")

  sp::proj4string(survey_extent) <- coordinate.system
  survey_extent <- sp::spTransform(survey_extent, CRSobj = crs.ll)
  survey_extent <- raster::extent(survey_extent)

  #---------------------------------------------
  # Converts segments to SpatialLines
  #---------------------------------------------

  if(!is.null(tracks)){

    if(latlon.tracks==1){crstracks <- crs.ll
    }else{crstracks <- coordinate.system}

    if(is.data.frame(tracks) | is.matrix(tracks)){

      if("TransectID"%in%names(tracks)){

        tracks <- tracks %>%
          split(.$ID) %>%
          purrr::map(.x = ., .f = ~dplyr::select(.x, x, y) %>%
                       as.matrix(.) %>%
                       raster::spLines(., crs = crstracks)) %>%
          do.call("rbind", .) %>%
          sp::SpatialLinesDataFrame(., data = tracks) %>%
          sp::spTransform(., CRSobj = crs.ll)

      }else{

        if(verbose) message("No transect ID detected: Linking all segments")

        tracks <- tracks %>%
          dplyr::select(x, y) %>%
          as.matrix(.) %>%
          raster::spLines(., crs = crstracks) %>%
          sp::SpatialLinesDataFrame(., data = tracks) %>%
          sp::spTransform(., CRSobj = crs.ll)

      }


    }else{tracks <- sp::spTransform(tracks, CRSobj = crs.ll)}
  }

  #---------------------------------------------
  # Sightings
  #---------------------------------------------

  if(!is.null(sightings)){

    if(latlon.sightings==1){crssightings <- crs.ll
    }else{crssightings <- coordinate.system}

    if(is.data.frame(sightings) | is.matrix(sightings)){

      sp::coordinates(sightings) <- ~ x + y
      sp::proj4string(sightings) <- crssightings

    }else{sightings <- sp::spTransform(x = sightings, CRSobj = crs.ll)}
  }

  #---------------------------------------------
  # Basemap
  #---------------------------------------------

  exleaf <- leaflet::leaflet(sightings) %>%
    leaflet::setView(lng = mean(c(survey_extent[1], survey_extent[2])),
                     lat = mean(c(survey_extent[3],survey_extent[4])),
                     zoom = 5) %>%
    leaflet::fitBounds(survey_extent[1],
                       survey_extent[3],
                       survey_extent[2],
                       survey_extent[4]) %>%

    leaflet::addProviderTiles(provider = baselyr)
    # leaflet::addProviderTiles(provider = providers[which(names(providers)==baselyr)][[1]])


  if(map.type == "extrapolation"){

    #---------------------------------------------
    # Project rasters to lat/lon for plotting
    #---------------------------------------------

    projr <- proj_rasters(ll = extrapolation.values$rasters$ExDet, coordinate.system = coordinate.system)

    #---------------------------------------------
    # Adds rasters and colour palettes
    #---------------------------------------------

    if("univariate"%in%names(projr)){

      pal.univariate <- leaflet::colorNumeric(c("#7f2704", "#fd8d3c", "#fee6ce"),
                                              raster::getValues(projr$univariate),na.color = "transparent")
suppressWarnings(
      exleaf <- exleaf %>%
        leaflet::addRasterImage(map = .,
                                x = projr$univariate,
                                colors = pal.univariate,
                                group = "Univariate",
                                project = FALSE,
                                opacity = 1) %>%

        #---------------------------------------------
        # legend
        #---------------------------------------------

        addLegend_decreasing(map = .,
                             pal = pal.univariate,
                             opacity = 1,
                             r.values = raster::getValues(projr$univariate),
                             decreasing = TRUE,
                             group = "Univariate",
                             title = "Univariate")
      )}

    if("combinatorial"%in%names(projr)){


      if(base.layer=="ocean"){
        pal.combinatorial <- leaflet::colorNumeric(c("#eadef7", "#a56bd6", "#47026d"),
                                                   raster::getValues(projr$combinatorial), na.color = "transparent")
      }else{
        pal.combinatorial <- leaflet::colorNumeric(c("#deebf7", "#6baed6", "#08519c"),
                                                   raster::getValues(projr$combinatorial), na.color = "transparent")}

      suppressWarnings(
      exleaf <- exleaf %>%
        leaflet::addRasterImage(map = .,
                                x = projr$combinatorial,
                                colors = pal.combinatorial,
                                group = "Combinatorial",
                                project = FALSE,
                                opacity = 1) %>%

        #---------------------------------------------
        # legend
        #---------------------------------------------

        addLegend_decreasing(map = .,
                             pal = pal.combinatorial,
                             opacity = 1,
                             decreasing = TRUE,
                             group="Combinatorial",
                             r.values = raster::getValues(projr$combinatorial),
                             title = "Combinatorial")
      )}

    if("analogue"%in%names(projr)){

      pal.analogue <- leaflet::colorNumeric(c("#e5f5e0", "#74c476", "#00441b"),
                                            raster::getValues(projr$analogue),
                                            na.color = "transparent")

      suppressWarnings(
      exleaf <- exleaf %>%
        leaflet::addRasterImage(map = .,
                                x = projr$analogue,
                                colors = pal.analogue,
                                group ="Analogue",
                                project = FALSE,
                                opacity = 1) %>%

        #---------------------------------------------
        # legend
        #---------------------------------------------

        addLegend_decreasing(map = .,
                             pal = pal.analogue,
                             opacity = 1,
                             decreasing = TRUE,
                             group="Analogue",
                             r.values = raster::getValues(projr$analogue),
                             title = "Analogue")
      )}

  }else if(map.type == "mic"){

    #---------------------------------------------
    # Adds rasters and colour palettes
    #---------------------------------------------

    pal <- leaflet::colorFactor(palette = c("#B2FF8C", dichromat::colorschemes$Categorical.12[c(1:4, 7:12)]),
                       domain = extrapolation.values$data$all$mic,
                       na.color = "transparent")

    micvars <- extrapolation.values$data$all %>%
      dplyr::count(mic)

    allvars <- tibble::tibble(vars = c("None", covariate.names), ID = seq(0, length(covariate.names)))

    mapvars <- dplyr::left_join(x = micvars, y = allvars, by = c("mic" = "ID"))

    suppressWarnings(
    exleaf <- exleaf %>% leaflet::addRasterImage(map = .,
                                                 x = raster::as.factor(extrapolation.values$rasters$mic$all),
                                                 colors = pal,
                                                 group = "MIC",
                                                 opacity = 1) %>%

      #---------------------------------------------
      # legend
      #---------------------------------------------

      addLegend_decreasing(map = .,
                           pal = pal,
                           labFormat = labelFormat(
                             prefix = "",
                             suffix = paste0(" - ", mapvars$vars)),
                           opacity = 1,
                           decreasing = TRUE,
                           group = "MIC",
                           r.values = raster::getValues(extrapolation.values$rasters$mic$all),
                           title = "MIC"))


  }else if(map.type == "nearby"){

    #---------------------------------------------
    # Adds rasters and colour palettes
    #---------------------------------------------

    pal.near <- leaflet::colorNumeric(pals::parula(100),
                                      raster::getValues(gower.values),
                                      na.color = "transparent")
    suppressWarnings(
    exleaf <- exleaf %>%
      leaflet::addRasterImage(map = .,
                              colors = pal.near,
                              x = gower.values,
                              group="% nearby",
                              opacity = 1) %>%

      #---------------------------------------------
      # legend
      #---------------------------------------------

      addLegend_decreasing(map = .,
                           pal = pal.near,
                           opacity = 1,
                           decreasing = TRUE,
                           group ="% nearby",
                           r.values = raster::getValues(gower.values),
                           title = "% nearby"))


  }


  #---------------------------------------------
  # Adds survey tracks and sightings
  #---------------------------------------------

  if(is.null(tracks)==FALSE)  exleaf <- exleaf %>%
    leaflet::addPolylines(data = tracks,
                          col = "black",
                          weight = 1,
                          group = "Tracks")

  if(is.null(sightings)==FALSE){

    if("size"%in%names(sightings)){

      # Create html labels for group sizes

      sizelist <- mapply(function(x, y) {
        htmltools::HTML(sprintf("%s: %s",
                                htmltools::htmlEscape(x),
                                htmltools::htmlEscape(y)))},
        "Group size", sightings$size, SIMPLIFY = FALSE)
      sizelist <- purrr::set_names(sizelist, NULL)

      exleaf <- exleaf %>%
        leaflet::addCircleMarkers(lng = sp::coordinates(sightings)[,1],
                                  lat = sp::coordinates(sightings)[,2],
                                  radius=~(log(sightings$size)+1)*3,
                                  label = sizelist,
                                  labelOptions = lapply(1:nrow(sightings),
                                                        function(x) {
                                                          labelOptions(direction='auto')}),
                                  stroke = FALSE,
                                  color="black",
                                  fillOpacity = 0.5,
                                  group = "Sightings")
    }else{

      if(verbose) message("No 'size' column detected: Group size not shown")

      exleaf <- exleaf %>%
        leaflet::addCircleMarkers(lng = sp::coordinates(sightings)[,1],
                                  lat = sp::coordinates(sightings)[,2],
                                  stroke = FALSE,
                                  radius = 4,
                                  color="black",
                                  fillOpacity = 0.5,
                                  group = "Sightings")

    }
  }

  #---------------------------------------------
  # Layer controls
  #---------------------------------------------

  toggles <- paste0(switch(map.type,
                           "extrapolation" = "1",
                           "mic" = "2",
                           "nearby" = "3"),
                    paste0(as.character(as.numeric(c(!is.null(tracks),
                                                     !is.null(sightings)))), collapse=""))

  lyr.controls <- switch(toggles,
                         "001" = c("Sightings"),
                         "010" = c("Tracks"),
                         "011" = c("Tracks", "Sightings"),
                         "100" = c(types),
                         "101" = c(types, "Sightings"),
                         "110" = c(types, "Tracks"),
                         "111" = c(types, "Tracks", "Sightings"),
                         "200" = c("MIC"),
                         "201" = c("MIC", "Sightings"),
                         "210" = c("MIC", "Tracks"),
                         "211" = c("MIC", "Tracks", "Sightings"),
                         "300" = c("% nearby"),
                         "301" = c("% nearby", "Sightings"),
                         "310" = c("% nearby", "Tracks"),
                         "311" = c("% nearby", "Tracks", "Sightings"))

  exleaf <- exleaf %>%
    addLayersControl(position = "topleft",
                     overlayGroups = lyr.controls,
                     options = layersControlOptions(collapsed = FALSE))

  if(verbose) warning('map_extrapolation relies on the leaflet package, which is built around a Web Mercator projection (EPSG:3857), and therefore requires rasters to be reprojected for plotting. As a result, minor discrepancies may  occur between the interactive maps shown in the viewer, and the underlying raw data. The latter can be accessed directly from extrapolation object returned by <compute_extrapolation> and visualised using alternative packages such as ggplot2.')

  return(exleaf)
}
