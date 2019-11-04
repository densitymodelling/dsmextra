## ----eval=FALSE, message=FALSE, warning=FALSE, include=TRUE--------------
#  install.packages("dsmextra")

## ----eval=TRUE, message=FALSE, warning=FALSE, include=TRUE---------------
library(dsmextra)
library(raster)
library(magrittr)

## ----eval=FALSE, include=TRUE--------------------------------------------
#  devtools::install_github("densitymodelling/dsmextra", dependencies = TRUE, build_vignettes = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE-------------------------

#'--------------------------------------------------------------------
# Set tibble options
#'--------------------------------------------------------------------

options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

#'--------------------------------------------------------------------
# Set knitr options
#'--------------------------------------------------------------------

knitr::opts_chunk$set(echo = TRUE)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
#'---------------------------------------------
# Load and extract the data
#'---------------------------------------------
data("spermwhales")

segs <- spermwhales$segs
predgrid <- spermwhales$predgrid

#'---------------------------------------------
# Quick inspection
#'---------------------------------------------
knitr::kable(head(segs[,c(1, 3:12)]), format = "pandoc")
knitr::kable(head(predgrid), format = "pandoc")

## ----echo=TRUE, message=FALSE, warning=FALSE, include=FALSE--------------
load('../R/sysdata.rda')

## ----echo=TRUE, message=FALSE, warning=FALSE, fig.cap="<b>Fig. 1:</b> Distribution of sperm whale (<em>Physeter macrocephalus</em>) sightings within the study area off the US East coast. Surveyed transects are shown as solid lines."----
#'---------------------------------------------
# Create an outline of the study area boundaries
#'---------------------------------------------

study_area <- predgrid[, c("x", "y")]
study_area$value <- 1
study_area <- raster::rasterFromXYZ(study_area, crs = sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
study_area <- raster::rasterToPolygons(study_area, dissolve = TRUE)
study_area <- sp::spTransform(study_area, CRSobj = sp::proj4string(transects))
study_area <- smoothr::smooth(study_area, method = "ksmooth", smoothness = 5)

#'---------------------------------------------
# Produce a simple plot
#'---------------------------------------------

plot(study_area, col = "lightskyblue1") # Region boundary
plot(transects, add = TRUE, col = "skyblue3") # Survey tracks
maps::map("world", fill = TRUE, col = "grey",
          xlim = range(obs$coords.x1),
          ylim = range(obs$coords.x2), add = TRUE)
pts <- obs # Sightings
coordinates(pts) <- ~ coords.x1 + coords.x2
axis(1); axis(2); box()
points(pts, pch = 16)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
#'---------------------------------------------
# Define projected coordinate system
#'---------------------------------------------

aftt_crs <- sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

#'---------------------------------------------
# Define environmental covariates of interest
#'---------------------------------------------

covariates.spermwhale <- c("Depth", "SST", "NPP", "DistToCAS", "EKE")

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
spermwhale.extrapolation <- compute_extrapolation(segments = segs, 
                                   covariate.names = covariates.spermwhale, 
                                   prediction.grid = predgrid,
                                   coordinate.system = aftt_crs,
                                   print.summary = TRUE,
                                   save.summary = FALSE,
                                   print.precision = 2)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
# Number of cells subject to univariate extrapolation (see below for definition)
predgrid %>% 
  dplyr::filter(!dplyr::between(Depth, min(segs$Depth), max(segs$Depth)) | 
                  !dplyr::between(SST, min(segs$SST), max(segs$SST)) |
                  !dplyr::between(NPP, min(segs$NPP), max(segs$NPP)) |
                  !dplyr::between(DistToCAS, min(segs$DistToCAS), max(segs$DistToCAS)) |
                  !dplyr::between(EKE, min(segs$EKE), max(segs$EKE))) %>% 
  nrow()

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
tibble::glimpse(spermwhale.extrapolation)

## ----ExDet, echo=FALSE, out.width = '50%', fig.align="center", fig.cap="<b>Fig. 2:</b> Schematic representation of extrapolation in multivariate environmental space, based on two hypothetical covariates. Reference data points are represented as grey circles. Shaded areas correspond to different types of extrapolation outside the envelope of the reference data. Univariate extrapolation occurs beyond the range of individual covariates. Combinatorial extrapolation occurs within this range, but outside the reference hyperspace/hypervolume. <u>Source</u>: @Bouchet2019, adapted from @Mesgaran2014."----
knitr::include_graphics("dsmextra_Figure2.png")

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
# Example from combinatorial extrapolation
head(spermwhale.extrapolation$data$combinatorial)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
summary.values <- summarise_extrapolation(extrapolation.object = spermwhale.extrapolation,
                        covariate.names = covariates.spermwhale,
                        extrapolation = TRUE,
                        mic = TRUE,
                        print.precision = 2)
summary.values

## ----echo=TRUE, message=TRUE, warning=FALSE------------------------------
compare_covariates(extrapolation.type = "both",
                   covariate.names = covariates.spermwhale,
                   n.covariates = NULL,
                   segments = segs, 
                   prediction.grid = predgrid,
                   coordinate.system = aftt_crs,
                   create.plots = TRUE,
                   display.percent = TRUE)

## ----neighbour, echo=FALSE, out.width = '85%', fig.align="center", fig.cap="<b>Fig. 3:</b> Conceptual representation of two key extrapolation metrics. (A) Distance from the sampled environmental space. A target point far from the envelope of the reference data (eg. falling in the yellow area) is arguably 'more of an extrapolation' than one close to it (eg. falling in the purple area).  (B) Neighbourhood (syn. percentage of data nearby). Due to the shape of the reference data cloud, the amount of sample information available to 'inform' predictions made at target points can vary considerably. <u>Source</u>: @Bouchet2019."----
knitr::include_graphics("dsmextra_Figure3.png")

## ----echo=TRUE, message=TRUE, warning=FALSE------------------------------
#'---------------------------------------------
# Calculate Gower's distances and %N
#'---------------------------------------------

spermwhale.nearby <- compute_nearby(segments = segs, 
                                    prediction.grid = predgrid, 
                                    coordinate.system = aftt_crs,
                                    covariate.names = covariates.spermwhale,
                                    nearby = 1)

## ----echo=TRUE, message=FALSE, warning=FALSE-----------------------------
#'---------------------------------------------
# Rename coordinates and convert to SpatialPointsdf
#'---------------------------------------------

obs.sp <- obs %>% 
  dplyr::rename(., x = coords.x1, y = coords.x2) %>% 
  sp::SpatialPointsDataFrame(coords = cbind(.$x, .$y), data = ., proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
  sp::spTransform(., CRSobj = aftt_crs)


## ----echo=TRUE, message=FALSE, warning=FALSE, fig.cap=" "----------------
map_extrapolation(map.type = "extrapolation",
                  extrapolation.values = spermwhale.extrapolation,
                  covariate.names = covariates.spermwhale,
                  prediction.grid = predgrid, 
                  coordinate.system = aftt_crs,
                  sightings = obs.sp, 
                  tracks = transects)

## ----echo=TRUE, message=FALSE, warning=FALSE, fig.cap=" "----------------

map_extrapolation(map.type = "mic",
                  extrapolation.values = spermwhale.extrapolation,
                  covariate.names = covariates.spermwhale,
                  prediction.grid = predgrid, 
                  coordinate.system = aftt_crs,
                  sightings = obs.sp, 
                  tracks = transects)

## ----echo=TRUE, message=FALSE, warning=FALSE, fig.cap=" "----------------
map_extrapolation(map.type = "nearby",
                  gower.values = spermwhale.nearby,
                  covariate.names = covariates.spermwhale,
                  prediction.grid = predgrid, 
                  coordinate.system = aftt_crs,
                  sightings = obs.sp, 
                  tracks = transects)

## ----echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE-----------------
#  #'---------------------------------------------
#  # One function to rule them all, one funtion to find them,
#  # One function to bring them all, and in the darkness bind them.
#  #'---------------------------------------------
#  
#  spermwhale.analysis <- extrapolation_analysis(segments = segs,
#                                            covariate.names = covariates.spermwhale,
#                                            prediction.grid = predgrid,
#                                            coordinate.system = aftt_crs,
#                                            summarise.extrapolation = TRUE,
#                                            summary.print.precision = 2,
#                                            compare.covariates = TRUE,
#                                            compare.extrapolation.type = "both",
#                                            compare.n.covariates = NULL,
#                                            compare.create.plots = TRUE,
#                                            compare.display.percent = TRUE,
#                                            nearby.compute = TRUE,
#                                            nearby.nearby = 1,
#                                            map.generate = TRUE,
#                                            map.sightings = obs.sp,
#                                            map.tracks = transects)

