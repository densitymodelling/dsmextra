#' Compare different combinations of covariates
#'
#' Summarises the extent of  univariate (Type I) and combinatorial (Type II) extrapolation associated with different combinations of input covariates.
#'
#' The extent and magnitude of extrapolation naturally vary with the type and number of covariates considered. It may be useful, therefore, to test different combinations of covariates to inform their selection \emph{a priori}, i.e. before model fitting, thereby supporting model parsimony.
#' @import ggplot2
#' @param extrapolation.type Character string. Type of extrapolation to be assessed. Can be one of \code{univariate}, \code{combinatorial}, or \code{both} (default).
#' @param n.covariates Maximum number of covariates. The function will compare all combinations of 1 to \code{n.covariates} covariates.
#' @param create.plots Logical, defaults to \code{TRUE}. Whether to produce summary plots.
#' @param display.percent Logical. If \code{TRUE} (default), scales the y-axis of the summary plots as a percentage of the total number of grid cells in \code{prediction.grid}.
#'
#' @inheritParams compute_extrapolation
#'
#' @return Prints a summary table in the R console. Also generates summary boxplots if \code{create.plots} is set to \code{TRUE}.
#'
#' @seealso \code{\link{compute_extrapolation}}, \code{\link{summarise_extrapolation}}
#'
#' @export
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
#' # Define covariates
#' my_cov <- c("Depth", "SST", "NPP", "DistToCAS", "EKE")
#'
#' # Compare the extent of univariate and combinatorial
#' # extrapolation for all combinations of 1 to 5 covariates
#' compare_covariates(extrapolation.type = "both",
#'                   covariate.names = my_cov,
#'                   n.covariates = NULL,
#'                   segments = segs,
#'                   prediction.grid = predgrid,
#'                   coordinate.system = my_crs,
#'                   create.plots = TRUE,
#'                   display.percent = TRUE)
#' @author Phil J. Bouchet
compare_covariates <- function(extrapolation.type = "both",
                               segments,
                               covariate.names,
                               n.covariates = NULL,
                               prediction.grid,
                               coordinate.system,
                               create.plots = TRUE,
                               display.percent = TRUE,
                               resolution = NULL){

  #---------------------------------------------
  # Total number of prediction grid cells
  #---------------------------------------------

  ntot <- nrow(prediction.grid)

  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------

  if(!extrapolation.type%in%c("both", "univariate", "combinatorial"))
    stop("Unknown extrapolation type")

  if(!is.null(n.covariates)){
    if(max(n.covariates) > length(covariate.names))
      stop("n.covariates exceeds the number of covariates available")}

  coordinate.system <- check_crs(coordinate.system = coordinate.system)

  segments <- na.omit(segments)
  prediction.grid <- na.omit(prediction.grid)

  #---------------------------------------------
  # Determine all possible combinations of covariates
  #---------------------------------------------

  # message("Identifying covariate combinations ... ")
  Sys.sleep(time = 0.5)

  if(is.null(n.covariates)){

    combs <- purrr::map2(.x = rep(list(covariate.names), length(covariate.names)),
                         .y = as.list(1:length(covariate.names)),
                         .f = ~utils::combn(x = .x, m = .y, simplify = FALSE)) %>%
      purrr::flatten(.)

  }else if(length(n.covariates)>1){

    combs <- purrr::map2(.x = rep(list(covariate.names), length(n.covariates)),
                         .y = as.list(n.covariates),
                         .f = ~combn(x = .x, m = .y, simplify = FALSE)) %>%
      purrr::flatten(.)

  }else{

    combs <- purrr::map(.x = list(covariate.names),
                        .f = ~combn(x = .x, m = n.covariates, simplify = FALSE)) %>%
      purrr::flatten(.)
  }

  message("Preparing the data ...")

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
  # Carry out extrapolation analysis for each combination of covariates
  #---------------------------------------------

  message("Computing ...")

  pb <- dplyr::progress_estimated(length(combs))

  extrap <- suppressMessages(purrr::map(.x = combs,
                                        .f = ~{
                                          pb$tick()$print()
                                          compute_extrapolation(segments = segments,
                                                                covariate.names = .x,
                                                                prediction.grid = prediction.grid,
                                                                coordinate.system = coordinate.system)},
                                        .pb = pb))

  #---------------------------------------------
  # Summarise extrapolation results for each combination
  #---------------------------------------------

  exsum <- purrr::map(.x = extrap,
                      .f = ~summarise_extrapolation(extrapolation.object = .x,
                                                    extrapolation = FALSE,
                                                    mic = FALSE)) %>%
    purrr::flatten(.)

  extype <- switch(extrapolation.type,
                   "both" = c("univariate", "combinatorial"),
                   "univariate" = "univariate",
                   "combinatorial" = "combinatorial")

  #---------------------------------------------
  # Add zeroes if no extrapolation occurred
  #---------------------------------------------

  for(k in 1:length(exsum)){

    temp <- exsum[k] %>% purrr::flatten(.)
    temp.names <- gsub(pattern = ".n", replacement = "", x = names(temp), fixed = TRUE) %>%
      gsub(pattern = ".p", replacement = "", x = ., fixed = TRUE) %>%
      unique(.)

    if(!"univariate"%in%temp.names) temp <- append(temp, list(univariate.n = 0))
    if(!"combinatorial"%in%temp.names) temp <- append(temp, list(combinatorial.n = 0))

    exsum[[k]] <- temp
  }

  #---------------------------------------------
  # Retrieve numbers of cells
  #---------------------------------------------

  vars <- exsum %>%
    purrr::map(.x = ., .f = ~.x[names(.x)%in%paste0(extype, ".n")]) %>%
    unlist(., use.names = FALSE)

  #---------------------------------------------
  # Build text string of variables
  #---------------------------------------------

  message("\n")
  message("Creating summaries ...")

  if(extrapolation.type=="both"){

    vars.univariate <- vars[seq(1, length(vars), by = 2)]
    vars.combinatorial <- vars[seq(2, length(vars), by = 2)]
    vars.both <- vars.univariate + vars.combinatorial

    #---------------------------------------------
    # Variables minimising/maxisming univariate extrapolation
    #---------------------------------------------

    min.univariate <- vars.univariate[which(vars.univariate==min(vars.univariate))]
    max.univariate <- vars.univariate[which(vars.univariate==max(vars.univariate))]

    varmin.univariate <- combs[which(vars.univariate==min(vars.univariate))] %>%
      collapse.vars(.)

    varmax.univariate <- combs[which(vars.univariate==max(vars.univariate))] %>%
      collapse.vars(.)

    #---------------------------------------------
    # Variables minimising/maxisming combinatorial extrapolation
    #---------------------------------------------

    min.combinatorial <- vars.combinatorial[which(vars.combinatorial==min(vars.combinatorial))]
    max.combinatorial <- vars.combinatorial[which(vars.combinatorial==max(vars.combinatorial))]

    varmin.combinatorial <- combs[which(vars.combinatorial==min(vars.combinatorial))] %>%
      collapse.vars(.)

    varmax.combinatorial <- combs[which(vars.combinatorial==max(vars.combinatorial))] %>%
      collapse.vars(.)

    #---------------------------------------------
    # Variables minimising/maxisming both types of extrapolation
    #---------------------------------------------

    min.both <- vars.both[which(vars.both==min(vars.both))]
    max.both <- vars.both[which(vars.both==max(vars.both))]

    varmin <- combs[which(vars.both==min(vars.both))] %>%
      collapse.vars(.)

    varmax <- combs[which(vars.both==max(vars.both))] %>%
      collapse.vars(.)

  }else{

    #---------------------------------------------
    # Variables minimising/maxisming extrapolation
    #---------------------------------------------

    min.ex <- vars[which(vars==min(vars))]
    max.ex <- vars[which(vars==max(vars))]

    varmin <- combs[which(vars==min(vars))] %>%
      collapse.vars(.)

    varmax <- combs[which(vars==max(vars))] %>%
      collapse.vars(.)

  }

  #---------------------------------------------
  # Manipulate text strings to avoid repetition
  #---------------------------------------------

  if(extrapolation.type=="both"){

    if(!length(varmin.univariate)==length(varmax.univariate)){

      if(length(varmin.univariate)<length(varmax.univariate)){

        varmin.univariate <- c(varmin.univariate, rep("-", times = length(varmax.univariate)-length(varmin.univariate)))
        min.univariate <- c(min.univariate, rep("-", times = length(max.univariate)-length(min.univariate)))

      }else{

        varmax.univariate <- c(varmax.univariate, rep("-", times = length(varmin.univariate)-length(varmax.univariate)))
        max.univariate <- c(max.univariate, rep("-", times = length(min.univariate)-length(max.univariate)))
      }
    }

    if(!length(varmin.combinatorial)==length(varmax.combinatorial)){

      if(length(varmin.combinatorial)<length(varmax.combinatorial)){

        varmin.combinatorial <- c(varmin.combinatorial, rep("-", times = length(varmax.combinatorial)-length(varmin.combinatorial)))
        min.combinatorial <- c(min.combinatorial, rep("-", times = length(max.combinatorial)-length(min.combinatorial)))

      }else{

        varmax.combinatorial <- c(varmax.combinatorial, rep("-", times = length(varmin.combinatorial)-length(varmax.combinatorial)))
        max.combinatorial <- c(max.combinatorial, rep("-", times = length(min.combinatorial)-length(max.combinatorial)))
      }
    }

    if(!length(varmin)==length(varmax)){

      if(length(varmin)<length(varmax)){

        varmin <- c(varmin, rep("-", times = length(varmax)-length(varmin)))
        min.both <- c(min.both, rep("-", times = length(max.both)-length(min.both)))

      }else{

        varmax <- c(varmax, rep("-", times = length(varmin)-length(varmax)))
        max.both <- c(max.both, rep("-", times = length(min.both)-length(max.both)))
      }
    }

    #---------------------------------------------
    # Summarise all results in a table
    #---------------------------------------------

    restxt <- data.frame(Extrapolation = c("Univariate",
                                           rep("", max(length(varmin.univariate), length(varmax.univariate))-1),
                                           "Combinatorial",
                                           rep("", max(length(varmin.combinatorial), length(varmax.combinatorial))-1),
                                           "Both",
                                           rep("", max(length(varmin), length(varmax))-1)))


    restxt$Minimum <-  c(varmin.univariate, varmin.combinatorial, varmin)
    restxt$n_min <-  c(min.univariate, min.combinatorial, min.both)
    restxt$Maximum <-  c(varmax.univariate, varmax.combinatorial, varmax)
    restxt$n_max <-  c(max.univariate, max.combinatorial, max.both)

  }else{

    if(!length(varmin)==length(varmax)){

      if(length(varmin)<length(varmax)){

        varmin <- c(varmin, rep("-", times = length(varmax)-length(varmin)))

      }else{

        varmax <- c(varmax, rep("-", times = length(varmin)-length(varmax)))
      }
    }

    #---------------------------------------------
    # Summarise all results in a table
    #---------------------------------------------

    restxt <- data.frame(Extrapolation = switch(extrapolation.type,
                                                "univariate" = c("Univariate",
                                                                 rep("", max(length(varmin),
                                                                             length(varmax))-1)),
                                                "combinatorial" = c("Combinatorial",
                                                                    rep("", max(length(varmin),
                                                                                length(varmax))-1))
    ))

    restxt$Minimum = varmin
    restxt$n_min = min.ex
    restxt$Maximum = varmax
    restxt$n_max = max.ex
  }

  #---------------------------------------------
  # Visualise results as boxplots
  #---------------------------------------------

  if(create.plots){

    #---------------------------------------------
    # Number of covariate permutations
    #---------------------------------------------

    Ax <- purrr::map_dbl(.x = combs,
                         .f = ~length(.x)) %>%
      table()%>%
      data.frame(.)

    names(Ax) <- c("samp", "Freq")
    Ax$samp <- as.numeric(as.character(Ax$samp))

    #---------------------------------------------
    # Summary by number of variables
    #---------------------------------------------

    if(extrapolation.type=="both"){

      res <- data.frame(nvars = purrr::map_dbl(combs, ~length(.)),
                        extrap = c(vars.univariate, vars.combinatorial),
                        type = rep(c("Univariate", "Combinatorial"),
                                   each = length(vars)/2))
    }else{

      res <- data.frame(nvars = purrr::map_dbl(combs, ~length(.)),
                        extrap = vars)
    }

    if(display.percent) res$extrap <- 100*res$extrap/ntot

    res <- dplyr::left_join(res, Ax, by = c("nvars" = "samp"))

    #---------------------------------------------
    # Compute the number of permutations for each covariate
    #---------------------------------------------

    var.ind <- purrr::map(.x = as.list(covariate.names),
                          .f = ~.x %>%
                            purrr::map2(.x = rep(., length(combs)),
                                        .y = combs,
                                        .f = ~match(.x, .y))) %>%
      purrr::map(.x = ., .f = ~!is.na(.x)) %>%
      purrr::set_names(., covariate.names)

    Bx <- purrr::map_df(var.ind, ~sum(.)) %>%
      tidyr::gather()

    #---------------------------------------------
    # Summary by variable name
    #---------------------------------------------

    if(extrapolation.type=="both"){

      var.ind.univariate <- purrr::map_df(.x = var.ind, .f = ~vars.univariate[.]) %>%
        tidyr::gather(data = .)

      var.ind.combinatorial <- purrr::map_df(.x = var.ind, .f = ~vars.combinatorial[.]) %>%
        tidyr::gather(data = .)

      var.ind.univariate <- var.ind.univariate %>% dplyr::mutate(type = "Univariate")
      var.ind.combinatorial <- var.ind.combinatorial %>% dplyr::mutate(type = "Combinatorial")

      varplot <- data.frame(extrap = rbind(var.ind.univariate, var.ind.combinatorial))

      if(display.percent) varplot$extrap.value <- 100*varplot$extrap.value/ntot

      names(varplot) <- c("var", "extrap", "type")
      varplot <- dplyr::left_join(varplot, Bx, by = c("var" = "key"))

    }else{

      var.ind <- purrr::map(.x = var.ind, .f = ~vars[.])

      varplot <- data.frame(var = rep(covariate.names,
                                      each = length(var.ind[[1]])),
                            extrap = unlist(var.ind))

      if(display.percent) varplot$extrap <- 100*varplot$extrap/ntot

      varplot$var <- as.character(varplot$var)
      varplot <- dplyr::left_join(varplot, Bx, by = c("var" = "key"))

    }

    #---------------------------------------------
    # Plot settings
    #---------------------------------------------

    if(extrapolation.type=="both"){

      colsgroup.a <- c("nvars", "type")
      colsgroup.b <- c("var", "type")

    }else{

      colsgroup.a <- c("nvars")
      colsgroup.b <- c("var")

    }

    #---------------------------------------------
    # Plot labels
    #---------------------------------------------

    labelres <- res %>%
      dplyr::mutate(label = paste0("Nc = ",Freq)) %>%
      dplyr::group_by_at(colsgroup.a) %>%
      dplyr::summarize(ypos = max(res$extrap)+0.15*max(res$extrap),
                       label = unique(label))

    labelvar <- varplot %>%
      dplyr::mutate(label = paste0("Nc = ",value)) %>%
      dplyr::group_by_at(colsgroup.b) %>%
      dplyr::summarize(ypos = max(varplot$extrap)+0.15*max(varplot$extrap),
                       label = unique(label))

    #---------------------------------------------
    # Plots
    #---------------------------------------------

    p1 <- ggplot2::ggplot(data = res, ggplot2::aes(x = factor(nvars), y = extrap))+

      {if(extrapolation.type=="both") geom_boxplot(aes(fill=type))}+
      {if(extrapolation.type=="univariate") geom_boxplot(fill='#F0B039')}+
      {if(extrapolation.type=="combinatorial") geom_boxplot(fill='#2E86A6')}+

      geom_text(data = labelres, ggplot2::aes(x = factor(nvars), y = ypos, label=label))+
      xlab("Number of variables")+
      {if(display.percent) ylab("Extrapolation (%)")}+
      {if(display.percent==FALSE) ylab("Extrapolation")}+
      ylim(c(0, max(res$extrap)+0.15*max(res$extrap)))+

      theme(legend.title = element_blank())+
      theme_minimal()+
      {if(extrapolation.type=="both") scale_fill_manual(values = c('#2E86A6', '#F0B039'))}

    p2 <- ggplot2::ggplot(data = varplot, ggplot2::aes(x = factor(var), y = extrap))+

      {if(extrapolation.type=="both") geom_boxplot(aes(fill=type))}+
      {if(extrapolation.type=="univariate") geom_boxplot(fill='#F0B039')}+
      {if(extrapolation.type=="combinatorial") geom_boxplot(fill='#2E86A6')}+

      geom_text(data = labelvar, ggplot2::aes(x = factor(var), y = ypos, label=label))+
      xlab("Variable")+
      {if(display.percent) ylab("Extrapolation (%)")}+
      {if(display.percent==FALSE) ylab("Extrapolation")}+
      ylim(c(0, max(varplot$extrap)+0.15*max(varplot$extrap)))+

      ggplot2::theme(legend.title = element_blank())+
      theme_minimal()+
      {if(extrapolation.type=="both") scale_fill_manual(values = c('#2E86A6', '#F0B039'))}

    p3 <- cowplot::ggdraw() +
      cowplot::draw_plot(p1, x = 0, y = 0.5, width = 1, height = 0.5) +
      cowplot::draw_plot(p2, x = 0, y = 0, width = 1, height = 0.5)

    print(p3)
  }

  message("Done!")

  print(knitr::kable(restxt, format = "pandoc"))

}
