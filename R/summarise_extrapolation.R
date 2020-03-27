#' Summary of extrapolation
#'
#' Displays a tabular summary of the geographic extent of extrapolation. This is calculated as the number (and proportion) of prediction locations (i.e. grid cells) subject to extrapolation.
#'
#' @param extrapolation.object Output object from a run of \link{compute_extrapolation}.
#' @param extrapolation Logical. Whether to return a summary of univariate/combinatorial extrapolation. Defaults to TRUE.
#' @param mic Logical. Whether to return a summary of the most influential covariates (MIC) - see \link{compute_extrapolation}. Defaults to TRUE.
#' @inheritParams compute_extrapolation
#' @return Prints a summary table in the R console. In addition, if assigned to an object, returns a list with the table values (.n = number of locations, .p = corresponding percentage).
#' @author Phil J. Bouchet
#' @seealso \code{\link{compute_extrapolation}}
#' @references Bouchet PJ, Miller DL, Roberts JJ, Mannocci L, Harris CM and Thomas L (2019). From here and now to there and then: Practical recommendations for extrapolating cetacean density surface models to novel conditions. CREEM Technical Report 2019-01, 59 p. \href{https://research-repository.st-andrews.ac.uk/handle/10023/18509}{https://research-repository.st-andrews.ac.uk/handle/10023/18509}
#'
#' Mesgaran MB, Cousens RD, Webber BL (2014). Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity & Distributions, 20: 1147-1159, DOI: \href{https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209}{10.1111/ddi.12209}
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
#' # Assess extrapolation in the multivariate space defined by five covariates
#' spermw.extrapolation <- compute_extrapolation(segments = segs,
#'       covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
#'       prediction.grid = predgrid,
#'       coordinate.system = my_crs,
#'       print.summary = FALSE,
#'       save.summary = TRUE,
#'       print.precision = 2)
#'
#' # Summarise extrapolation
#' spermw.summary <- summarise_extrapolation(extrapolation.object = spermw.extrapolation,
#'      covariate.names = c("Depth", "DistToCAS", "SST", "EKE", "NPP"),
#'      extrapolation = TRUE,
#'      mic = TRUE,
#'      print.precision = 2)
#'
#' print(spermw.summary)

summarise_extrapolation <- function(extrapolation.object,
                                    covariate.names = NULL,
                                    extrapolation = TRUE,
                                    mic = TRUE,
                                    print.precision = 2){

  #.............................................
  # Extract extrapolation values
  #---------------------------------------------

  ex.data <- extrapolation.object$data$all$ExDet

  #---------------------------------------------
  # Extract output values
  #---------------------------------------------

  res <- n_and_p(ex.data)

  #---------------------------------------------
  # Format as nice-looking table
  #---------------------------------------------

  # Which type(s) extrapolation did not occur?

  zeroes <- purrr::map_lgl(.x = res, .f = ~.x==0)
  zeroes.names <- names(zeroes[!zeroes])

  tb.names <- gsub(pattern = ".n", replacement = "", x = zeroes.names, fixed = TRUE) %>%
    gsub(pattern = ".p", replacement = "", x = ., fixed = TRUE) %>%
    unique(.) %>%
    tools::toTitleCase(.)

  #---------------------------------------------
  # Filter accordingly
  #---------------------------------------------

  res <- res[zeroes.names]

  #---------------------------------------------
  # Begin formatting
  #---------------------------------------------

  resdf <- purrr::map2(.x = res[grepl(pattern = ".n", x = names(res), fixed = TRUE)],
                       .y = res[grepl(pattern = ".p", x = names(res), fixed = TRUE)],
                       .f = ~c(paste0("n = ", format(round(as.numeric(.x), 0),
                                                     nsmall=0, big.mark=",")), paste0(round(.y, print.precision), " %")))

  #---------------------------------------------
  # Convert to matrix form
  #---------------------------------------------

  resdf <- purrr::set_names(resdf, tb.names) %>%
    data.frame(.) %>%
    t(.)

  resdf <- cbind(resdf, row.names(resdf))

  #---------------------------------------------
  # Calculate totals
  #---------------------------------------------

  Total.n <- res[grepl(pattern = ".n", x = names(res), fixed = TRUE)] %>%
    unlist(.) %>%
    sum(.) %>%
    format(round(., 0), nsmall=0, big.mark=",") %>%
    paste0("n = ", .)

  Total.p <- res[grepl(pattern = ".p", x = names(res), fixed = TRUE)] %>%
    unlist(.) %>%
    sum(.) %>%
    round(., print.precision) %>%
    paste0(., " %")

  resdf <- data.frame(resdf)
  names(resdf) <- c("Count", "Percentage", "Type")
  row.names(resdf) <- NULL

  #---------------------------------------------
  # Rearrange rows and reorder columns
  #---------------------------------------------

  resdf <- resdf %>%
    dplyr::mutate(., Type = factor(Type,
                                   levels = c("Analogue", "Univariate", "Combinatorial", "-----------", "Total"))) %>%
    dplyr::arrange(., Type) %>%
    dplyr::select(., Type, Count, Percentage)


  #---------------------------------------------
  # Add separator beneath "Analogue" if present
  #---------------------------------------------

  if("Analogue"%in%as.character(resdf$Type)) add.sep <- TRUE else add.sep <- FALSE

  resdf <- purrr::map_dfr(resdf, as.character)

  #---------------------------------------------
  # Add sub-totals, if necessary
  #---------------------------------------------

  if("Univariate"%in%resdf[,1] & "Combinatorial"%in%resdf[,1]){

    Totalex.n <- res[!grepl(pattern = "analogue", x = names(res))] %>%
      .[grepl(pattern = ".n", x = names(.), fixed = TRUE)] %>%
      unlist(.) %>%
      sum(.) %>%
      format(round(., 0), nsmall=0, big.mark=",") %>%
      paste0("n = ", .)

    Totalex.p <- res[!grepl(pattern = "analogue", x = names(res))] %>%
      .[grepl(pattern = ".p", x = names(.), fixed = TRUE)] %>%
      unlist(.) %>%
      sum(.) %>%
      round(., print.precision) %>%
      paste0(., " %")

    resdf <- rbind(resdf, rep("-----------",3))
    resdf <- rbind(resdf, c("  Sub-total", Totalex.n, Totalex.p))

  }

  #---------------------------------------------
  # Add totals to matrix
  #---------------------------------------------

  resdf <- rbind(resdf, rep("-----------",3))
  resdf <- rbind(resdf, c("Total", Total.n, Total.p))

  if(add.sep) resdf <- rbind(resdf[1,], rep("-----------",3), resdf[2:nrow(resdf),])

  colnames(resdf) <- NULL

  #---------------------------------------------
  # Most influential covariates - by extrapolation type
  #---------------------------------------------

  if(mic){

    mic_data_univariate <- covariate.names[extrapolation.object$data$all$mic_univariate]
    mic_data_combinatorial <- covariate.names[extrapolation.object$data$all$mic_combinatorial]

    #---------------------------------------------
    # Tabulate the results
    #---------------------------------------------

    mic_data <- list()

    #---------------------------------------------
    # Univariate extrapolation
    #---------------------------------------------

    if(all(is.na(mic_data_univariate))){

      mic_data_univariate <- list()

    }else{

      mic_data_univariate <- mic_data_univariate %>%
        table(.) %>%
        as.data.frame(.) %>%
        dplyr::mutate(type = "Univariate") %>%
        dplyr::mutate(perc = 100*Freq/length(extrapolation.object$data$all$ExDet))

    }

    #---------------------------------------------
    # Combinatorial extrapolation
    #---------------------------------------------

    if(all(is.na(mic_data_combinatorial))){

      mic_data_combinatorial <- character(0)

    }else{

      mic_data_combinatorial <- mic_data_combinatorial %>%
        table(.) %>%
        as.data.frame(.) %>%
        dplyr::mutate(type = "Combinatorial") %>%
        dplyr::mutate(perc = 100*Freq/length(extrapolation.object$data$all$ExDet))
    }

    #---------------------------------------------
    # Combine into single list
    #---------------------------------------------

    if(!purrr::is_empty(mic_data_univariate)) mic_data <- append(mic_data, list(mic_data_univariate))
    if(!purrr::is_empty(mic_data_combinatorial)) mic_data <- append(mic_data, list(mic_data_combinatorial))

    #---------------------------------------------
    # Rename columns
    #---------------------------------------------

    if(!purrr::is_empty(mic_data)){

      mic_data <- purrr::map(.x = mic_data,
                             .f = ~purrr::set_names(.x, c("covariate", "freq", "Type", "perc")))

      #---------------------------------------------
      # Format list data
      #---------------------------------------------

      mic_res <- purrr::map(.x = mic_data,
                            .f = ~strsplit(paste(.x$freq, .x$perc), " ")) %>%
        purrr::map2(.x = ., .y = mic_data, .f = ~set_names(.x, sort(.y$covariate))) %>%
        purrr::map(.x = ., .f = ~purrr::map(.x = ., .f = ~as.numeric(.x))) %>%
        purrr::map(.x = ., .f = ~purrr::map(.x = ., .f = ~list(.n = .x[1], .p = .x[2]))) %>%
        purrr::flatten()

      #---------------------------------------------
      # Format for output in console
      #---------------------------------------------

      mic_resdf <- purrr::map(.x = mic_data,
                              .f = ~dplyr::arrange(.x, desc(perc), covariate) %>%
                                dplyr::mutate(freq = paste0("n = ", format(round(as.numeric(freq),
                                                                                 print.precision), nsmall=0, big.mark=","))) %>%
                                dplyr::mutate(perc = paste0(round(as.numeric(perc), print.precision), " %")))

      mic_data <- purrr::set_names(mic_data, tb.names[which(!tb.names=="Analogue")])
      mic_resdf <- purrr::set_names(mic_resdf, tb.names[which(!tb.names=="Analogue")])

      #---------------------------------------------
      # Compact into one tibble
      #---------------------------------------------

      mic_resdf <- do.call(rbind, mic_resdf)

      #---------------------------------------------
      # Re-order columns
      #---------------------------------------------

      mic_resdf <- mic_resdf %>%
        dplyr::select(., Type, covariate, freq, perc)

      row.names(mic_resdf) <- NULL

      #---------------------------------------------
      # Calculate totals
      #---------------------------------------------

      mic_data_total <- do.call(rbind, mic_data)
      mic_subtotals <- mic_data_total %>%
        dplyr::group_by(Type) %>%
        dplyr::summarise(sum(freq), sum(perc))

      names(mic_subtotals) <- c("Type", "freq", "perc")


      if("Univariate"%in%as.character(mic_resdf$Type) & "Combinatorial"%in%as.character(mic_resdf$Type)){

        sub.univ <- as.matrix(mic_subtotals[mic_subtotals$Type=="Univariate",])
        sub.univ <- c("  Sub-total", "", paste0("n = ", format(round(as.numeric(sub.univ[2]),
                                                                     print.precision),
                                                               nsmall=0, big.mark=",")),
                      paste0(round(as.numeric(sub.univ[3]), print.precision), " %"))
        sub.univ <- data.frame(t(matrix(sub.univ)))
        names(sub.univ) <- names(mic_resdf)

        sub.comb <- as.matrix(mic_subtotals[mic_subtotals$Type=="Combinatorial",])
        sub.comb <- c("  Sub-total", "", paste0("n = ", format(round(as.numeric(sub.comb[2]),
                                                                     print.precision),
                                                               nsmall=0, big.mark=",")),
                      paste0(round(as.numeric(sub.comb[3]), print.precision), " %"))
        sub.comb <- data.frame(t(matrix(sub.comb)))
        names(sub.comb) <- names(mic_resdf)

        mic_resdf <- rbind(mic_resdf[mic_resdf$Type=="Univariate",],
                           data.frame(Type = "-----------", covariate = "-----------", freq = "-----------", perc = "-----------"),
                           sub.univ,
                           data.frame(Type = "-----------", covariate = "-----------", freq = "-----------", perc = "-----------"),
                           mic_resdf[mic_resdf$Type=="Combinatorial",],
                           data.frame(Type = "-----------", covariate = "-----------", freq = "-----------", perc = "-----------"),
                           sub.comb)

      }

      mic_resdf <- purrr::map_dfr(mic_resdf, as.character)

      mic_resdf <- rbind(mic_resdf, rep("-----------", ncol(mic_resdf)))
      mic_resdf <- rbind(mic_resdf, c("Total",
                                      "",
                                      paste0("n = ", format(round(sum(mic_data_total$freq),
                                                                  print.precision),
                                                            nsmall=0, big.mark=",")),
                                      paste0(round(as.numeric(sum(mic_data_total$perc)), print.precision), " %")))


      colnames(mic_resdf) <- NULL
      rownames(mic_resdf) <- NULL

      # If all covariates contribute equally, then use alphabetical order

      # cols_ordered <- purrr::map(.x = mic_resdf, ~c("Type", as.character(.x$covariate)))
      #
      # mic_resdf <- purrr::map2(.x = mic_resdf,
      #                          .y = cols_ordered,
      #                          .f = ~tidyr::gather(data = .x, key = item,
      #                                              value = result, -Type, -covariate) %>%
      #                            tidyr::spread(data = ., key = covariate, value = result) %>%
      #                            .[, .y])

    } # End empty mic_data

  } # End mic


  #---------------------------------------------
  # Print in console
  #---------------------------------------------

  if(extrapolation){
    print(knitr::kable(resdf,
                       format = "pandoc",
                       caption = "Extrapolation"))}

  if(mic){
    print(knitr::kable(mic_resdf,
                       format = "pandoc",
                       caption = "Most influential covariates (MIC)"))
  }

  #---------------------------------------------
  # Return output
  #---------------------------------------------

  if(mic){

    invisible(list(extrapolation = res,
                   mic = mic_res))
  }else{

    invisible(list(extrapolation = res))
  }





} # End summarise_extrapolation
