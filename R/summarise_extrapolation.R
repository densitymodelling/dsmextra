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
#' @keywords internal

summarise_extrapolation <- function(extrapolation.object,
                                    covariate.names = NULL,
                                    extrapolation = TRUE,
                                    mic = TRUE){
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
    }

  }
  #---------------------------------------------
  # Return output
  #---------------------------------------------

  if(mic){

    invisible(list(extrapolation = res,
                   mic = mic_data))
  }else{

    invisible(list(extrapolation = res))
  }


}
