#' Print extrapolation results summary
#'
#' Print the summary of extrapolation results (already calculated by \code{\link{compute_extapolation}}).
#'
#' @export
#' @author David L Miller
#' @param x \code{summary} element of the result of \code{\link{compute_extapolation}}
#' @param digits precision of results
#' @param \dots for S3 compatability
#' @return invisibly returns the summary part of the object only, printing the results
print.extrapolation_results_summary <- function(x, digits=2, ...){

  class(x) <- class(x)[-1]

  if(!is.null(x$extrapolation)){
    #---------------------------------------------
    # Begin formatting
    #---------------------------------------------

    res <- x$extrapolation
    zeroes <- purrr::map_lgl(.x = res, .f = ~.x==0)
    zeroes.names <- names(zeroes[!zeroes])

    tb.names <- gsub(pattern = ".n", replacement = "", x = zeroes.names,
                     fixed = TRUE) %>%
      gsub(pattern = ".p", replacement = "", x = ., fixed = TRUE) %>%
      unique(.) %>%
      tools::toTitleCase(.)

    resdf <- purrr::map2(.x = res[grepl(pattern = ".n", x = names(res),
                                  fixed = TRUE)],
                         .y = res[grepl(pattern = ".p", x = names(res),
                                  fixed = TRUE)],
                         .f = ~c(.x, .y))

    #---------------------------------------------
    # Convert to data.frame form
    #---------------------------------------------

    resdf <- t(data.frame(resdf))
    resdf <- signif(resdf, digits)
    row.names(resdf) <- tb.names

    resdf <- cbind.data.frame(row.names(resdf), resdf)
    names(resdf) <- c("Type", "Count", "Percentage")

    #---------------------------------------------
    # Calculate totals
    #---------------------------------------------

    Total.n <- res[grepl(pattern = ".n", x = names(res), fixed = TRUE)] %>%
      unlist(.) %>%
      sum(.)

    Total.p <- res[grepl(pattern = ".p", x = names(res), fixed = TRUE)] %>%
      unlist(.) %>%
      sum(.)

    resdf <- data.frame(resdf, stringsAsFactors=FALSE)

    #---------------------------------------------
    # Add separator beneath "Analogue" if present
    #---------------------------------------------

    if("Analogue" %in% as.character(resdf$Type)){
      add.sep <- TRUE
    }else{
      add.sep <- FALSE
    }

    resdf <- purrr::map_dfr(resdf, as.character)

    #---------------------------------------------
    # Add sub-totals, if necessary
    #---------------------------------------------

    if("Univariate" %in% resdf[,1] &
       "Combinatorial" %in% resdf[,1]){

      Totalex.n <- res[!grepl(pattern = "analogue", x = names(res))] %>%
        .[grepl(pattern = ".n", x = names(.), fixed = TRUE)] %>%
        unlist(.) %>%
        sum(.)

      Totalex.p <- res[!grepl(pattern = "analogue", x = names(res))] %>%
        .[grepl(pattern = ".p", x = names(.), fixed = TRUE)] %>%
        unlist(.) %>%
        sum(.)

      resdf <- rbind(resdf, rep("-----------",3))
      resdf <- rbind(resdf, c("  Sub-total", Totalex.n, Totalex.p))

    }

    #---------------------------------------------
    # Add totals to matrix
    #---------------------------------------------

    resdf <- rbind(resdf, rep("-----------",3))
    resdf <- rbind(resdf, c("Total", Total.n, Total.p))

    if(add.sep) resdf <- rbind(resdf[1,], rep("-----------",3), resdf[2:nrow(resdf),])

    print(knitr::kable(resdf,
                       format = "pandoc",
                       caption = "Extrapolation"))
  }

  if(!is.null(x$mic)){
    mic_data <- x$mic
    #---------------------------------------------
    # Format for output in console
    #---------------------------------------------

    mic_resdf <- purrr::map(.x = mic_data,
                            .f = ~dplyr::arrange(.x, desc(perc), covariate))
    mic_data <- purrr::set_names(mic_data,
                                 tb.names[which(!tb.names=="Analogue")])
    mic_resdf <- purrr::set_names(mic_resdf,
                                  tb.names[which(!tb.names=="Analogue")])

    #---------------------------------------------
    # Compact into one tibble
    #---------------------------------------------

    mic_resdf <- do.call(rbind, mic_resdf)

    #---------------------------------------------
    # Re-order columns
    #---------------------------------------------

    mic_resdf <- mic_resdf %>%
      dplyr::mutate(perc = signif(perc, digits)) %>%
      dplyr::select(., Type, covariate, freq, perc)

    row.names(mic_resdf) <- NULL

    #---------------------------------------------
    # Calculate totals
    #---------------------------------------------

    mic_data_total <- do.call(rbind, mic_data)
    mic_subtotals <- mic_data_total %>%
      dplyr::group_by(Type) %>%
      dplyr::summarise(sum(freq), signif(sum(perc), digits))

    names(mic_subtotals) <- c("Type", "freq", "perc")


    if("Univariate" %in% as.character(mic_resdf$Type) &
       "Combinatorial" %in% as.character(mic_resdf$Type)){

      sub.univ <- as.matrix(mic_subtotals[mic_subtotals$Type=="Univariate",])
      sub.univ <- c("  Sub-total", "", sub.univ[2], sub.univ[3])
      sub.univ <- data.frame(t(matrix(sub.univ)))
      names(sub.univ) <- names(mic_resdf)

      sub.comb <- as.matrix(mic_subtotals[mic_subtotals$Type=="Combinatorial",])
      sub.comb <- c("  Sub-total", "", sub.comb[2], sub.comb[3])
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
                                    sum(mic_data_total$freq),
                                    signif(sum(mic_data_total$perc), digits)))

    # get column names
    colnames(mic_resdf) <- c("Type", "Covariate", "Count", "Percentage")
    rownames(mic_resdf) <- NULL

    print(knitr::kable(mic_resdf,
                       format = "pandoc",
                       caption = "Most influential covariates (MIC)"))
  }

}
