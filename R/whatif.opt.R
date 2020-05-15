#' Counterfactual evaluation, optimised for dsmextra
#'
#' Implements the methods described in King and Zeng (2007) for evaluating counterfactuals.
#'
#' The \code{\link[WhatIf]{whatif}} function from the \href{https://CRAN.R-project.org/package=WhatIf}{Whatif} package (Gandrud et al. 2017) may not run on very large datasets. To circumvent this problem, \code{whatif.opt} sets the calculations performed by \code{\link[WhatIf]{whatif}} to run on partitions of the data instead, for greater efficiency. \code{whatif.opt} can be called internally within compute_nearby by using two additional arguments, namely:
#' \tabular{ll}{
#'   \code{max.size} \tab Threshold above which partitioning will be triggered \cr
#'   \code{no.partitions} \tab Number of required partitions \cr
#'  }
#'  In practice, a run of \code{\link{compute_nearby}} begins with a quick assessment of the dimensions of the input data, i.e. the reference and target data.frames. If the product of their dimensions (i.e. number of segments multiplied by number of prediction grid cells) exceeds the value set for \code{max.size}, then \code{no.partitions} subsets of the data will be created and the computations run on each using map functions from the \code{\link{purrr}} package (Henry and Wickham 2019). This means that a smaller \code{max.size} will trigger partitioning on correspondingly smaller datasets. By default, \code{max.size} is set to \code{1e7}. This value was chosen arbitrarily, and should be sufficiently large as to obviate the need for partitioning on most datasets.
#'
#' @param formula An optional formula without a dependent variable, allowing transformations of combinations of the variables in both data and cfact.
#' @param data Either a model output object, or a n-by-k non-character (logical or numeric) matrix or data frame of observed covariate data with n data points or units and k covariates.
#' @param cfact Counterfactuals.
#' @param nearby Scalar indicating which reference data points are considered to be 'nearby' (i.e. within ‘nearby’ mean geometric Gower's distance of) prediction points.
#' @param miss Optional string indicating the strategy for dealing with missing data.
#' @param no.partitions Integer. Number of desired partitions of the data (default of 10).
#'
#' @importFrom stats model.frame model.matrix na.fail na.omit terms quantile update.formula
#'
#' @return A list object containing extrapolation values in both data.frame and raster format.
#'
#' @author Phil J. Bouchet
#'
#' @seealso \code{\link{compute_nearby}}, \code{\link[WhatIf]{whatif}}
#'
#' @references Gandrud C, King G, Stoll H, Zeng L (2017). WhatIf: Evaluate Counterfactuals. R package version 1.5-9. \href{https://CRAN.R-project.org/package=WhatIf}{https://CRAN.R-project.org/package=WhatIf}.
#'
#' King G, Zeng L (2007). When can history be our guide? The pitfalls of counterfactual inference. International Studies Quarterly 51, 183–210. DOI: \href{https://www.jstor.org/stable/pdf/4621707.pdf?seq=1#page_scan_tab_contents}{10.1111/j.1468-2478.2007.00445.x}
#'
#' Henry L, Wickham H (2019). purrr: Functional Programming Tools. R package version 0.3.2. \href{https://CRAN.R-project.org/package=purrr}{https://CRAN.R-project.org/package=purrr}.
#'
whatif.opt <- function (formula = NULL,
                        data, cfact,
                        nearby = 1,
                        miss = "list",
                        no.partitions)
{

  message("Preprocessing data ...")

  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------

  # if (grepl("Zelig*", class(data)) & missing(cfact))
  #   cfact <- Zelig::zelig_setx_to_df(data)
  #
  # if (grepl("Zelig*", class(data)) & !missing(cfact)) {
  #   formula <- formula(stats::delete.response(stats::terms(data$formula)))
  #   data <- data$zelig.out$z.out[[1]]$model
  # }

  if (!((is.character(cfact) && is.vector(cfact) && length(cfact) ==
         1) || is.data.frame(cfact) || (is.matrix(cfact) && !is.character(cfact)))) {
    stop("'cfact' must be either a string, a R data frame, or a R non-character matrix")
  }

  if (is.character(cfact)) {
    cfact <- utils::read.table(cfact)
  }

  if (dim(cfact)[1] == 0) {
    stop("no counterfactuals supplied: 'cfact' contains zero rows")
  }

  if (!any(stats::complete.cases(cfact))) {
    stop("there are no cases in 'cfact' without missing values")
  }

  if ("(Intercept)" %in% dimnames(cfact)[[2]]) {
    cfact <- cfact[, -(which(dimnames(cfact)[[2]] == "(Intercept)"))]
  }

  if (is.list(data) && !(is.data.frame(data))) {
    if (!((("formula" %in% names(data)) || ("terms" %in%
                                            names(data))) && (("data" %in% names(data)) || ("model" %in%
                                                                                            names(data))))) {
      stop("the list supplied to 'data' is not a valid output object")
    }

    tt <- terms(data)
    attr(tt, "intercept") <- rep(0, length(attr(tt, "intercept")))
    if ("data" %in% names(data)) {
      if (is.data.frame(data$data)) {
        data <- model.matrix(tt, model.frame(tt, data = data$data,
                                             na.action = NULL))
      }else {
        data <- model.matrix(tt, model.frame(tt, data = eval(data$data,
                                                             envir = .GlobalEnv), na.action = NULL))
      }
    }else {
      data <- model.matrix(tt, data = data$model)
    }
    if (!(is.matrix(data))) {
      stop("observed covariate data could not be extracted from output object")
    }
    rm(tt)
  }else {
    if (!((is.character(data) && is.vector(data) && length(data) ==
           1) || is.data.frame(data) || (is.matrix(data) &&
                                         !is.character(data)))) {
      stop("'data' must be either a string, a R data frame, a R non-character matrix, or an output object")
    }
    if (is.character(data)) {
      data <- utils::read.table(data)
    }
  }

  if (dim(data)[1] == 0) {
    stop("no observed covariate data supplied: 'data' contains zero rows")
  }

  if (!any(stats::complete.cases(data))) {
    stop("there are no cases in 'data' without missing values")
  }

  if (!(is.null(formula))) {
    if (identical(class(formula), "formula")) {
      if (!(is.data.frame(as.data.frame(data)))) {
        stop("'data' must be coercable to a data frame in order to use 'formula'")
      }
      if (!(is.data.frame(as.data.frame(cfact)))) {
        stop("'cfact' must be coercable to a data frame in order to use 'formula'")
      }
      formula <- update.formula(formula, ~. - 1)
      ttvar <- all.vars(formula)
      for (i in 1:length(ttvar)) {
        if (!(ttvar[i] %in% dimnames(data)[[2]])) {
          stop("variables in 'formula' either unlabeled or not present in 'data'")
        }
        if (!(ttvar[i] %in% dimnames(cfact)[[2]])) {
          stop("variable(s) in 'formula' either unlabeled or not present in 'cfact'")
        }
      }
      rm(ttvar)
      data <- model.matrix(formula, data = model.frame(formula,
                                                       as.data.frame(data), na.action = NULL))
      cfact <- model.matrix(formula, data = model.frame(formula,
                                                        as.data.frame(cfact), na.action = NULL))
    }else {
      stop("'formula' must be of class 'formula'")
    }
  }

  if (!(identical(stats::complete.cases(cfact), rep(TRUE, dim(cfact)[1])))) {
    cfact <- na.omit(cfact)
    message("Note:  counterfactuals with missing values eliminated from cfact")
  }

  if (is.data.frame(data)) {
    if (is.character(as.matrix(data))) {
      stop("observed covariate data not coercable to numeric matrix due to character column(s)")
    }
    data <- suppressWarnings(data.matrix(data))
  }else {
    data <- data.matrix(as.data.frame(data))
  }

  if (is.data.frame(cfact)) {
    if (is.character(as.matrix(cfact))) {
      stop("counterfactual data not coercable to numeric matrix due to character column(s)")
    }
    cfact <- suppressWarnings(data.matrix(cfact))
  }else{
    cfact <- data.matrix(as.data.frame(cfact))
  }

  if (!(is.matrix(data) && is.numeric(data))) {
    stop("observed covariate data not coercable to numeric matrix")
  }

  if (!(is.matrix(cfact) && is.numeric(cfact))) {
    stop("counterfactual data not coercable to numeric matrix")
  }
  na.fail(cfact)

  if (!identical(ncol(cfact), ncol(data))) {
    stop("number of columns of 'cfact' and 'data' are not equal")
  }


  if (!(is.null(nearby))) {
    if (!(is.numeric(nearby) && is.vector(nearby) && length(nearby) ==
          1 && nearby >= 0)) {
      stop("'nearby' must be numeric, greater than or equal to 0, and a scalar")
    }
  }

  if (!(identical(miss, "list") || identical(miss, "case"))) {
    stop("'miss' must be either ''case'' or ''list''")
  }

  n = nrow(data)

  #---------------------------------------------
  # Define functions
  #---------------------------------------------

  # Original functions

  calc.gd <- function(dat, cf, range) {
    n <- nrow(dat)
    m <- nrow(cf)
    dat = t(dat)
    dist = matrix(0, m, n, dimnames = list(1:m, 1:n))
    for (i in 1:m) {
      temp <- abs(dat - cf[i, ])/range
      if (any(range == 0)) {
        temp[is.nan(temp)] <- 0
        temp[temp == Inf] <- NA
      }
      dist[i, ] <- colMeans(temp, na.rm = T)
    }
    return(t(dist))
  }

  geom.var <- function(dat, rang) {
    n <- nrow(dat)
    dat <- t(dat)
    ff <- function(x) {
      temp <- abs(dat - x)/rang
      if (any(rang == 0)) {
        temp[is.nan(temp)] <- 0
        temp[temp == Inf] <- NA
      }
      tmp <- sum(colMeans(temp, na.rm = TRUE))
      return(tmp)
    }
    sum.gd.x <- sum(apply(dat, 2, ff), na.rm = TRUE)
    gv.x <- (0.5 * sum.gd.x)/(n^2)
    return(gv.x)
  }


  calcgd <- function(dat, cf, range, split.factor = no.partitions) {

    # Split matrices into smaller chunks

    nlist <- split(1:nrow(dat),
                   cut(seq_along(1:nrow(dat)),
                       split.factor, labels = FALSE))
    mlist <- split(1:nrow(cf),
                   cut(seq_along(1:nrow(cf)),
                       split.factor, labels = FALSE))

    chunkdat <- purrr::map(.x = nlist, .f = ~dat[.x,])
    chunkcf <- purrr::map(.x = mlist, .f = ~cf[.x,])

    # split segments then rbind
    # split predgrid then cbind

    pb <- dplyr::progress_estimated(split.factor, 0) # Progress bar

    chunk.results <- purrr::map(.x = chunkdat,
                                .f = function(x) {
                                  pb$tick()$print()
                                  purrr::map(.x = chunkcf,
                                             function(y) calc.gd(dat = x, cf = y, range = range)) %>%
                                    do.call(cbind, .)})

    chunk.results <- do.call(rbind, chunk.results)
    return(chunk.results)

  } # End calc.gd

  geomvar <- function(dat, rang) {

    n <- nrow(dat)
    dat <- t(dat)

    pbb <- dplyr::progress_estimated(ncol(dat))
    # assign(x = 'pbb', value = dplyr::progress_estimated(ncol(dat)), envir = .GlobalEnv)

    fff <- function(x, dat, rang){
      # pbb$tick()$print()
      return(colMeans(abs(dat - dat[,x])/rang))
    }

    temp <- purrr::map(1:ncol(dat),
                       ~{pbb$tick()$print()
                         fff(x = .x, dat = dat, rang = rang) %>%
                         sum(.)})

    temp <- Reduce('+', temp)
    gv.x <- (0.5 * temp)/(n^2)

  }

  if (identical(miss, "list")) {
    data <- na.omit(data)
    n <- nrow(data)
  }

  #---------------------------------------------
  # Perform calculations
  #---------------------------------------------

  message("Calculating distances ....")

  samp.range <- apply(data, 2, max, na.rm = TRUE) - apply(data, 2, min, na.rm = TRUE)


  if (identical(TRUE, any(samp.range == 0))) {
    message("Note:  range of at least one variable equals zero")
  }

  dist <- calcgd(dat = data,
                 cf = cfact,
                 range = samp.range,
                 split.factor = no.partitions)

  gc()

  message("\n")
  message("Calculating the geometric variance ...")

  gv.x <- geomvar(dat = data, rang = samp.range)

  gc()

  summary <- colSums(dist <= nearby * gv.x) * (1/n)

  #---------------------------------------------
  # Wrap up
  #---------------------------------------------

  message("\n")
  message("Finishing up ...")

  out <- list(call = match.call(), geom.var = gv.x,
              sum.stat = summary)


  class(out) <- "whatif"
  return(invisible(out))
}
