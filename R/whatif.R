#' Counterfactual evaluation
#'
#' Implements the methods described in King and Zeng (2007) for evaluating counterfactuals. This function is a duplicate of the \code{whatif} function from the \href{https://CRAN.R-project.org/package=WhatIf}{Whatif} package (GPL >= 3 license).
#'
#' @param formula An optional formula without a dependent variable that is of class "formula" and that follows standard \code{R} conventions for formulas, e.g. ~ x1 + x2.  Allows you to transform or otherwise re-specify combinations of the variables in both \code{data} and \code{cfact}. To use this parameter, both \code{data} and \code{cfact} must be coercable to a \code{data.frame}; the variables of both \code{data} and \code{cfact} must be labelled; and all variables appearing in \code{formula} must also appear in both \code{data} and \code{cfact}. Otherwise, errors are returned. The intercept is automatically dropped. Default is \code{NULL}.
#' @param data Either a model output object, or a n-by-k non-character (logical or numeric) matrix or data frame of observed covariate data with n data points or units and k covariates.
#' @param cfact A \code{R} object or a string.  If a \code{R} object, a \eqn{m}-by-\eqn{k} non-character matrix or data frame of \emph{counterfactuals} with \eqn{m} counterfactuals and the same \eqn{k} covariates (in the same order) as in \code{data}.  However, if \code{formula} is used to select a subset of the \eqn{k} covariates, then \code{cfact} may contain either only these \eqn{j \leq k}{j <= K} covariates or the complete set of \eqn{k} covariates. An intercept should not be included as one of the covariates. Data frames will again be coerced to their internal numeric values if possible. If a string, either the complete path (including file name) of the file containing the counterfactuals or the path relative to your working directory. This file should be a white space delimited text file.
#' @param nearby An optional scalar indicating which observed data points are considered to be nearby (i.e., withing ‘nearby’ geometric variances of) the counterfactuals. Used to calculate the summary statistic returned by the function: the fraction of the observed data nearby each counterfactual. By default, the geometric variance of the covariate data is used. For example, setting nearby to 2 will identify the proportion of data points within two geometric variances of a counterfactual. Default is \code{NULL}.
#' @param choice An optional string indicating which analyses to undertake. The options are either "hull", only perform the convex hull membership test; "distance", do not perform the convex hull test but do everything else, such as calculating the distance between each counterfactual and data point; or "both", undertake both the convex hull test and the distance calculations (i.e., do everything). Default is "both".
#' @param verbose Logical. Show or hide possible warnings and messages.
#'
#' @importFrom utils read.table setTxtProgressBar txtProgressBar
#' @importFrom lpSolve lp
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats complete.cases delete.response model.frame model.matrix na.fail na.omit terms update.formula
#'
#' @return An object of class "whatif".
#'
#' @seealso \code{\link{compute_nearby}}
#'
#' @references Gandrud C, King G, Stoll H, Zeng L (2017). WhatIf: Evaluate Counterfactuals. R package version 1.5-9. \href{https://CRAN.R-project.org/package=WhatIf}{https://CRAN.R-project.org/package=WhatIf}.
#'
#' King G, Zeng L (2007). When can history be our guide? The pitfalls of counterfactual inference. International Studies Quarterly 51, 183–210. DOI: \href{https://www.jstor.org/stable/pdf/4621707.pdf?seq=1#page_scan_tab_contents}{10.1111/j.1468-2478.2007.00445.x}
#'
#'@keywords internal

whatif <- function (formula = NULL, data, cfact, range = NULL, freq = NULL,
          nearby = 1, distance = "gower", miss = "list", choice = "both",
          return.inputs = FALSE, return.distance = FALSE, mc.cores = 2, verbose = TRUE,
          ...)
{

  if (mc.cores <= 0)
    stop("mc.cores must be an integer greater than 0.", call. = FALSE)
  if(verbose) message("Preprocessing data ...")

  # if (grepl("Zelig*", class(data)) & missing(cfact))
  #   cfact <- zelig_setx_to_df(data)
  # if (grepl("Zelig*", class(data)) & !missing(cfact)) {
  #   formula <- formula(delete.response(terms(data$formula)))
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
  if (!any(complete.cases(cfact))) {
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
      }
      else {
        data <- model.matrix(tt, model.frame(tt, data = eval(data$data,
                                                             envir = .GlobalEnv), na.action = NULL))
      }
    }
    else {
      data <- model.matrix(tt, data = data$model)
    }
    if (!(is.matrix(data))) {
      stop("observed covariate data could not be extracted from output object")
    }
    rm(tt)
  }
  else {
    if (!((is.character(data) && is.vector(data) && length(data) ==
           1) || is.data.frame(data) || (is.matrix(data) &&
                                         !is.character(data)))) {
      stop("'data' must be either a string, a R data frame, a R non-character matrix, or an output object")
    }
    if (is.character(data)) {
      data <- read.table(data)
    }
  }
  if (dim(data)[1] == 0) {
    stop("no observed covariate data supplied: 'data' contains zero rows")
  }
  if (!any(complete.cases(data))) {
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
    }
    else {
      stop("'formula' must be of class 'formula'")
    }
  }
  if (!(identical(complete.cases(cfact), rep(TRUE, dim(cfact)[1])))) {
    cfact <- na.omit(cfact)
    if(verbose) message("Note:  counterfactuals with missing values eliminated from cfact")
  }
  if (is.data.frame(data)) {
    if (is.character(as.matrix(data))) {
      stop("observed covariate data not coercable to numeric matrix due to character column(s)")
    }
    data <- suppressWarnings(data.matrix(data))
  }
  else {
    data <- data.matrix(as.data.frame(data))
  }
  if (is.data.frame(cfact)) {
    if (is.character(as.matrix(cfact))) {
      stop("counterfactual data not coercable to numeric matrix due to character column(s)")
    }
    cfact <- suppressWarnings(data.matrix(cfact))
  }
  else {
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
  if (!(is.null(range))) {
    if (!(is.vector(range) && is.numeric(range))) {
      stop("'range' must be a numeric vector")
    }
    if (!identical(length(range), ncol(data))) {
      stop("length of 'range' does not equal number of columns of 'data'")
    }
  }
  if (!(is.null(freq))) {
    if (!(is.vector(freq) && is.numeric(freq))) {
      stop("'freq' must be a numeric vector")
    }
    na.fail(freq)
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
  if (!(identical(distance, "gower") || identical(distance,
                                                  "euclidian"))) {
    stop("'distance' must be either ''gower'' or ''euclidian''")
  }
  if (!(identical(choice, "both") || identical(choice, "hull") ||
        identical(choice, "distance"))) {
    stop("'choice' must be either ''both'', ''hull'', or ''distance''")
  }
  if (!(is.logical(return.inputs))) {
    stop("'return.inputs' must be logical, i.e. either TRUE or FALSE")
  }
  if (!(is.logical(return.distance))) {
    stop("'return.distance' must be logical, i.e. either TRUE or FALSE")
  }
  n = nrow(data)
  convex.hull.test <- function(x, z, mc.cores = mc.cores) {
    one_core_pb <- mc.cores == 1
    n <- nrow(x)
    k <- ncol(x)
    m <- nrow(z)
    if (one_core_pb && m == 1)
      one_core_pb <- FALSE
    if (one_core_pb)
      pb <- txtProgressBar(min = 1, max = m, style = 3)
    A <- rbind(t(x), rep(1, n))
    C <- c(rep(0, n))
    D <- c(rep("=", k + 1))
    in_ch <- function(i, one_core_pb = FALSE) {
      B <- c(z[i, ], 1)
      lp.result <- lpSolve::lp(objective.in = C, const.mat = A,
                      const.dir = D, const.rhs = B)
      if (one_core_pb)
        setTxtProgressBar(pb, i)
      if (lp.result$status == 0)
        return(TRUE)
      else return(FALSE)
    }
    if (one_core_pb) {
      hull <- sapply(1:m, in_ch, one_core_pb = one_core_pb)
    }
    else {
      if (.Platform$OS.type == "windows")
        hull <- parallel::mclapply(1:m, in_ch, mc.cores = mc.cores)
      else hull <- pbmclapply(1:m, in_ch, mc.cores = mc.cores)
      hull <- unlist(hull)
    }
    if (one_core_pb)
      close(pb)
    return(hull)
  }
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
  calc.ed <- function(dat, cf) {
    n <- nrow(dat)
    m <- nrow(cf)
    dat <- t(dat)
    dist = matrix(0, m, n, dimnames = list(1:m, 1:n))
    for (i in 1:m) {
      temp <- (dat - cf[i, ])^2
      dist[i, ] <- (colSums(temp))
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
  calc.cumfreq <- function(freq, dist) {
    m <- length(freq)
    n <- ncol(dist)
    res <- matrix(0, n, m)
    for (i in 1:m) res[, i] <- (colSums(dist <= freq[i]))/nrow(dist)
    return(res)
  }
  if (identical(miss, "list")) {
    data <- na.omit(data)
    n <- nrow(data)
  }
  if ((choice == "both") | (choice == "hull")) {
    if(verbose) message("Performing convex hull test ...")
    test.result <- convex.hull.test(x = na.omit(data), z = cfact,
                                    mc.cores = mc.cores)
  }
  if ((choice == "both") | (choice == "distance")) {
    if(verbose) message("Calculating distances ....")
    if (identical(distance, "gower")) {
      samp.range <- apply(data, 2, max, na.rm = TRUE) -
        apply(data, 2, min, na.rm = TRUE)
      if (!is.null(range)) {
        w <- which(!is.na(range))
        samp.range[w] <- range[w]
      }
      if (identical(TRUE, any(samp.range == 0))) {
        if(verbose) message("Note:  range of at least one variable equals zero")
      }
      dist <- calc.gd(dat = data, cf = cfact, range = samp.range)
    }
    else {
      dist <- calc.ed(dat = na.omit(data), cf = cfact)
    }
    if(verbose) message("Calculating the geometric variance...")
    if (identical(distance, "gower")) {
      gv.x <- geom.var(dat = data, rang = samp.range)
    }
    else {
      gv.x <- 0.5 * mean(calc.ed(dat = na.omit(data), cf = na.omit(data)))
    }
    if (identical(miss, "case") && identical(distance, "euclidian")) {
      summary <- colSums(dist <= nearby * gv.x) * (1/nrow(na.omit(data)))
    }
    else {
      summary <- colSums(dist <= nearby * gv.x) * (1/n)
    }
    if(verbose) message("Calculating cumulative frequencies ...")
    if (is.null(freq)) {
      if (identical(distance, "gower")) {
        freqdist <- seq(0, 1, by = 0.05)
      }
      else {
        min.ed <- min(dist)
        max.ed <- max(dist)
        freqdist <- round(seq(min.ed, max.ed, by = (max.ed -
                                                      min.ed)/20), 2)
      }
    }
    else {
      freqdist <- freq
    }
    cumfreq <- calc.cumfreq(freq = freqdist, dist = dist)
    dimnames(cumfreq) <- list(seq(1, nrow(cfact), by = 1),
                              freqdist)
  }
  if(verbose) message("Finishing up ...")
  if (return.inputs) {
    if (choice == "both") {
      if (return.distance) {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), in.hull = test.result, dist = t(dist),
                    geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), in.hull = test.result, geom.var = gv.x,
                    sum.stat = summary, cum.freq = cumfreq)
      }
    }
    if (choice == "distance") {
      if (return.distance) {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), dist = t(dist), geom.var = gv.x,
                    sum.stat = summary, cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), geom.var = gv.x, sum.stat = summary,
                    cum.freq = cumfreq)
      }
    }
    if (choice == "hull") {
      out <- list(call = match.call(), inputs = list(data = data,
                                                     cfact = cfact), in.hull = test.result)
    }
  }
  else {
    if (choice == "both") {
      if (return.distance) {
        out <- list(call = match.call(), in.hull = test.result,
                    dist = t(dist), geom.var = gv.x, sum.stat = summary,
                    cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), in.hull = test.result,
                    geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
      }
    }
    if (choice == "distance") {
      if (return.distance) {
        out <- list(call = match.call(), dist = t(dist),
                    geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), geom.var = gv.x,
                    sum.stat = summary, cum.freq = cumfreq)
      }
    }
    if (choice == "hull") {
      out <- list(call = match.call(), in.hull = test.result)
    }
  }
  class(out) <- "whatif"
  return(invisible(out))
}
