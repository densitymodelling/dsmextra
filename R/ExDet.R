#' EXtrapolation DETection tool
#'
#' Assesses univariate (Type I) and combinatorial (Type II) extrapolation between a reference system (ref) and a projection system (p). See Mesgaran et al. (2014) for an explanation. This function is an updated version of some original code from the \code{\link[ecospat]{ecospat}} package (Broennimann et al. 2016).
#'
#' @param ref Reference data. A data.frame with the values of the variables (i.e. columns) for each sample unit (e.g. segments).
#' @param tg Target data. A data.frame with the values of the variables (i.e. columns) for each point of the prediction extent.
#' @param xp Character string. Names of the covariates of interest.
#' @return Returns a tibble with four columns. (1) \code{ExDet}: Value of the ExDet metric (negative for univariate extrapolation, >1 for combinatorial extrapolation, within the range 0-1 for analogue conditions). (2) \code{mic_univariate}: Integer indicating the covariate with the largest contribution to univariate extrapolation (most influential covariate, MIC).  (3) \code{mic_combinatorial}: Integer indicating the covariate with the largest contribution to combinatorial extrapolation. (4) \code{mic}: most influential covariate.
#'
#' @author Phil J. Bouchet
#'
#' @references Mesgaran, M.B., Cousens, R.D. & Webber, B.L. (2014) \href{https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209}{Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models}. Diversity & Distributions, 20: 1147-1159, DOI: 10.1111/ddi.12209
#'
#' Broennimann O, Di Cola V, Guisan A (2016). ecospat: Spatial Ecology Miscellaneous Methods. R package version 2.1.1. \href{https://CRAN.R-project.org/package=ecospat}{https://CRAN.R-project.org/package=ecospat}.
#' @keywords internal

ExDet <- function(ref, tg, xp){

  #---------------------------------------------
  # Converts data to matrix form
  #---------------------------------------------

  tg <- as.matrix(tg)
  ref <- as.matrix(ref)

  #---------------------------------------------
  # Min/Max for each covariate in reference system
  #---------------------------------------------

  a <- apply(ref, 2, min, na.rm = TRUE)
  b <- apply(ref, 2, max, na.rm = TRUE)

  #---------------------------------------------
  # Matrices of min and max
  #---------------------------------------------

  minref <- matrix(a, nrow = nrow(tg), ncol = ncol(tg), byrow = TRUE)
  maxref <- matrix(b, nrow = nrow(tg), ncol = ncol(tg), byrow = TRUE)

  #---------------------------------------------
  # Use matrix algebra to calculate univariate extrapolation (NT1)
  #---------------------------------------------

  nt1.df <- data.frame(apply(array(data = c(tg - minref, maxref - tg,
                                            rep(0, nrow(tg) * ncol(tg))),
                                   dim = c(dim(tg), 3)),
                             c(1, 2), min)/(maxref - minref))
  names(nt1.df) <- xp

  #---------------------------------------------
  # Compute total NT1
  #---------------------------------------------

  nt1 <- rowSums(nt1.df)

  #---------------------------------------------
  # Identify most influential covariates (MIC) for NT1
  #---------------------------------------------

  mic_nt1 <- apply(nt1.df, 1, FUN = function(x) base::which.min(x))

  #---------------------------------------------
  # Set MIC(nt1) to NA within univariate range
  #---------------------------------------------

  # Zero indicates within univariate range
  univ.rge <- which(nt1 == 0)
  mic_nt1[univ.rge] <- NA

  #---------------------------------------------
  # Extract associated data from the target system
  #---------------------------------------------

  tg.univ <- matrix(tg[univ.rge, ], ncol = ncol(tg))

  aa <- apply(ref, 2, mean, na.rm = TRUE) # col means
  bb <- stats::var(ref, na.rm = TRUE) # covariance matrix

  #---------------------------------------------
  # Mahalanobis distance from centre of environmental space in reference system
  #---------------------------------------------

  mah.ref <- stats::mahalanobis(x = ref, center = aa, cov = bb)
  mah.pro <- stats::mahalanobis(x = tg.univ, center = aa, cov = bb)

  #---------------------------------------------
  # Maximum Mahalanobis distance
  #---------------------------------------------

  mah.max <- max(mah.ref[is.finite(mah.ref)])

  #---------------------------------------------
  # Combinatorial extrapolation (NT2) as % of that distance
  #---------------------------------------------

  nt2 <- mah.pro/mah.max

  #---------------------------------------------
  # Save values
  #---------------------------------------------

  nt1[univ.rge] <- nt2

  #---------------------------------------------
  # Identify MIC for NT2
  #---------------------------------------------

  # All combinations of covariates when one is dropped and the others remain

  if(length(xp) == 1) cov.combs <- matrix(1)
  if(length(xp) > 1) cov.combs <- utils::combn(x = 1:ncol(tg.univ), m = length(xp)-1)

  cov.combs <- as.list(data.frame(cov.combs))

  #---------------------------------------------
  # Means and variances
  #---------------------------------------------

  if(length(xp) == 1){

    cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[,.]), 2, mean))
    cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[,.])))

  } else {

    cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[,.]), 2, mean, na.rm = TRUE))
    cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[,.]), na.rm = TRUE))

  }

  #---------------------------------------------
  # Calculate Mahalanobis distance for each combination of covariates
  #---------------------------------------------

  # Need min of two analogue observations to calculate Mahalanobis distance
  if(nrow(tg.univ) < 2){

    warning("Only one prediction point within analogue conditions. Mahalanobis distances cannot be calculated.")
    mah_nt2 <- vector(mode = "list", length = length(cov.combs))

  } else {

    mah_nt2 <- purrr::pmap(.l = list(cov.combs, cov.aa, cov.bb),
                           .f = function(a, b, c) stats::mahalanobis(x = as.matrix(tg.univ[,a]),
                                                                     center = b,
                                                                     cov = c))
  }

  # if(length(xp) == 1){
  #
  #   mah_nt2 <- purrr::pmap(list(cov.combs, cov.aa, cov.bb),
  #                          function(a, b, c) stats::mahalanobis(x = as.matrix(tg.univ[,a]),
  #                                                               center = b,
  #                                                               cov = c))
  #
  # } else {
  #
  #   mah_nt2 <- purrr::pmap(list(cov.combs, cov.aa, cov.bb),
  #                          function(a, b, c) stats::mahalanobis(x = as.matrix(tg.univ[,a]),
  #                                                               center = b,
  #                                                               cov = c))
  # }


  #---------------------------------------------
  # Retrieve missing covariate names
  #---------------------------------------------

  # if(length(xp) == 1){

  mah_nt2 <- mah_nt2 %>% purrr::set_names(., xp)

  # }else{mah_nt2 <- mah_nt2 %>%
  #   purrr::set_names(.,cov.combs %>%
  #                      purrr::map(., ~colnames(tg)[which(!1:ncol(tg)%in%.)]))
  # }

  mah_nt2 <- mah_nt2 %>% purrr::map_df(., cbind)
  mah_nt2 <- as.matrix(mah_nt2)

  mic_nt2 <- 100 * (mah.pro - mah_nt2)/mah_nt2

  #---------------------------------------------
  # MIC (NT2)
  #---------------------------------------------

  mic_nt2 <- apply(mic_nt2, 1, FUN = function(x) base::which.max(x))

  #---------------------------------------------
  # Combine results
  #---------------------------------------------

  results <- tibble::tibble(ExDet = nt1, mic_univariate = mic_nt1, mic_combinatorial = NA)
  if(nrow(tg.univ) > 1)  results$mic_combinatorial[univ.rge] <- mic_nt2

  # Analog conditions have no MIC

  results <- results %>%
    dplyr::mutate(mic_combinatorial = ifelse(ExDet>=0 & ExDet<=1, NA, mic_combinatorial))

  # Combined MIC column

  results <- results %>%
    dplyr::mutate(mic = rowSums(.[2:3], na.rm = TRUE))

  return(results)

}
