#' Summarize extrapolation results
#'
#' Print the summary of extrapolation results (already calculated by \code{\link{compute_extrapolation}}).
#'
#' @export
#' @author David L Miller
#' @param object result of running \code{\link{compute_extrapolation}}
#' @param digits precision of results
#' @param \dots for S3 compatability
#' @return invisibly returns the summary part of the object only, printing the results
#' @keywords internal

summary.extrapolation_results <- function(object, digits, ...){

  object$summary

}
