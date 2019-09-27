#'---------------------------------------------
# Function to tally the number/% of cells subject to extrapolation
#'---------------------------------------------

n_and_p <- function(x){

  exl <- list(univariate.n = length(data[x < 0]),
              univariate.p = 100 * length(x[x < 0])/length(x),
              combinatorial.n = length(x[x > 1]),
              combinatorial.p = 100 * length(x[x > 1])/length(x),
              analogue.n = length(x[x >= 0 & x <=1]),
              analogue.p = 100 * length(x[x >= 0 & x <=1])/length(x))
  return(exl)
}

#'---------------------------------------------
# Take mean of raster values - used in rasterize
#'---------------------------------------------

mean_ras <- function(x, ...){mean(x, na.rm = TRUE)}
