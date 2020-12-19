#' @title bootstrapCI
#' @description compute the confidence interval of kernel density estimation through bootstrap method.
#' @param data the data from which the estimate is to be computed. For the default method a numeric vector
#' @param n integer:the number of replications
#' @param alpha significance level
#' @param border the color to draw the border of the confidence band,Use border = NA to omit borders
#' @param col the color for filling the confidence band,The default,NA
#' @return a matrix for storing the confidence interval for the density estimation \code{n}
#' @examples
#' \dontrun{
#' bootstrapCI(faithful$eruptions,10000,0.05,border = TRUE,col="blue")
#' }
#' @importFrom stats density
#' @importFrom graphics polygon lines
#' @import knitr
#' @import Rcpp
#' @import scales
#' @import stats
#' @import RANN
#' @import energy
#' @import Ball
#' @import boot
#' @import parallel
#' @import microbenchmark
#' @useDynLib StatComp20040
#' @export
bootstrapCI <- function(data,n,alpha,border = NULL, col = NA){
x <- data
est0 <- density(x)
set.seed(1)
  resm <- replicate(n, {
    x1 <- sample(x, replace=TRUE)
    density(x1, from=min(est0$x),to=max(est0$x))$y
  })
  CI <- apply(resm, 1, quantile, c(alpha/2, 1-alpha/2))
  plot(est0, ylim=range(CI), type='n')
  polygon(c(est0$x, rev(est0$x)),
          c(CI[1,], rev(CI[2,])),
          col=col, border=FALSE)
  lines(est0, lwd=2)
  return(as.matrix(CI))
} 