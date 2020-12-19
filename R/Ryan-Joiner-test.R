#' @title Ryan-Joiner-test
#' @description  Ryan-Joiner test for normality test
#' @param x the data from which the estimate is to be computed. For the default method a numeric vector
#' @param alpha 	significance level
#' @return Ryan-Joiner test statistic,marginal value,the result of normality test\code{n}
#' @examples
#' \dontrun{
#' x <- rpois(100,2)
#' RJtest(x,0.01)
#' y <- rnorm(1000)
#' RJtest(y,0.01)
#' }
#' @importFrom stats qnorm rnorm rpois
#' @useDynLib StatComp20040 
#' @export
RJtest <- function(x,alpha){
  n <- length(x)
  p <- numeric(n)
  b <- numeric(n)
  y <- sort(x)
  cv <- numeric(1)
  for (i in seq_along(x)) {
    i <- rank(y)[i]    #the rank of the data
    p[i] <- (i-3/8)/(n+1/4)  
    p <- p[i]
    b[i] <- qnorm(p,0,1)  #the sorted data's normality scores
  }
  R <- sum(b*y)/sqrt(var(y)*(n-1)*sum(b^2))
  if(alpha == 0.1){
    cv <-1.0071-0.1371/sqrt(n)-0.3682/n+0.7780/n^2  #marginal value for different alpha
  }
  if(alpha == 0.05){
    cv <-1.0063-0.1288/sqrt(n)-0.6118/n+1.3505/n^2
  } 
  if(alpha == 0.01){
    cv <-0.9963-0.0211/sqrt(n)-1.4106/n+3.1791/n^2
  } 
  res <- c(R,cv,ifelse(R<cv,"nonnormal","normal"))
  return(res)
}


