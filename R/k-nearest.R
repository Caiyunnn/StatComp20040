#' @title knnde
#' @description  the k-nearest neighbor density estimation 
#' @param data the data from which the estimate is to be computed. For the default method a numeric vector
#' @param k 	the maximum number of nearest neighbors to search. The default value is set to 10.
#' @param dim the dimension of data
#' @param col the color(s) of the surface facets.
#' @return the estimated kernel density \code{n}
#' @examples
#' \dontrun{
#' fit_knnde2 <- knnde(faithful, k=5, dim =2,col='red')
#' fit_knnde1 <- knnde(faithful$eruptions, k=5, dim =1)
#' }
#' @importFrom FNN knnx.dist
#' @importFrom graphics persp
#' @export
knnde <- function(data,k,dim,col=NA){
  x <- as.matrix(data)
  if(dim==2){
    xrange <- seq(from = floor(min(x[,1])), to = ceiling(max(x[,1])), length.out=100)
    yrange <- seq(from = floor(min(x[,2])), to = ceiling(max(x[,2])), length.out=100)
    p <- ncol(x)
    n <- nrow(x)
    est_pt <- expand.grid(xrange, yrange)
    distance <- knnx.dist(x, est_pt, k)
    est_de <- matrix(k/(2*n*distance[,k]), nrow = length(xrange))
    persp(xrange, yrange, est_de, phi = 30, theta = 45, col = col, border=0)
  }
  
  if(dim==1){
    xrange <- seq(from = floor(min(x)), to = ceiling(max(x)), length.out=50)
    n <- nrow(x)
    distance <- knnx.dist(x,xrange, k)
    est_de <- matrix(k/(2*n*distance[,k]), nrow = length(xrange))
    plot(xrange,est_de,type='l')
  }
  return(est_de)
}



