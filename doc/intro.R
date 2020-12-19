## ----eval=FALSE---------------------------------------------------------------
#  bootstrapCI <- function(data,n,alpha,border = NULL, col = NA){
#  x <- data
#  est0 <- density(x)
#  set.seed(1)
#    resm <- replicate(n, {
#      x1 <- sample(x, replace=TRUE)
#      density(x1, from=min(est0$x),to=max(est0$x))$y
#    })
#    CI <- apply(resm, 1, quantile, c(alpha/2, 1-alpha/2))  #confidence interval
#    plot(est0, ylim=range(CI), type='n')
#    polygon(c(est0$x, rev(est0$x)),
#            c(CI[1,], rev(CI[2,])),
#            col=col, border=FALSE)
#    lines(est0, lwd=2)        #estimated kernel density
#    return(as.matrix(CI))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  bootstrapCI(faithful$eruptions,10000,0.05,border = TRUE,col="blue")
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(FNN)
#  knnde <- function(data,k,dim,col=NA){
#    x <- as.matrix(data)
#    if(dim==2){
#      xrange <- seq(from = floor(min(x[,1])), to = ceiling(max(x[,1])), length.out=100)
#      yrange <- seq(from = floor(min(x[,2])), to = ceiling(max(x[,2])), length.out=100)
#      p <- ncol(x)
#      n <- nrow(x)
#      est_pt <- expand.grid(xrange, yrange)
#      distance <- knnx.dist(x, est_pt, k)
#      est_de <- matrix(k/(2*n*distance[,k]), nrow = length(xrange))
#      persp(xrange, yrange, est_de, phi = 30, theta = 45, col = col, border=0)
#    }
#  
#    if(dim==1){
#      xrange <- seq(from = floor(min(x)), to = ceiling(max(x)), length.out=50)
#      n <- nrow(x)
#      distance <- knnx.dist(x,xrange, k)
#      est_de <- matrix(k/(2*n*distance[,k]), nrow = length(xrange))
#      plot(xrange,est_de,type='l')
#    }
#    return(est_de)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  fit_knnde2 <- knnde(faithful, k=5, dim =2, col = 'red')
#  fit_knnde1 <- knnde(faithful$eruptions, k=5, dim =1)

## ----eval=FALSE---------------------------------------------------------------
#  RJtest <- function(x,alpha){
#    n <- length(x)
#    p <- numeric(n)
#    b <- numeric(n)
#    y <- sort(x)
#    cv <- numeric(1)
#    for (i in seq_along(x)) {
#      i <- rank(y)[i]    #the rank of the data
#      p[i] <- (i-3/8)/(n+1/4)
#      p <- p[i]
#      b[i] <- qnorm(p,0,1)  #the sorted data's normality scores
#    }
#    R <- sum(b*y)/sqrt(var(y)*(n-1)*sum(b^2))
#    if(alpha == 0.1){
#      cv <-1.0071-0.1371/sqrt(n)-0.3682/n+0.7780/n^2  #marginal value for different alpha
#    }
#    if(alpha == 0.05){
#      cv <-1.0063-0.1288/sqrt(n)-0.6118/n+1.3505/n^2
#    }
#    if(alpha == 0.01){
#      cv <-0.9963-0.0211/sqrt(n)-1.4106/n+3.1791/n^2
#    }
#    res <- c(R,cv,ifelse(R<cv,"nonnormal","normal"))
#    return(res)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  x <- rpois(100,2)
#  RJtest(x,0.01)
#  y <- rnorm(1000)
#  RJtest(y,0.01)

