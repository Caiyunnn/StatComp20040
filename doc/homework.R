## -----------------------------------------------------------------------------
library(car)
attach(sleep)                                                        #load data                                                                                     
data=sleep                                                
data$ID=as.numeric(sleep$ID)
scatterplot(extra~ID|group,data=data, main= "two drugs'effect",      #Create a scatter diagram
            legend = list(coords="topleft"),boxplot='y',             #set the parameters of graph
            xlab = "patientsID",regLine=FALSE,                       
            ylab = "increase in hours of sleep")
abline(h=0, col = "green",lwd=2,lty=1)                               #add line:Y=0
detach(sleep)

## -----------------------------------------------------------------------------
data(women)
attach(women)                                             #load data 
library(knitr)                                            #load package
lm.D1 = lm(height~weight,data = women)                    #regression
knitr::kable(anova(lm.D1),caption = "analysis of variance table")   #create a form
detach(women)

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)               
x <- 2/(1-u)^{1/2}  #the inverse transform method:F(x) = 1-(2/x)^2, x>=2
hist(x, prob = TRUE, main = expression(f(x)==frac(8,x^3))) #Graph the density histogram of the sample
y <- seq(2, 1000, .01)
lines(y, 8/y^3, col="red")

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1e5;k<-0;y <- numeric(n)
U1 <- runif(n,-1,1)
U2 <- runif(n,-1,1)
U3 <- runif(n,-1,1)                #U1,U2,U3~U(0,1)
while (k < n) {
  k <- k + 1
  if (abs(U1[k])<=abs(U3[k]) && abs(U2[k])<=abs(U3[k])) {
  U3[k] <- U2[k]
  }
} #when |U1|<=|U3|and |U2|<=|U3|,replace U3 with U2
hist(U3,  prob = TRUE,main = expression(f(x)==frac(3,4)*(1-x^2)),ylim = c(0,0.8))
# estimate of a large simulated random sample
y <- seq(-1,1, .01)
lines(y, (3/4)*(1-y^2),col="red")


## -----------------------------------------------------------------------------
set.seed(123)
n <- 1e5; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n, lambda) # the length of lambda = n
hist(x, prob = TRUE, main = expression(f(x)==frac(64,(2+y)^5)),ylim = c(0,0.5))
# the density histogram of the sample 
y <- seq(0, 1000, .01)
lines(y, 64/(2+y)^5, col="red")


## -----------------------------------------------------------------------------
library(knitr)
n <- 1e5
x <- runif(n, min=0, max=pi/3)         #generate random x~u(0,pi/3)
theta.hat <- mean(sin(x))*pi/3         #Monte Carlo estimate

 #creat a table to compare estimation with exact value
knitr::kable(data.frame(c(theta.hat,-cos(pi/3) + cos(0)),row.names=c('estimate','exact')),row.names = TRUE,col.names = "value", caption = "Comparison of estimation and exact value", align = "c")


## -----------------------------------------------------------------------------
library(knitr)
library(scales)
#generate random x~u(0,1)
n <- 1e5
x <- runif(n, min=0, max=1)

# simple Monte Carlo method
theta1.hat <- mean(exp(x))
theoreticalValue <- exp(1)-exp(0)

#antithetic variate approach
y <- x[1:(n/2)]
theta2.hat <- mean(exp(y)+exp(1-y))/2

#estimate two kinds methods' variation
var1 <- var(exp(x))/n
var2 <- var(exp(y)+exp(1-y))/(2*n)
reduction <- percent((var1-var2)/var1)

#create a table to campare
z<-data.frame(c(theta1.hat,theoreticalValue,var1),c(theta2.hat,theoreticalValue,var2),c("\\","\\",reduction),row.names = c("theta.hat","theta","variance"))
colnames(z)= c("simpleMC","antithetic","redutionPercent")
knitr::kable(z,row.names = TRUE, caption = "Comparison of two estimations",align = 'c')


## -----------------------------------------------------------------------------
library(knitr)
n <- 1e5
x <- runif(n)
y <- 1-log(1-x)     

#generate random numbers whose density is f2
u <- runif(n)
v <- 1/(1-u)
g <- function(x){(x^2 *exp(-x^2/2))/sqrt(2*pi)}
f1 <- function(x){exp(1-x)}  #set the first influence function
f2 <- function(x){1/x^2}  #set the second influence function
theta1 <- mean(g(y)/f1(y))
theta2 <- mean(g(v)/f2(v))
var1 <- var(g(y)/f1(y))/length(y)
var2 <- var(g(v)/f2(v))/length(u)
z <- seq(1,4,0.01)

#generate a graph to compare
plot(z, g(z)/f1(z), type = "l", ylab = "",ylim = c(0,3.2), lwd = 2, lty = 2,col=2,main='(ratio of g and f)')
lines(z,g(z)/f2(z), lty = 3, lwd = 2,col=3)
legend("topright",legend=c(expression(f[1](x)==e^{-x}),expression(f[2](x)==4/((1+x^2)*pi))),lty = 2:3, lwd = 2, inset = 0.02,col=2:3)


 #generate a table to compare
knitr::kable(data.frame(c(theta1,var1),c(theta2,var2),row.names=c('thetahat','thetahat.var')),row.names = TRUE,col.names = c("f1","f2"), caption = "Comparison of two estimations", align = "c")


## -----------------------------------------------------------------------------
n <- 1e4
r <- n/5         #replicates per stratum
N <- 50          #number of times to repeat the estimation
T2 <- numeric(5)
est <- matrix(0, N, 2)
g1<-function(x){(1/(1+x^2))*(x>0)*(x<1)}

for (j in 1:N){
  y <- rexp(n,rate = 1)
  est[j, 1] <- mean(g1(y))
for (i in 1:5){
    u <- runif(r)
    v <- -log(exp(-(i-1)/5)*(1-u)+exp(-i/5)*u)
    g <- function(x){((exp(-(i-1)/5)-exp(-i/5))/(1+x^2))*(x>(i-1)/5)*(x<i/5)}
    T2[i] <- mean(g(v))
}
  est[j, 2] <- sum(T2)
}
apply(est,2,mean)


## -----------------------------------------------------------------------------
alpha <- .05
UCL <- numeric(1000)
LCL <- numeric(1000)

#generate a function to compute the confidential interval
CI <- function(n,alpha){
 x <- rlnorm(n, meanlog = 0,sdlog = 1) 
 y <- log(x)
 UCL <-mean(y)+ sd(y)*qt(1-alpha/2, df = n-1)/sqrt(n)
 LCL <-mean(y)+ sd(y)*qt(alpha/2, df = n-1)/sqrt(n)
 return(c(LCL,UCL))
}

k <- replicate(1000,CI(20,0.05))
sum((k[1,]<0)*(k[2,]>0))/1000




## -----------------------------------------------------------------------------
alpha <- .05
UCL <- numeric(1000)
LCL <- numeric(1000)

#generate a function to compute the confidential interval
CI <- function(n,alpha){
 x <- rchisq(n,2) 
 UCL <-mean(x)+ sd(x)*qt(1-alpha/2, df = n-1)/sqrt(n)
 LCL <-mean(x)+ sd(x)*qt(alpha/2, df = n-1)/sqrt(n)
 return(c(LCL,UCL))
}

#x~chisquare distribution,df=2,E(x)=2
k <- replicate(1000,CI(20,0.05))
sum((k[1,]<2)*(k[2,]>2))/1000


## -----------------------------------------------------------------------------
#pure distribution
alpha <-0.05
#creat a function to compute the sample skewness statistic.
sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return(m3/m2^1.5)
}
n <- 40
cv <- qnorm(1-alpha/2, 0,sqrt(6*(n-2)/((n+1)*(n+3))))
m <- 4000
t <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(t)
test1 <- test2 <- numeric(m)
a <- b <- numeric(n)
POWER <- matrix(0,N,2)
for (i in 1:N) {
  p <- t[i]
  for (j in 1:m) {
    a <- sample(c(2, 10), replace = TRUE, size = n, prob = c(1-p,p))
    x <- rbeta(n,a,a)
    test1[j] <- as.integer(abs(sk(x)) >= cv)     #mixed beta distribution
    b <- sample(c(2, 5), replace = TRUE, size = n, prob = c(1-p,p))
    y <- rt(n,df=b)
    test2[j] <- as.integer(abs(sk(y)) >= cv)    #mixed t distribution
  }
  POWER[i,1] <- mean(test1)
  POWER[i,2] <- mean(test2)
}             
# graph to compare
plot(t, POWER[,1], xlab = bquote(t), ylim = c(0,1),col=2,type = "b")
lines(t,POWER[,2],xlab = bquote(t),ylim = c(0,1),col=3,type = "b")
legend("topright", legend=c("beta ditribution","t distribution"), col=c(2,3)) 


## -----------------------------------------------------------------------------
k <-matrix(0,2,3)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}

power <- function(n,alpha){
  x <- rnorm(n, 0, sigma1)
  y <- rnorm(n, 0, sigma2)
  pow1 <- count5test(x, y)
  pow2 <- as.integer(var.test(x,y,ratio = 1,conf.level = 1-alpha)$p.value < alpha)
  return(c(pow1,pow2))
}

# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 2
alpha  <- 0.055
m <- 1000
n <- c(20,100,1000)
for (i in 1:length(n)) {
  k[,i] <- apply(replicate(m,power(n[i],alpha)),MARGIN = 1,mean)
}
result <- as.data.frame(k,row.names = c("count five test","F test"))
colnames(result) <- c("small","medium","large")
print(result)

## -----------------------------------------------------------------------------
#matrix x is n*d dimension,(x1,x2....xn),every sample has d variables,and s is covariance matrix

b_1d <- function(x,d){
  xbar <- apply(x, margin=1, mean)
   y <- x-xbar
   s <- solve(cov(x))
   b_1d <- (t(y)%*%s%*%y)^3/n^2
}


## -----------------------------------------------------------------------------
maxout <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}
   

## -----------------------------------------------------------------------------
set.seed(1)
m<-1e3 #number of cycles
B<-1e3 #The number of bootstrap replicates
library(boot)#for boot and boot.ci function
#generates 4 different types of equi-tailed two-sided nonparametric confidence intervals
for(i in 1:m){
  boot.obj <- boot(data=aircondit,statistic=function(x,i){mean(x[i,])}, R=B )
  CI <- boot.ci(boot.obj,conf=0.95,type=c("norm","basic","perc","bca"))

}
print(CI)


## -----------------------------------------------------------------------------
vartest <- function(x, y) {
x0 <- x - mean(x)
y0 <- y - mean(y)
outx <- sum(x0 > max(y0)) + sum(x0 < min(y0))
outy <- sum(y0 > max(x0)) + sum(y0 < min(x0))
return(max(c(outx, outy)))
}

set.seed(1234)
R <- 999
n1 <- 20
n2 <- 40
n <- n1+n2  
m <- 1000
p <- numeric(R)
p.hat <- numeric(m)
test <- replicate(m,expr = {
x <- rnorm(n1,0,1)
y <- rnorm(n2,0,1)
z <- c(x, y)
temp <- vartest(x,y)
for (i in 1:R) {
ind <- sample(n, size = n1, replace = FALSE)
x1 <- z[ind]
y1 <- z[-ind] 
p[i] <-vartest(x1, y1)
}
p.hat <- (1+sum(as.integer(p>temp)))/1000
})
print(round(mean(test),6))

   

## -----------------------------------------------------------------------------
library(RANN)  #for NN test
library(energy) #for energy test
library(Ball)  #for ball test
library(boot)
Tn <- function(z, ix, sizes,k) {                  
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) 
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
  sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
m <- 1e3
p.values <- matrix(NA,m,3)




#Unequal variances and equal expectations
k<-3; p<-2; set.seed(12345)
n1 <- n2 <- 20; R<-999; n <- n1+n2; N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,1),ncol=p);
y <- matrix(rnorm(n2*p,0,2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
print(pow)




## -----------------------------------------------------------------------------
#Unequal variances and unequal expectations
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,1),ncol=p);
y <- matrix(rnorm(n2*p,1,2),ncol=p);
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
print(pow)


## -----------------------------------------------------------------------------
#Non-normal distributions
for(i in 1:m){
x <- matrix(rt(n1*p,df=1),ncol=p);
y <-cbind(rnorm(n2,0,1),rnorm(n2,1,2))
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value#performance of NN method
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#performance of energy method
p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value#performance of ball method
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
print(pow)


## -----------------------------------------------------------------------------
#Unbalanced samples
n1<-20;n2<-100;N=c(n1,n2)
for(i in 1:m){
  x <- matrix(rt(n1*p,1),ncol=p);
  y <- matrix(rnorm(n2*p,1,2),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#performance of NN method
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#performance of energy method
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed=i*12345)$p.value#performance of ball method
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
print(pow)


## -----------------------------------------------------------------------------
library(knitr)
set.seed(12345)
Ld <-function(x){exp(-abs(x))/2}
rw.Metropolis <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= (Ld(y)/ Ld(x[i-1])))
                    x[i] <- y  
                else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
    }

    N <- 2000
    sigma <- c(.05, .5,2,16)

    x0 <- 25
    rw1 <- rw.Metropolis(sigma[1], x0, N)
    rw2 <- rw.Metropolis(sigma[2], x0, N)
    rw3 <- rw.Metropolis(sigma[3], x0, N)
    rw4 <- rw.Metropolis(sigma[4], x0, N)
    

#number of candidate points rejected
no.reject <- data.frame(sigma=sigma,no.reject=c(rw1$k, rw2$k, rw3$k, rw4$k))
knitr::kable(no.reject,format='latex',caption = "Number of rejection per 2000")
matrix(c(no.reject[,1],1-no.reject[,2]/2000),ncol = 2,nrow=4,byrow=FALSE,dimnames=list(c(1, 2, 3, 4), c("sigma", "acceptance rates")))
   

   

## ----fig.width=7,fig.height=4-------------------------------------------------
par(mfrow=c(2,2))  
rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
for (j in 1:4) {
plot(rw[,j], type="l",
xlab=bquote(sigma == .(round(sigma[j],3))),
ylab="X", ylim=range(rw[,j]))
abline(h=-log(2))
abline(h=log(2))
}
par(mfrow=c(1,1)) #reset to default

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)

psi.means <- rowMeans(psi)     #row means
B <- n * var(psi.means)        #between variance est.
psi.w <- apply(psi, 1, "var")  #within variances
W <- mean(psi.w)               #within est.
v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
r.hat <- v.hat / W             #G-R statistic
return(r.hat)
}

Ld.chain <- function(sigma, N, X1) {
        #generates a Metropolis chain for Normal(0,1)
        #with Normal(X[t], sigma) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rnorm(1, xt, sigma)     #candidate point
            if (u[i] <= (Ld(y)/Ld(xt))) x[i] <- y
            else  x[i] <- xt
            }
        return(x)
}

sigma <- 2     #parameter of proposal distribution
k <- 4          #number of chains to generate
n <- 15000      #length of chains
b <- 500       #burn-in length

#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)

#generate the chains
set.seed(12345)
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k){
X[i, ] <- Ld.chain(sigma,n,x0[i])
}
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))

#plot psi for the four chains
# par(mfrow=c(2,2))for (i in 1:k)
 for (i in 1:k)
 if(i==1){
        plot((b+1):n,psi[i, (b+1):n],ylim=c(-0.2,0.2), type="l",
            xlab='Index', ylab=bquote(phi))
      }else{
        lines(psi[i, (b+1):n], col=i)
    }
par(mfrow=c(1,1)) #restore default


## ----echo=FALSE,fig.width=7,fig.height=4--------------------------------------
    
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)


## -----------------------------------------------------------------------------
b <- c(4:25,100,500,1000)
f <- function(a){
  s1 <- 1-pt((a^2*k/(k+1-a^2))^(1/2),df=k)
  s2 <- 1-pt((a^2*(k-1)/(k-a^2))^(1/2),df=k-1)
  s <- s1-s2
  return(s)
}
e <- 1e-6
root <- numeric(25)
f.root <- numeric(25)
h <- numeric(25)
x <- matrix(0,25,3)
for (i in 1:25) {
  h[i] <- b[i]
  k=h[i]
  root[i] <- uniroot(f,c(0+e,h[i]^(1/2)-e))$root
  f.root[i] < uniroot(f,c(0+e,h[i]^(1/2)-e))$f.root
}
x <- cbind(h,root,f.root)
knitr::kable(x,col.names=c("k","root","f.root"))


## -----------------------------------------------------------------------------
library("knitr")
nA <- 444
nB <- 132
nO <- 361
nAB <- 63
n <- sum(c(nA,nB,nO,nAB))
pqr <- matrix(0,100,3)
pql <- matrix(0,10,3)
ll <- numeric(100)

#estimate initial p0,q0,r0
l <- function(x){
  p <- x[1]
  q <- x[2]
  r <- 1-p-q
  -(nA*log(p^2 + 2*p*r)+nB*log(q^2 + 2*q*r) + 2*nO*log(r)+ nAB*log(2*p*q))
}
pq <- optim(c(0.3,0.1),l,method = "CG")$par
p0 <- pq[1]
q0 <- pq[2]
r0 <- 1-p0-q0
pqr[1,] <- c(p0,q0,r0)
ll[1] <- -l(pq)

#the conditional likelihood
El <- function(x,m){
  p <- x[1]
  q <- x[2]
  r <- 1-p-q
  base <- c(p^2,2*p*r,q^2,2*q*r,r^2,2*p*q)
  return(-sum(log(base)*m))
      
}

#EM estimatimate
i <- 1
while((i == 1)||(abs(ll[i]-ll[i-1]) >.Machine$double.eps^0.5)){
p <- pqr[i,1]
q <- pqr[i,2]
r <- 1-p-q
nAA <- nA * p/(p+2*r)
nAO <- nA - nAA
nBB <- nB * q/(q+2*r)
nBO <- nB - nBB
m <- c(nAA,nAO,nBB,nBO,nO,nAB)
res <- optim(c(p,q),El,method="CG",m=m) 
i <- i+1
pqr[i,] <- c(res$par,1-res$par[1]-res$par[2])
ll[i] <- -l(pqr[i,1:2])   #the corresponding log-maximum likelihood values (for observed data)
}

#create table and graph
pql <- cbind(pqr[1:i,1:2],ll[1:i])  #i is the iteration times
round(pql,6)
knitr::kable(pql,col.names = c("p","q","loglikValue"),align = "c",caption="EM estimates")
plot(pql[,3])
   

## -----------------------------------------------------------------------------
attach(mtcars)
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
lapply(formulas, lm)

for (i in 1:4) {
  lm(formulas[[i]])
}
detach(mtcars)


## -----------------------------------------------------------------------------
#Use sapply() and an anonymous function
pvalue1 <- numeric(100)
pvalue2 <- numeric(100)
results <- numeric(100)
i <- 1
trials <- replicate(100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
pvalue1 <- sapply(trials, function(x) x$p.value)

#Extra challenge
results <- replicate(100,expr = {
pvalue2 <- t.test(rpois(10, 10), rpois(7, 10))$p.value
},simplify = FALSE)
pvalue2 <- unlist(results)

## -----------------------------------------------------------------------------
set.seed(1)
library(parallel)

Mapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE){return(simplify2array(out))}
  out
}

#example
data1 <- matrix(rnorm(16, 0, 1), nrow = 4)
data2 <- matrix(rnorm(16, 0, 4), nrow = 4)
data <- cbind(data1,data2)
Mapply(data,mean,numeric(1))

## ----fig.width=7,fig.height=4-------------------------------------------------
library(knitr)
library(Rcpp)
library(StatComp20040)
N = 2000
sigma = c(.05, .5, 2, 16)
x0 = 25
rw1 = rwc(sigma[1],x0,N)
rw2 = rwc(sigma[2],x0,N)
rw3 = rwc(sigma[3],x0,N)
rw4 = rwc(sigma[4],x0,N)
#number of candidate points rejected
Rej = cbind(rw1[[2]], rw2[[2]], rw3[[2]], rw4[[2]])
Acc = round((N-Rej)/N,4)
rownames(Acc) = "Accept rates"
colnames(Acc) = paste("sigma",sigma)
knitr::kable(Acc)
#plot
par(mfrow=c(2,2))  #display 4 graphs together
rw = cbind(rw1[[1]], rw2[[1]], rw3[[1]], rw4[[1]])
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
    }


## ----fig.width=7,fig.height=4-------------------------------------------------
library(StatComp20040)
set.seed(1234)
N = 2000
sigma = c(.05, .5, 2, 16)
x0 = 25
rw1 = rwc(sigma[1],x0,N)
rw2 = rwc(sigma[2],x0,N)
rw3 = rwc(sigma[3],x0,N)
rw4 = rwc(sigma[4],x0,N)

rwr1 = rw.Metropolis(sigma[1],x0,N)
rwr2 = rw.Metropolis(sigma[2],x0,N)
rwr3 = rw.Metropolis(sigma[3],x0,N)
rwr4 = rw.Metropolis(sigma[4],x0,N)

rw = cbind(rw1[[1]], rw2[[1]], rw3[[1]], rw4[[1]])
rwr = cbind(rwr1$x, rwr2$x, rwr3$x,  rwr4$x)
par(mfrow=c(2,2))  #display 4 graphs together
for (j in 1:4) {
qqplot(rw[,j], rwr[,j],main=bquote(sigma == .(round(sigma[j],3))))
    }



## -----------------------------------------------------------------------------
library(microbenchmark)
library(StatComp20040)
N = 2000
sigma = c(.05, .5, 2, 16)
x0 = 25
ts <- microbenchmark(rwc1 = rwc(sigma[1],x0,N),rwr1 = rw.Metropolis(sigma[1],x0,N),times = 4)
summary(ts)[,c(1,3,5,6)]



