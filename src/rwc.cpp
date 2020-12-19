#include <Rcpp.h>
using namespace Rcpp;

double lap_f(double x){
double z = exp(-abs(x));
return z;}
//' @title Implement a random walk Metropolis sampler.
//' @description Implement a random walk Metropolis sampler for generating the standard Laplace distribution
//' @param sigma the standard deviation of a normal distribution
//' @param x0 the initial value of a random walk sequence
//' @param N the length of a random walk sequence
//' @return  a random walk Metropolis sequence and number of candidate points rejected \code{n}
//' @examples
//' \dontrun{
//' set.seed(1234)
//' N = 2000
//' sigma = 1
//' x0 = 25
//' rw = rwc(sigma,x0,N)
//' }
//' @useDynLib StatComp20040
//' @export
// [[Rcpp::export]]
List rwc(double sigma,double x0, int N){
  NumericVector x(N);
  x[1] = x0;
  NumericVector u = runif(N);
  int k = 0;
  for (int i=1;i<N;i++) {
    double y = rnorm(1, x[i-1], sigma)[0];
    if(u[i] <= (lap_f(y)/lap_f(x[i-1]))){
      x[i] = y;}else{
        x[i] = x[i-1];
        k++;
      }
  }
  return(List::create(x = x, k = k));
}


