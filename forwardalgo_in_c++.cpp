#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double forwardalgo(arma::rowvec x, arma::mat phi, arma::mat Gamma, arma::mat u, double l,
                   arma::rowvec shape, arma::rowvec scale,int nrows) {
        for(int t=2;t<nrows;t++) {
                u = phi*Gamma;
                arma::mat probvec = R::diag(dgamma(x[t],shape=shape,scale=scale,0));
                u = u * probvec;
                l = l+log(sum(u));
                phi = u/sum(u);
        }
        return l;
}


/*** R
x <- c(24.98643, 23.42482, 22.25363, 24.98658, 28.93724, 33.67242, 34.72468, 28.47550, 35.85791, 33.74849,
       35.77754, 12.02943, 14.94568, 21.50725, 21.50729, 41.63528, 42.27572, 30.74591, 20.49724, 43.89852)
phi <- matrix(c(7.341532e-05,0.9999266),nrow=1)
u <- matrix(c(7.607699e-07, 0.01036179),nrow=1)
Gamma <- matrix(c(0.7091270, 0.2908730, 0.2825269, 0.7174731),nrow=2,byrow=T)
l <- -4.569557
mu <- c(117.53181, 38.37499)
sig <- c(31.56159, 23.46791)
result <- forwardalgo(x,phi,Gamma,u,l,scale = mu^2/sig^2,shape = sig^2/mu,20))
result
*/
