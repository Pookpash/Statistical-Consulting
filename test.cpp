#include <Rcpp.h>
#include <Eigen/Dense>
using namespace Rcpp;

// [[Rcpp::export]]
double forwardalgo(arma::rowvec x, arma::mat phi, arma::mat Gamma, arma:mat u, double l, arma::rowvec mu, arma::rowvec sig,int nrows) {
        for(int t=2;t<nrows;t++)) {
                u = phi*Gamma;
                probvec = dgamma(x[t],shape=mu^2/sig^2,scale=sig^2/mu); //this is the last part that has to be rewritten I think (hope)
                u = u * probvec.asDiagonal();
                l = l+log(sum(u));
                phi = u/sum(u);
        }
        return l;
}


/*** R
forwardalgo(c(1,2,3,4,5,6,phi,Gamma,u,l,mu,sig,nrows))
*/
