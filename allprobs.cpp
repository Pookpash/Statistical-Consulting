#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <valarray>
using namespace Rcpp;
using namespace std;

//'Parameters for the pdfs (non self explanatory)
//'@params Numericvector x is the data vector
//'@params mu is the expectation E(X) based on the distribution of X
//'@params sigma is the Variance Var(X) based on the distribution of X
//'@params kappa is the ... for the von Mises dist
//'@params phi is the concentration parameter for the beta distribution. his definition is phi = alpha + beta, resp. shape1 + shape2

//'Parameters for the allprobs algorithm
//'@params nStates - number of States
//'@params nObs - number of Observations
//'@params data - data matrix of the 5 variables. For different variables needed, the code has to be changed.
//'@params mumat - matrix of mu's, mu's of one state for different variables in one column!
//'@params sigmat - matrix of secondary parameters for the dists. see mumat!

// [[Rcpp::export]]
//  density of gamma dist
arma::colvec dgamma_rcpp(NumericVector x, double mu, double sigma){
        
        arma::colvec res(x.size());
        
        // convert mean and sd to shape and scale
        double shape = pow(mu,2)/pow(sigma,2);
        double scale = pow(sigma,2)/mu;
        
        for(int i=0; i<x.size(); i++) {
                if(!arma::is_finite(x(i)))
                        res(i) = 1; // if missing observation
                else
                        res(i) = R::dgamma(x(i),shape,scale,0);
        }
        
        return res;
}

// [[Rcpp::export]]
// density of von Mises dist
arma::colvec dvm_rcpp(NumericVector x, double mu, double kappa){
        
        arma::colvec res(x.size());
        double b = R::bessel_i(kappa,0,2);
        
        for(int i=0; i<x.size(); i++) {
                if(!arma::is_finite(x(i)))
                        res(i) = 1; // is missing observation
                else
                        res(i) = 1/(2*M_PI*b)*pow((exp(cos(x(i)-mu)-1)),kappa);
        }
        
        return res;
}


// [[Rcpp::export]]
// compute joint probabilities
arma::mat allprobs_rcpp(int nStates, int nObs, arma::mat data, arma::mat mumat, arma::mat sigmat){
        
        arma::mat allProbs(nObs, nStates);
        allProbs.ones();
        
        NumericVector step(nObs);
        NumericVector angle(nObs);
        NumericVector surf(nObs);
        NumericVector dur(nObs);
        NumericVector dive(nObs);
        
        step = data.cols(0,0);
        angle = data.cols(1,1);
        surf = data.cols(2,2);
        dur = data.cols(3,3);
        dive = data.cols(4,4);
        
        for(int i=0; i<nStates; i++){
                
                arma::colvec stepProb(nObs);
                arma::colvec angleProb(nObs);
                arma::colvec surfProb(nObs);
                arma::colvec durProb(nObs);
                arma::colvec diveProb(nObs);
                
                
                
                //compute probs here
                stepProb = dgamma_rcpp(step, mumat(i,4), sigmat(i,4));
                angleProb = dvm_rcpp(angle, mumat(i,3), sigmat(i,3));
                surfProb = dgamma_rcpp(surf, mumat(i,2), sigmat(i,2));
                durProb = dgamma_rcpp(dur, mumat(i,1), sigmat(i,1));
                diveProb = dgamma_rcpp(dive, mumat(i,0), sigmat(i,0));
                
                for(int j=0; j<nObs; j++)
                {
                        float temp;
                        temp = stepProb[j]*angleProb[j]*surfProb[j]*durProb[j]*diveProb[j];
                        allProbs(j,i) = temp;
                }
                
        }
        
        return(allProbs);
}

// time comparisons (on Till's home desktop, so performance can vary depending on your setup)
//
// WIP! WIP! not done yet since no access to desktop!
//
// |   n   | states | time in R | time in C++ | times faster |
// |-------|--------|-----------|-------------|--------------|
// |   500 |        |     3.33s |       0.78s |          4.3 |
// |  2000 |        |    17.85s |       3.66s |          4.9 |
// |  5000 |        |    22.67s |       5.76s |          3.9 |
// |  5000 |        |   106.47s |      28.10s |          3.8 |
// | 20000 |        |    98.04s |      24.07s |          4.1 |


/*** R

*/
