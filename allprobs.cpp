#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <valarray>
using namespace Rcpp;
using namespace std;

//'Parameters for the pdfs (non self explanatory)
//'@params Numericvector x is the data vector

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
//density of beta dist
arma::colvec dbeta_rcpp(NumericVector x, double mu, double sigma){
        
        arma::colvec res(x.size());
        
        //convert mean and sd to shape1 and shape2
        double shape1 = (((1-mu)/pow(sigma,2))-(1/mu))*pow(mu,2);
        double shape2 = shape1*(1/mu-1);
        
        for(int i=0; i<x.size(); i++) {
                if(!arma::is_finite(x(i)))
                        res(i) = 1; // is missing observation
                else
                        res(i) = R::dbeta(x(i),shape1,shape2,0);
        }
        
        return res;
}

// [[Rcpp::export]]
// compute joint probabilities
arma::mat allprobs_rcpp(int nStates, int nObs, arma::mat data, arma::mat mumat, arma::mat sigmat){
        
        arma::mat allProbs(nObs, nStates);
        allProbs.ones();
        for(int i=0; i<nStates; i++){
                
                arma::colvec stepProb(nObs);
                arma::colvec angleProb(nObs);
                arma::colvec surfProb(nObs);
                arma::colvec durProb(nObs);
                arma::colvec diveProb(nObs);
                
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
               
                //compute probs here
                stepProb = dgamma_rcpp(step, mumat(i,0), sigmat(i,0));
                angleProb = dvm_rcpp(angle, mumat(i,1), sigmat(i,1));
                surfProb = dgamma_rcpp(surf, mumat(i,2), sigmat(i,2));
                durProb = dgamma_rcpp(dur, mumat(i,3), sigmat(i,3));
                diveProb = dbeta_rcpp(dive, mumat(i,4), sigmat(i,4));
                
                //cout << stepProb;
                //cout << angleProb;
                //cout << surfProb;
                //cout << durProb;
                //cout << diveProb;
                
                for(int j=0; j<nObs; j++)
                {
                        float temp;
                        temp = stepProb[j]*angleProb[j]*surfProb[j]*durProb[j]*diveProb[j];
                        //allProbs.col(i)[j] = temp;
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
gammatestdata <- c(0.041187346, 4.84, 0.039264483)
gammatestmu <- c(0.5)
gammatestsig <- c(2)
dgamma(gammatestdata, shape= gammatestmu^2/gammatestsig^2, scale = gammatestsig^2/gammatestmu)
dgamma_rcpp(gammatestdata, gammatestmu, gammatestsig)

betatestdata <- c(0.041187346, 0.8, 0.039264483)
betatestmu <- c(0.5)
betatestsig <- c(sqrt(0.05))
dbeta(betatestdata,shape1=2,shape2=2)
dbeta_rcpp(betatestdata, betatestmu, betatestsig)

N <- 2
nObs <- 7
data <- matrix(c(3,179,5,4,0.5,2,160,3,1,0.6,5,181,4,5,0.7,4,170,8,3,0.5,
                 9,182,4,2,0.4,8.8,181,2.2,2.5,0.43,6,186,2.1,0.9,0.45),nrow=7,byrow=T)

mumat <- matrix(c(4,180,5,2,0.5,6,165,3,4,0.6),ncol=5,byrow=T)
sigmat <- matrix(c(2,2,2.5,1,0.2,2.2,3,2.1,1.9,0.25),ncol=5,byrow=T)

allprobs_rcpp(N,nObs,data,mumat,sigmat)

allprobsR <- matrix(rep(1,N*nObs),nrow=nObs)
for (j in 1:2){
        surf.prob  <- dive.prob <- dep.prob <- step.prob <- rep(1,nObs)
        step.prob  <- dgamma(data[,1],  shape=mumat[j,1]^2/sigmat[j,1]^2, scale=sigmat[j,1]^2/mumat[j,1])
        angle.prob <- dvm_rcpp(data[,2],mumat[j,2],sigmat[j,2])
        surf.prob  <- dgamma(data[,3],  shape=mumat[j,3]^2/sigmat[j,3]^2, scale=sigmat[j,3]^2/mumat[j,3])
        dive.prob  <- dgamma(data[,4],  shape=mumat[j,4]^2/sigmat[j,4]^2, scale=sigmat[j,4]^2/mumat[j,4])
        shape1 <- ((1-mumat[j,5])/(sigmat[j,5]^2)-(1/mumat[j,5]))*(mumat[j,5]^2)
        dep.prob   <- dbeta(data[,5],   shape1=shape1, shape2=shape1*(1/mumat[j,5]-1),0)
        allprobsR[,j] <- surf.prob*angle.prob*dive.prob*dep.prob*step.prob
}
allprobsR

*/
