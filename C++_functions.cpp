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

//'Parameters for the forward algorithm w/o covariates
//'@params foo - 
//'@params gamma - t.p.m matrix
//'@params allprobs - matrix computed from the allprobs algorithm
//'@params lscale - scaled log likelihood
//'@params nrows - number of observations

//'Parameters for the create gamma cube when covariates are included
//'@params nbStates - number of States
//'@params beta - matrix with beta estimates
//'@params covs - matrix with covariates

//'Parameters for the forward algorithm with covariates
//'@params foo - 
//'@params gamma - t.p.m cube
//'@params allprobs - matrix computed from the allprobs algorithm
//'@params lscale - scaled log likelihood
//'@params nrows - number of observations

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


// [[Rcpp::export]]
double forwardalgo(arma::rowvec foo, arma::mat gamma, arma::mat allprobs, double lscale, int nrows) {
        for(int i=0;i<nrows;i++) {
                foo = foo*gamma%allprobs.row(i);
                double sumfoo = sum(foo);
                lscale = lscale+log(sumfoo);
                foo = foo/sumfoo;
        }
        return lscale;
}
// time comparisons (on Till's home desktop, so performance can vary depending on your setup)
//
// |   n   | states | time in R | time in C++ | times faster |
// |-------|--------|-----------|-------------|--------------|
// |   500 |      3 |     3.33s |       0.78s |          4.3 |
// |  2000 |      3 |    17.85s |       3.66s |          4.9 |
// |  5000 |      3 |    22.67s |       5.76s |          3.9 |
// |  5000 |      4 |   106.47s |      28.10s |          3.8 |
// | 20000 |      3 |    98.04s |      24.07s |          4.1 |


// [[Rcpp::export]]
// t.p.m. with covariates
arma::cube trMatrix_rcpp(int nbStates, arma::mat beta, arma::mat covs)
{
        int nbObs = covs.n_rows;
        arma::cube trMat(nbStates,nbStates,nbObs);
        trMat.zeros();
        arma::mat rowSums(nbStates,nbObs);
        rowSums.zeros();
        
        arma::mat g(nbObs,nbStates*(nbStates-1));
        g = covs*beta;
        
        for(int k=0;k<nbObs;k++) {
                int cpt=0;
                for(int i=0;i<nbStates;i++) {
                        for(int j=0;j<nbStates;j++) {
                                if(i==j) {
                                        trMat(i,j,k)=1;
                                        cpt++;
                                }
                                else trMat(i,j,k) = exp(g(k,i*nbStates+j-cpt));
                                rowSums(i,k)=rowSums(i,k)+trMat(i,j,k);
                        }
                }
        }
        
        // normalization
        for(int k=0;k<nbObs;k++)
                for(int i=0;i<nbStates;i++)
                        for(int j=0;j<nbStates;j++)
                                trMat(i,j,k) = trMat(i,j,k)/rowSums(i,k);
        
        return trMat;
}

// [[Rcpp::export]]
// forward algo with covariates
double forwardalgo_w_cov(arma::rowvec foo, arma::cube gamma, arma::mat allprobs, double lscale, int nrows) {
        for(int i=0;i<nrows;i++) {
                foo = foo*gamma.slice(i)%allprobs.row(i);
                double sumfoo = sum(foo);
                lscale = lscale+log(sumfoo);
                foo = foo/sumfoo;
        }
        return lscale;
}

/*** R
#check forwardalgo w/o covariates
foo <- c(0.8,0.2)
allprobs <- matrix(c(0.041187346, 4.840003e-04, 0.039264483, 3.077633e-04, 0.037015958, 2.121437e-04, 0.041187472,
                     4.840205e-04, 0.040742014, 1.264151e-03, 0.033085791, 3.036417e-03, 0.030864256, 3.570524e-03,
                     0.041144891, 1.143672e-03, 0.028404224, 4.202477e-03, 0.032928408, 3.073300e-03, 0.028579731, 
                     4.155762e-03, 0.005797453, 1.253671e-06, 0.013722984, 8.651436e-06, 0.035245663, 1.646933e-04,
                     0.035245781, 1.646960e-04, 0.016583294, 8.208760e-03, 0.015461981, 8.713352e-03, 0.038462354, 
                     1.822016e-03, 0.032477476, 1.144015e-04, 0.012844118, 1.002045e-02),nrow=20,byrow=T)
gamma <- matrix(c(0.8,0.2,0.8,0.2),nrow=2,byrow=T)
lscale <- 0
result <- forwardalgo(foo,gamma,allprobs,lscale,20)
result # should be -76.01029, then forwardalgo w/o covariates works

# check tr_matrix
N <- 2
beta <- matrix(c(1,2,1,3),byrow=T,nrow=2)
betat <- t(beta)
covs <- matrix(c(rep(1,10),1,2,1.2,4,7,4,3,2.3,3,9),byrow=F,ncol=2,nrow=10)
trmat <- trMatrix_rcpp(N,betat,covs)
trmat[,,3]
# should be a 2x2 Matrix with rowsums equal to 1 

#check forward algo with covariates

*/