#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double forwardalgo(arma::rowvec foo, arma::mat gamma, arma::mat allprobs, double lscale,int nrows) {
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

/*** R
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
result
*/
