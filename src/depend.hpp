
#ifndef depend_hpp
#define depend_hpp

 
#include <RcppArmadillo.h>
#include <stdio.h>

using namespace Rcpp;
using namespace arma;



double logpexp(double x);

fvec  xtx_diag(char* Xtmp, int n, int p, arma::fmat xmean, arma::fvec xsd);

vec VdotC (char* x, double c, int n, float xmean_j, float xsd_j);


float VdotProd (char* x, arma::mat y, int n, float xmean_j, float xsd_j);


fmat MdotProd(char* Xtmp, arma::mat Y, int p, arma::fmat xmean, arma::fvec xsd);


mat XdotB(char* Xtmp, arma::mat B, int n, arma::fmat xmean, arma::fvec xsd);

mat Vouter(arma::rowvec x, char* y, int n, float mean_j, float sd_j);

void Centering(char* Xtmp, float* meanOut, float* sdOut, int n, int p);

double LogLikInd(arma::fvec& xx, arma::mat& E, 
           int p,      int    K,
           arma::mat& mu,     arma::mat&    s, 
           arma::vec& sige,   arma::vec&    sigb, 
		   arma::mat& Alpha,  arma::vec& alpha);

double uThree(arma::mat& Theta, arma::rowvec mu_j, arma::rowvec Alpha_j, float XtX_jj);

double mu_jkThree(arma::rowvec theta_k, arma::rowvec Alpha_j, arma::rowvec mu_j, float XtX_jj, int k);

double LogLikMultiple(arma::fvec& xx, arma::mat& E, 
           int p,      int    K,
           arma::mat& mu,     arma::mat&    s, 
           arma::mat& Theta,  arma::vec&    sigb, 
		   arma::vec& Lambda, double lambda,
		   mat& Alpha,        arma::vec alpha);

#endif /* depend_hpp */
