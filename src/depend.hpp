
#ifndef depend_hpp
#define depend_hpp

 
#include <RcppArmadillo.h>
#include <stdio.h>

using namespace Rcpp;
using namespace arma;



double logpexp(double x);

vec  xtx_diag(char* Xtmp, int n, int p, arma::mat xmean, arma::vec xsd);

vec VdotC (char* x, double c, int n, double xmean_j, double xsd_j);


double VdotProd (char* x, arma::mat y, int n, double xmean_j, double xsd_j);


mat MdotProd(char* Xtmp, arma::mat Y, int p, arma::mat xmean, arma::vec xsd);


mat XdotB(char* Xtmp, arma::mat B, int n, arma::mat xmean, arma::vec xsd);

mat Vouter(arma::rowvec x, char* y, int n, double mean_j, double sd_j);

void Centering(char* Xtmp, double* meanOut, double* sdOut, int n, int p);

double LogLikInd(arma::vec& xx, arma::mat& E, 
           int p,      int    K,
           arma::mat& mu,     arma::mat&    s, 
           arma::vec& sige,   arma::vec&    sigb, 
		   arma::mat& Alpha,  arma::vec& alpha);

double uThree(arma::mat& Theta, arma::rowvec mu_j, arma::rowvec Alpha_j, double XtX_jj);

double mu_jkThree(arma::rowvec theta_k, arma::rowvec Alpha_j, arma::rowvec mu_j, double XtX_jj, int k);

double LogLikMultiple(arma::vec& xx, arma::mat& E, 
           int p,      int    K,
           arma::mat& mu,     arma::mat&    s, 
           arma::mat& Theta,  arma::vec&    sigb, 
		   arma::vec& Lambda, double lambda,
		   mat& Alpha,        arma::vec alpha);

#endif /* depend_hpp */
