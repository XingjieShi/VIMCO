#ifndef vb_hpp
#define vb_hpp

#include <stdio.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
 


List emInd (arma::fmat X, 
			arma::mat Y, 
			Rcpp::Nullable<Rcpp::NumericMatrix> mu0,
			Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0,
            Rcpp::Nullable<Rcpp::NumericVector> sigb0,
		 	Rcpp::Nullable<Rcpp::NumericVector> sige0,
            int maxit, 
			double epsStopLogLik);
 


List emMultiple (arma::fmat X,
                arma::mat Y,
            	Rcpp::Nullable<Rcpp::NumericMatrix> mu0,
               	Rcpp::Nullable<Rcpp::NumericVector> sigb0,
		 	    Rcpp::Nullable<Rcpp::NumericMatrix> Theta0,
			    Rcpp::Nullable<Rcpp::NumericVector> Lambda0,
		 	    Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0,
		 	    bool   fixLambda,
		        int    maxit, 
		        double epsStopLogLik);



			
#endif /* vb_hpp */
