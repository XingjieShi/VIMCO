// We can now use the BH package
// [[Rcpp::depends(BH, RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>


#include "plinkfun.hpp"
#include "vb.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
Rcpp::List vimco(std::string stringname, int nPheno, 
					bool   fit           = true, 
					bool   IndAsInit     = true,
	            	Rcpp::Nullable<Rcpp::NumericMatrix> mu0     = R_NilValue,
               		Rcpp::Nullable<Rcpp::NumericVector> sigb0   = R_NilValue,
		 	    	Rcpp::Nullable<Rcpp::NumericMatrix> Theta0  = R_NilValue,
			    	Rcpp::Nullable<Rcpp::NumericVector> Lambda0 = R_NilValue,
		 	   		Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0  = R_NilValue,
		 	    	bool   fixlambda     = true,
		        	int    maxit         = 10^3, 
		        	double epsStopLogLik = 10^-5){  
	// plink file prefix: stringname

	cout << "## Start loading fam file ... " << endl;
	// load phenotype (fam file)
	string famfile = stringname;
	famfile += ".fam";
	int N = getLineNum(famfile);

	IntegerVector sex(N);
	NumericMatrix pheno(N,nPheno);
	CharacterVector FID(N), IID(N);
	ReadPlinkFamFile3(famfile, FID, IID, pheno, N, nPheno);

	cout << "## Start loading genotype file ... " << endl;
	//read plink file (plink bim bed fam files)
	string bimfile = stringname;
	bimfile += ".bim";
	int P =  getLineNum(bimfile);
	unsigned* Xtmp = new unsigned[ N * P];

	//cout << " break 1: " << N1 << "X" << P1 << endl;
	readPlink(stringname, N, P, Xtmp);
	arma::Mat<unsigned> X(Xtmp, N, P, false, false);
	X.replace(3, 0);
	
	vector<string> snps = read_snpnames(bimfile, P);
	
	arma::mat Y = as<arma::mat>(pheno);
	int K = Y.n_cols;
	fmat Xf = conv_to<fmat>::from(X);
	X.clear();
	// add your algorithm here for X and Y.
	List vi, vi_Ind;
	if (fit) {
		if (IndAsInit) {
			vi_Ind = emInd(Xf, Y, 
				 				R_NilValue, R_NilValue,
				 				R_NilValue, R_NilValue,
				 				maxit, epsStopLogLik);

			Rcpp::NumericVector sigb = as<Rcpp::NumericVector>(vi_Ind["sigb"]);
			Rcpp::NumericMatrix mu = as<Rcpp::NumericMatrix>(vi_Ind["mu"]);
            Rcpp::NumericMatrix Theta(K, K);
			NumericVector Theta_diag = 1.0/as<Rcpp::NumericVector>(vi_Ind["sige"]);
			for (int i = 0; i < K; i++)
				Theta(i, i) = Theta_diag(i);

            Rcpp::NumericVector Lambda(P, 1.0);
            Rcpp::NumericMatrix Alpha = as<Rcpp::NumericMatrix>(vi_Ind["Alpha"]);

            vi = emMultiple(Xf, Y, mu, sigb, Theta, 
		                    Lambda, Alpha, fixlambda, 
		                    maxit, epsStopLogLik);

		}	else {
			vi = emMultiple(Xf, Y, mu0, sigb0, Theta0, 
		                    Lambda0, Alpha0, fixlambda, 
		                    maxit, epsStopLogLik);
		}
	} 
  
	// save the snp names.
	List output = List::create(Rcpp::Named("X")    = Xf,
						       Rcpp::Named("Y")    = Y,
							   Rcpp::Named("snps") = snps,
							   Rcpp::Named("vi_Ind")   = vi_Ind,
							   Rcpp::Named("vi")   = vi);

    Xf.clear();
	return output;
}
