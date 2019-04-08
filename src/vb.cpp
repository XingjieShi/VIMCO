
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <time.h>

#include "depend.hpp"
#include "plinkfun.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export()]]
List emInd (arma::fmat X,
            arma::mat Y,
            Rcpp::Nullable<Rcpp::NumericMatrix> mu0    = R_NilValue, 
            Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0 = R_NilValue,
            Rcpp::Nullable<Rcpp::NumericVector> sigb0  = R_NilValue,
            Rcpp::Nullable<Rcpp::NumericVector> sige0  = R_NilValue,
                  int maxit = 10^3,
		       double epsStopLogLik = 10^-5
        ) {
    // inputs
    const int n = Y.n_rows ;
	const int K = Y.n_cols ;
    const int p = X.n_cols ;

	int i=0, j=0, k=0;
	
    //centeralize X and Y
    fmat xmean = mean(X);
    for (i = 0; i < n; i++) {
  		for (j = 0; j < p; j++) {
  			X(i,j) -= xmean(j);
  		}	
    } 
    mat ymean = mean(Y);
    for (i = 0; i < n; i++) {
  		for (k = 0; k < K; k++) {
  			Y(i,k) -= ymean(k);
  		}	
    } 
  
	
    // containers
	mat E (n, K, fill::zeros);
    mat mu(p, K, fill::zeros) ;
    if (mu0.isNotNull())  {
		mu = as<mat>(mu0);
    }

    mat s (p, K, fill::ones) ;

    vec sigb(K, fill::ones);
	if (sigb0.isNotNull()) {
     	sigb = as<vec>(sigb0);
    }  else sigb = vectorise(var(Y, 0, 0)/2);

    vec sige(K, fill::ones);
    if (sige0.isNotNull())	{
    	sige = as<vec>(sige0);
    } else sige = sigb;

    mat Alpha(p, K);	
    if (Alpha0.isNotNull()) {
		Alpha = as<mat>(Alpha0);
    } else Alpha.fill(0.01);

	vec alpha(K);
	if (Alpha0.isNotNull()) {
		alpha = vectorise(mean(Alpha, 0));
    }  else alpha.fill(0.01);
   
    mat Beta(p, K, fill::zeros);
	if (mu0.isNotNull()) {
    	for (k=0; k < K; k++) {
			Beta.col(k) = Alpha.col(k) % mu.col(k);
		}
	}	

    // precomputation
    fvec xx(p); 
    for (j=0; j < p; j++) {
		 xx(j) = sum(X.col(j) % X.col(j));
	}   	

	fmat XY = X.t() * conv_to<fmat>::from(Y);
	mat ytilde (n, K);
	vec ytildekj(n, fill::zeros);

	int it = 0;
	double L0 = -INFINITY;
	vec ell(maxit, fill::zeros);
    // algorithm
    while (it < maxit) {
		
		
         // E step 	
		for (k = 0; k < K; k++) {
			for (j = 0; j < p; j++) { 
			s (j, k) = sige(k)/(xx(j) + sige(k)/sigb(k)); 
			}
		}	
		
		for (k = 0; k < K; k++) {
			ytilde.col(k) = X * Beta.col(k); 
			for (j = 0; j < p; j++) {
				ytildekj = ytilde.col(k) - Beta(j, k) * X.col(j);
				
				double temp = sum(X.col(j) % ytildekj);
				mu(j, k) = (XY(j, k) - temp)/(xx(j) + sige(k)/sigb(k));
			
				// update Alpha
			    double ujk = log(alpha(k)/(1-alpha(k))) + 
				             0.5*(pow(mu(j,k), 2.0)/s(j,k) + log(s(j,k)/sigb(k)));				 
			    Alpha(j, k) = 1/(1+exp(-ujk));	
			    
				Beta(j, k)   = Alpha(j,k) * mu(j,k);
				// update ytilde
				ytilde.col(k) = ytildekj +  Beta(j, k) * X.col(j);//X * (Alpha % mu.col(k));
		    }	
		}		

	    
		 // M step
		for (k = 0; k < K; k++) {
			// update sigma_{e_k}
		    E.col(k) = Y.col(k) - X * Beta.col(k); 
            vec temp = (Alpha.col(k) % (pow(mu.col(k), 2.0) + s.col(k)) - pow(Alpha.col(k) % mu.col(k), 2.0)) % xx;
			sige(k) = (sum(E.col(k) % E.col(k)) + sum(temp))/n;
			//if (sige(k) < 10E-6) sige(k) = 10E-6;
			// update sigma_{beta_k}
			sigb(k) = sum(Alpha.col(k) % (pow(mu.col(k), 2.0) + s.col(k)))/sum(Alpha.col(k));
			//if (sigb(k) < 10E-6) sigb(k) = 10E-6;
		    alpha(k) = sum(Alpha.col(k))/p; 
		}    
		
		ell(it) = LogLikInd(xx, E, p, K, mu, s, sige, sigb, Alpha, alpha);

		// convergence 
		if (ell(it) < L0){
            printf("Lowerbound decreasing, Error at iteration %d th iteration, diff=%g", it, ell(it) - L0);
            break;
        } else if ((ell(it) - L0) < epsStopLogLik) {
			printf("Converge at %d th iteration, LowerBound = %f \n", it, ell(it));
            break;
		}
		L0 = ell(it);

		it++;
    }
	if (it >= maxit) it = it - 1;
    vec LowerBound(it, fill::zeros);
    LowerBound = ell.subvec(0, it);   
    
    // intercepts
    mat Beta0(ymean - xmean*Beta);
    // returns
    List ret ;
    ret["N"] = n ;
    ret["mu"] = mu ;
	ret["s"]  = s;
	ret["Beta0"] = Beta0;
	ret["Beta"]  = Beta;
	ret["Alpha"] = Alpha ;
	ret["fdr"]   = 1 - Alpha;
    ret["sigb"]  = sigb;
	ret["sige"]  = sige;
	ret["alpha"] = alpha;
	ret["iter"]  = it+1;
	ret["eps"]   = ell(it) - ell(it-1); 
	ret["LowerBound"] = LowerBound;

    return(ret) ;
}




// [[Rcpp::export()]]
List emMultiple (arma::fmat X,
                arma::mat Y,
            	Rcpp::Nullable<Rcpp::NumericMatrix> mu0     = R_NilValue,
               	Rcpp::Nullable<Rcpp::NumericVector> sigb0   = R_NilValue,
		 	    Rcpp::Nullable<Rcpp::NumericMatrix> Theta0  = R_NilValue,
			    Rcpp::Nullable<Rcpp::NumericVector> Lambda0 = R_NilValue,
		 	    Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0  = R_NilValue,
		 	    bool   fixlambda     = true, //if Lambda0 is not NULL, fixlambda should be false
		        int    maxit         = 10^3, 
		        double epsStopLogLik = 10^-5
        ) {
    // inputs
    const int N = Y.n_rows ;
	const int K = Y.n_cols ;
    const int p = X.n_cols ;

	int i=0, j=0, k=0;
	
    //centeralize X and Y
    fmat xmean = mean(X);
    for (i = 0; i < N; i++) {
  		for (j = 0; j < p; j++) {
  			X(i,j) -= xmean(j);
  		}	
    } 
    mat ymean = mean(Y);
    for (i = 0; i < N; i++) {
  		for (k = 0; k < K; k++) {
  			Y(i,k) -= ymean(k);
  		}	
    } 
  

    // containers
	mat E (N, K, fill::zeros);

	mat mu(p, K, fill::zeros);
	if (mu0.isNotNull())  {
		mu = as<mat>(mu0);
    }

    rowvec mu_jNumer(K, fill::zeros);


    mat s (p, K, fill::ones) ;

    vec sigb(K, fill::ones);
    sigb = vectorise(var(Y, 0, 0)/2.0);
    if (sigb0.isNotNull()) {
     	sigb = as<vec>(sigb0);
    }
    
    mat Theta(K, K, fill::zeros);
	Theta.diag() = vectorise(2.0/var(Y, 0, 0));
    if (Theta0.isNotNull()) {
    	Theta = as<mat>(Theta0);
    } 
	
	mat Sigma(K, K, fill::zeros);

    // Alpha
    mat Alpha(p, K, fill::zeros);	
    Alpha.fill(0.2);
    if (Alpha0.isNotNull()) {
		Alpha = as<mat>(Alpha0);
    }
	
	vec alpha =  vectorise(mean(Alpha, 0));
	
	//Lambda
	vec Lambda(p, fill::ones);	
	bool fixLambda = false;
	double lambda = 1;
	if (Lambda0.isNotNull()) {
		Lambda = as<vec>(Lambda0);	
		lambda = mean(Lambda);
	} 
	if (fixlambda & (lambda == 1.0))
		fixLambda = true;	

    mat Beta(p, K, fill::zeros);
    // If mu_jk has initial values, then calculate. 
    if (mu0.isNotNull()) {
    	for (k=0; k < K; k++) {
			Beta.col(k) = Lambda % Alpha.col(k) % mu.col(k);
		}
	}	
    // precomputation
	fvec xx(p); 
    for (j=0; j < p; j++) {
		 xx(j) = sum(X.col(j) % X.col(j));
	}   	
	
	// initialization
	int it = 0;
	double L0 = -INFINITY;
	vec X_j(N, fill::zeros);
	vec ell(maxit, fill::zeros);
	mat ytilde = X * Beta;
	mat y_j(N, K, fill::zeros);
	mat y_j_new(N, K, fill::zeros);

    // algorithm
    while (it < maxit) {
         
		// E step
		double logitlambda = log(lambda/(1.0-lambda));
		vec logitalpha = log(alpha/(1.0-alpha));
			// update mu, s_jk and Alpha.
            for (j = 0; j < p; j++) {   
                X_j = conv_to<vec>::from(X.col(j));
                //y_j = ytilde - Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
                for (k = 0; k < K; k++) {
                	y_j = ytilde - Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);

                    // update s_jk
			        s (j, k) = 1.0/(xx(j)*Theta(k,k) + 1.0/sigb(k));
                    
                    // update mu_jk
                    vec residual_j = (Y - y_j) * Theta.col(k);
                    double mu12 = sum(X.col(j) % residual_j);
                    double mu3 = mu_jkThree(Theta.row(k), Alpha.row(j), mu.row(j), xx(j), k);			   	
                    mu(j, k) = (mu12 - mu3) * s(j, k);

                    // update Alpha_jk
                    double v_jk = logitalpha(k) + Lambda(j)/2.0 * (mu(j,k)*mu(j,k)/s(j,k) + log(s(j,k)/sigb(k)));
                    Alpha(j, k) = 1.0 /(1.0+exp(-v_jk));
                    ytilde = y_j + Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
				}
				
				// update Lambda_j
                if (!fixLambda) {
                	y_j_new = ytilde - Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
                	for (k = 0; k < K; k++)	{
                        vec residual_j = (Y - y_j_new) * Theta.col(k);
                    	double mu12 = sum(X.col(j) % residual_j);
                    	double mu3 = mu_jkThree(Theta.row(k), Alpha.row(j), mu.row(j), xx(j), k);			   	
                    	mu_jNumer(k) = (mu12 - mu3);
                    	//ytilde = y_j + Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
                    }
                    double u_j2 = sum(Alpha.row(j) % mu.row(j) % mu_jNumer) + 
                             sum(Alpha.row(j) % (-mu.row(j)%mu.row(j)/s.row(j) + log(s.row(j)/sigb.t())))/2.0;
                    double u_j3 = uThree(Theta, mu.row(j), Alpha.row(j), xx(j));
                	double u_j  = logitlambda + u_j2 + u_j3;
                	Lambda(j) = 1.0/(1.0+exp(-u_j));
          	    }
          	    ytilde = y_j + Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
			}
	    
        // M step
		for (k=0; k < K; k++) {
		 	E.col(k) = Y.col(k) - ytilde.col(k);
		}
		// update Theta^-1        		
		 for (k=0; k < K; k++) {
		 	vec Siama_kk2 = Lambda%Alpha.col(k)%(mu.col(k)%mu.col(k)+s.col(k)) - 
			                  Lambda%Lambda%Alpha.col(k)%Alpha.col(k)%mu.col(k)%mu.col(k);
			Sigma(k, k) = (sum(E.col(k) % E.col(k)) + sum(xx % Siama_kk2))/N;
			for (int t=k+1; t < K; t++) {
				vec Siama_kt2 = (Lambda - Lambda%Lambda) % Alpha.col(k) % mu.col(k) % Alpha.col(t) % mu.col(t);
			    Sigma(k, t) = (sum(E.col(k) % E.col(t)) + sum(xx % Siama_kt2))/N;
                Sigma(t, k) = Sigma(k, t);				
			}
		}
	    Theta = inv_sympd(Sigma);
       
		// update sigma_{beta_k}
		for (k = 0; k < K; k++) {
			vec la = Lambda % Alpha.col(k);
			vec mu_k2 = mu.col(k) % mu.col(k);
			sigb(k) = sum(la % (mu_k2 + s.col(k))) / sum(la) ;
			alpha(k) = mean(Alpha.col(k));
		} 
        if(!fixlambda) lambda = mean(Lambda);		
		
		
		ell(it) = LogLikMultiple(xx, E, p, K, mu, s, Theta, sigb, Lambda, lambda, Alpha, alpha);
		// convergence 
		if (ell(it) < L0){
            printf("Lowerbound decreasing, Error at iteration %d th iteration, diff=%g \n", it, ell(it) - L0);
            break;
        } else if ((ell(it) - L0) < epsStopLogLik) {
			printf("Converge at %d th iteration, LowerBound = %f \n", it, ell(it));
            break;
		}
		L0 = ell(it);
		
		it++;
    } 
	
		
	
    if (it >= maxit) it = it - 1;
    vec LowerBound(it, fill::zeros);
    LowerBound = ell.subvec(0, it);     
	
	// local FDR
    mat fdr(p, K, fill::zeros);
    fdr = 1 - Alpha.each_col() % Lambda;
	
	//for (k=0; k < K; k++) {
		Beta = (Lambda % Alpha.each_col()) % mu;
	//}
    // intercepts
	mat Beta0(ymean - xmean*Beta);
    // returns
    List ret ;
    ret["N"]      = N;
    ret["mu"]     = mu;
	ret["s"]      = s;
	ret["Beta0"]  = Beta0;
	ret["Beta"]   = Beta;
	ret["Alpha"]  = Alpha;
	ret["alpha"]  = alpha;
	ret["Lambda"] = Lambda;
	ret["lambda"] = lambda;
	ret["fdr"]    = fdr;
    ret["sigb"]   = sigb;
	ret["Sigma"]  = Sigma;
	ret["Theta"]  = Theta; 
	ret["iter"]   = it+1;
	ret["eps"]    = ell(it) - ell(it-1); 
	ret["LowerBound"] = LowerBound;
    
    return(ret) ;
}




List emIndPlink (char* Xtmp, unsigned long long p,
            arma::mat Y,
            Rcpp::Nullable<Rcpp::NumericMatrix> mu0    = R_NilValue, 
            Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0 = R_NilValue,
            Rcpp::Nullable<Rcpp::NumericVector> sigb0  = R_NilValue,
            Rcpp::Nullable<Rcpp::NumericVector> sige0  = R_NilValue,
                  int maxit = 10^3,
		       double epsStopLogLik = 10^-5
        ) {
    // inputs
    const int n = Y.n_rows ;
	const int K = Y.n_cols ;

	int i=0, j=0, k=0;
	
	//fvec tmp = vectorise(X);
	//float* Xtmp = &tmp[0];
   
    float* mean_x = new float[p];
    float* sd_x = new float[p];
    Centering(Xtmp, mean_x, sd_x, n, p);
    fmat xmean(mean_x, 1, p, false, true); 
    fvec xsd(sd_x, p, false, true);
    //centeralize X and Y
    //fvec xmean(p), xsd(p);
    //Centering(Xtmp, xmean, xsd, n, p);
    /*fmat xmean = mean(X);
    for (i = 0; i < n; i++) {
  		for (j = 0; j < p; j++) {
  			X(i,j) -= xmean(j);
  		}	
    } */
    mat ymean = mean(Y);
    for (i = 0; i < n; i++) {
  		for (k = 0; k < K; k++) {
  			Y(i,k) -= ymean(k);
  		}	
    } 
  
	
    // containers
	mat E (n, K, fill::zeros);
    mat mu(p, K, fill::zeros) ;
    if (mu0.isNotNull())  {
		mu = as<mat>(mu0);
    }

    mat s (p, K, fill::ones) ;

    vec sigb(K, fill::ones);
	if (sigb0.isNotNull()) {
     	sigb = as<vec>(sigb0);
    }  else sigb = vectorise(var(Y, 0, 0)/2);

    vec sige(K, fill::ones);
    if (sige0.isNotNull())	{
    	sige = as<vec>(sige0);
    } else sige = sigb;

    mat Alpha(p, K);	
    if (Alpha0.isNotNull()) {
		Alpha = as<mat>(Alpha0);
    } else Alpha.fill(0.01);

	vec alpha(K);
	if (Alpha0.isNotNull()) {
		alpha = vectorise(mean(Alpha, 0));
    }  else alpha.fill(0.01);
   
    mat Beta(p, K, fill::zeros);
	if (mu0.isNotNull()) {
    	for (k=0; k < K; k++) {
			Beta.col(k) = Alpha.col(k) % mu.col(k);
		}
	}	

	

    // precomputation
    fvec xx = xtx_diag(Xtmp, n, p, xmean, xsd);

	fmat XY = MdotProd(Xtmp, Y, p, xmean, xsd);
	mat ytilde (n, K, fill::zeros);
	mat ytilde_j(n, K, fill::zeros);

	int it = 0;
	double L0 = -INFINITY;
	vec ell(maxit, fill::zeros);
    // algorithm
    while (it < maxit) {
		
		
         // E step 	
		for (k = 0; k < K; k++) {
			for (j = 0; j < p; j++) { 
			s (j, k) = sige(k)/(xx(j) + sige(k)/sigb(k)); 
			}
		}	
		
		for (k = 0; k < K; k++) {
			//ytilde.col(k) = X * Beta.col(k); 
			for (j = 0; j < p; j++) {
				char* col_j = Xtmp + j * n;
				//ytilde_j.col(k) = ytilde.col(k) - Beta(j, k) * X.col(j);
				ytilde_j.col(k) = ytilde.col(k) - VdotC(col_j, Beta(j, k), n, xmean(0, j), xsd[j]);
				//double temp = sum(X.col(j) % ytilde_j.col(k));
				double temp = VdotProd(col_j, ytilde_j.col(k), n, xmean(0, j), xsd[j]);
				mu(j, k) = (XY(j, k) - temp)/(xx(j) + sige(k)/sigb(k));
			
				// update Alpha
			    double ujk = log(alpha(k)/(1-alpha(k))) + 
				             0.5*(pow(mu(j,k), 2.0)/s(j,k) + log(s(j,k)/sigb(k)));				 
			    Alpha(j, k) = 1/(1+exp(-ujk));	
			    
				Beta(j, k)   = Alpha(j,k) * mu(j,k);
				// update ytilde
				ytilde.col(k) = ytilde_j.col(k) +  VdotC(col_j, Beta(j, k), n, xmean(0, j), xsd[j]);//X * (Alpha % mu.col(k));
		    }	
		}		

	   	// M step
		for (k = 0; k < K; k++) {
			// update sigma_{e_k}
		    E.col(k) = Y.col(k) - ytilde.col(k);
            vec temp = (Alpha.col(k) % (pow(mu.col(k), 2.0) + s.col(k)) - pow(Alpha.col(k) % mu.col(k), 2.0)) % xx;
			sige(k) = (sum(E.col(k) % E.col(k)) + sum(temp))/n;
			//if (sige(k) < 10E-6) sige(k) = 10E-6;
			// update sigma_{beta_k}
			sigb(k) = sum(Alpha.col(k) % (pow(mu.col(k), 2.0) + s.col(k)))/sum(Alpha.col(k));
			//if (sigb(k) < 10E-6) sigb(k) = 10E-6;
		    alpha(k) = sum(Alpha.col(k))/p; 
		}    
		
		ell(it) = LogLikInd(xx, E, p, K, mu, s, sige, sigb, Alpha, alpha);

		// convergence 
		if (ell(it) < L0){
            printf("Lowerbound decreasing, Error at iteration %d th iteration, diff=%g", it, ell(it) - L0);
            break;
        } else if ((ell(it) - L0) < epsStopLogLik) {
			printf("Converge at %d th iteration, LowerBound = %f \n", it, ell(it));
            break;
		}
		L0 = ell(it);

		it++;
    }
	if (it >= maxit) it = it - 1;
    vec LowerBound(it, fill::zeros);
    LowerBound = ell.subvec(0, it);   
    
    // intercepts
    for (int j = 0; j < p; ++j)
	{
		Beta.row(j) /= xsd[j];
	}
    mat Beta0(ymean - xmean*Beta);
    // returns
    List ret ;
    ret["N"] = n ;
    ret["mu"] = mu ;
	ret["s"]  = s;
	ret["Beta0"] = Beta0;
	ret["Beta"]  = Beta;
	ret["Alpha"] = Alpha ;
	ret["fdr"]   = 1 - Alpha;
    ret["sigb"]  = sigb;
	ret["sige"]  = sige;
	ret["alpha"] = alpha;
	ret["iter"]  = it+1;
	ret["eps"]   = ell(it) - ell(it-1); 
	ret["LowerBound"] = LowerBound;

    return(ret) ;
}




List emMultiplePlink (char* Xtmp, unsigned long long p,
                arma::mat Y,
            	Rcpp::Nullable<Rcpp::NumericMatrix> mu0     = R_NilValue,
               	Rcpp::Nullable<Rcpp::NumericVector> sigb0   = R_NilValue,
		 	    Rcpp::Nullable<Rcpp::NumericMatrix> Theta0  = R_NilValue,
			    Rcpp::Nullable<Rcpp::NumericVector> Lambda0 = R_NilValue,
		 	    Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0  = R_NilValue,
		 	    bool   fixlambda     = true, //if Lambda0 is not NULL, fixlambda should be false
		        int    maxit         = 10^3, 
		        double epsStopLogLik = 10^-5
        ) {
    // inputs
    unsigned long long N = Y.n_rows ;
	const int K = Y.n_cols ;

	int i=0, j=0, k=0;
	
   
    //centeralize X and Y
    float* mean_x = new float[p];
    float* sd_x = new float[p];
    Centering(Xtmp, mean_x, sd_x, N, p);
    fmat xmean(mean_x, 1, p, false); 
    fvec xsd(sd_x, p, false);
    /*fmat xmean = mean(X);
    for (i = 0; i < N; i++) {
  		for (j = 0; j < p; j++) {
  			X(i,j) -= xmean(j);
  		}	
    } */
    mat ymean = mean(Y);
    for (i = 0; i < N; i++) {
  		for (k = 0; k < K; k++) {
  			Y(i,k) -= ymean(k);
  		}	
    } 
  

    // containers
	mat E (N, K, fill::zeros);

	mat mu(p, K, fill::zeros);
	if (mu0.isNotNull())  {
		mu = as<mat>(mu0);
    }

    rowvec mu_jNumer(K, fill::zeros);


    mat s (p, K, fill::ones) ;

    vec sigb(K, fill::ones);
    sigb = vectorise(var(Y, 0, 0)/2.0);
    if (sigb0.isNotNull()) {
     	sigb = as<vec>(sigb0);
    }
    
    mat Theta(K, K, fill::zeros);
	Theta.diag() = vectorise(2.0/var(Y, 0, 0));
    if (Theta0.isNotNull()) {
    	Theta = as<mat>(Theta0);
    } 
	
	mat Sigma(K, K, fill::zeros);

    // Alpha
    mat Alpha(p, K, fill::zeros);	
    Alpha.fill(0.2);
    if (Alpha0.isNotNull()) {
		Alpha = as<mat>(Alpha0);
    }
	
	vec alpha =  vectorise(mean(Alpha, 0));
	
	//Lambda
	vec Lambda(p, fill::ones);	
	bool fixLambda = false;
	double lambda = 1;
	if (Lambda0.isNotNull()) {
		Lambda = as<vec>(Lambda0);	
		lambda = mean(Lambda);
	} 
	if (fixlambda & (lambda == 1.0))
		fixLambda = true;	

    mat Beta(p, K, fill::zeros);
    // If mu_jk has initial values, then calculate. 
    if (mu0.isNotNull()) {
    	for (k=0; k < K; k++) {
			Beta.col(k) = Lambda % Alpha.col(k) % mu.col(k);
		}
	}	

    // precomputation
	fvec xx = xtx_diag(Xtmp, N, p, xmean, xsd);
	/*fvec xx(p); 
    for (j=0; j < p; j++) {
		 xx(j) = sum(X.col(j) % X.col(j));
	}   */
	
	// initialization
	int it = 0;
	double L0 = -INFINITY;
	vec ell(maxit, fill::zeros);
	//mat ytilde = X * Beta;
	mat ytilde = XdotB(Xtmp, Beta, N, xmean, xsd);
	mat y_j(N, K, fill::zeros);
	mat y_j_new(N, K, fill::zeros);
	mat residual_j(N, 1, fill::zeros);
 
    // algorithm
    while (it < maxit) {
         
		// E step
		double logitlambda = log(lambda/(1.0-lambda));
		vec logitalpha = log(alpha/(1.0-alpha));
			// update mu, s_jk and Alpha.
            for (j = 0; j < p; j++) {   
                //X_j = conv_to<vec>::from(X.col(j));
                char* col_j = Xtmp + j * N;
                //y_j = ytilde - Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
                for (k = 0; k < K; k++) {
                	y_j = ytilde - Lambda(j) * Vouter(Alpha.row(j)%mu.row(j), col_j, N, xmean(0, j), xsd[j]);

                    // update s_jk
			        s (j, k) = 1.0/(xx(j)*Theta(k,k) + 1.0/sigb(k));
                   
                    // update mu_jk
                    residual_j = (Y - y_j) * Theta.col(k);

                    //double mu12 = sum(X.col(j) % residual_j);
                    double mu12 = VdotProd(col_j, residual_j, N, xmean(0, j), xsd[j]);

                    double mu3 = mu_jkThree(Theta.row(k), Alpha.row(j), mu.row(j), xx(j), k);			   	
                    mu(j, k) = (mu12 - mu3) * s(j, k);

                    // update Alpha_jk
                    double v_jk = logitalpha(k) + Lambda(j)/2.0 * (mu(j,k)*mu(j,k)/s(j,k) + log(s(j,k)/sigb(k)));
                    Alpha(j, k) = 1.0 /(1.0+exp(-v_jk));
                    ytilde = y_j + Lambda(j) * Vouter(Alpha.row(j)%mu.row(j), col_j, N, xmean(0, j), xsd[j]);
				}
				
				// update Lambda_j
                if (!fixLambda) {
                	y_j_new = ytilde - Lambda(j) * Vouter(Alpha.row(j)%mu.row(j), col_j, N, xmean(0, j), xsd[j]);
                	for (k = 0; k < K; k++)	{
                        residual_j = (Y - y_j_new) * Theta.col(k);
                    	double mu12 = VdotProd(col_j, residual_j, N, xmean(0, j), xsd[j]);
                    	double mu3 = mu_jkThree(Theta.row(k), Alpha.row(j), mu.row(j), xx(j), k);			   	
                    	mu_jNumer(k) = (mu12 - mu3);
                    }
                    double u_j2 = sum(Alpha.row(j) % mu.row(j) % mu_jNumer) + 
                             sum(Alpha.row(j) % (-mu.row(j)%mu.row(j)/s.row(j) + log(s.row(j)/sigb.t())))/2.0;
                    double u_j3 = uThree(Theta, mu.row(j), Alpha.row(j), xx(j));
                	double u_j  = logitlambda + u_j2 + u_j3;
                	Lambda(j) = 1.0/(1.0+exp(-u_j));
          	    }
          	    ytilde = y_j + Lambda(j) * Vouter(Alpha.row(j)%mu.row(j), col_j, N, xmean(0, j), xsd[j]);
			}
        // M step

		for (k=0; k < K; k++) {
		 	E.col(k) = Y.col(k) - ytilde.col(k);
		}
		// update Theta^-1        		
		 for (k=0; k < K; k++) {
		 	vec Siama_kk2 = Lambda%Alpha.col(k)%(mu.col(k)%mu.col(k)+s.col(k)) - 
			                  Lambda%Lambda%Alpha.col(k)%Alpha.col(k)%mu.col(k)%mu.col(k);
			Sigma(k, k) = (sum(E.col(k) % E.col(k)) + sum(xx % Siama_kk2))/N;
			for (int t=k+1; t < K; t++) {
				vec Siama_kt2 = (Lambda - Lambda%Lambda) % Alpha.col(k) % mu.col(k) % Alpha.col(t) % mu.col(t);
			    Sigma(k, t) = (sum(E.col(k) % E.col(t)) + sum(xx % Siama_kt2))/N;
                Sigma(t, k) = Sigma(k, t);				
			}
		}
	    Theta = inv_sympd(Sigma);
       
		// update sigma_{beta_k}
		for (k = 0; k < K; k++) {
			vec la = Lambda % Alpha.col(k);
			vec mu_k2 = mu.col(k) % mu.col(k);
			sigb(k) = sum(la % (mu_k2 + s.col(k))) / sum(la) ;
			alpha(k) = mean(Alpha.col(k));
		} 
        if(!fixlambda) lambda = mean(Lambda);		
		
		
		ell(it) = LogLikMultiple(xx, E, p, K, mu, s, Theta, sigb, Lambda, lambda, Alpha, alpha);
		// convergence 
		if (ell(it) < L0){
            printf("Lowerbound decreasing, Error at iteration %d th iteration, diff=%g \n", it, ell(it) - L0);
            break;
        } else if ((ell(it) - L0) < epsStopLogLik) {
			printf("Converge at %d th iteration, LowerBound = %f \n", it, ell(it));
            break;
		}
		L0 = ell(it);
		
		it++;
    } 
	
	
	
    if (it >= maxit) it = it - 1;
    vec LowerBound(it, fill::zeros);
    LowerBound = ell.subvec(0, it);     
	
	// local FDR
    mat fdr(p, K, fill::zeros);
    fdr = 1 - Alpha.each_col() % Lambda;
	
	//for (k=0; k < K; k++) {
		Beta = (Lambda % Alpha.each_col()) % mu;
	//}
    // intercepts
	for (int j = 0; j < p; ++j)
	{
		Beta.row(j) /= xsd[j];
	}
	//Beta = diagmat(1/xsd) * Beta;

	mat Beta0(ymean - xmean*Beta);
    // returns
    List ret ;
    ret["N"]      = N;
    ret["mu"]     = mu;
	ret["s"]      = s;
	ret["Beta0"]  = Beta0;
	ret["Beta"]   = Beta;
	ret["Alpha"]  = Alpha;
	ret["alpha"]  = alpha;
	ret["Lambda"] = Lambda;
	ret["lambda"] = lambda;
	ret["fdr"]    = fdr;
    ret["sigb"]   = sigb;
	ret["Sigma"]  = Sigma;
	ret["Theta"]  = Theta; 
	ret["iter"]   = it+1;
	ret["eps"]    = ell(it) - ell(it-1); 
	ret["LowerBound"] = LowerBound;
    
    return(ret) ;
}



// [[Rcpp::export()]]
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
	clock_t t;
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
	long long size = (long long)N * (long long)P;

	char* Xtmp = new char[size];

	
	readPlink(stringname, N, P, Xtmp);
	//arma::Mat<char> X(Xtmp, N, P, false, false);
 	//X.replace(3, 0);
	//fvec Xf = conv_to<fvec>::from(X);
	//X.clear();
	//float* XX = &Xf[0];

	vector<string> snps = read_snpnames(bimfile, P);
	
	arma::mat Y = as<arma::mat>(pheno);
	int K = Y.n_cols;
	

	//cout << "I am still fine from this line!" << endl;
	// add your algorithm here for X and Y.
	List vi, vi_Ind;
	t = clock();
	if (fit) {
		if (IndAsInit) {
			vi_Ind = emIndPlink(Xtmp, P, Y, 
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

            vi = emMultiplePlink(Xtmp, P, Y, mu, sigb, Theta, 
		                    Lambda, Alpha, fixlambda, 
		                    maxit, epsStopLogLik);

		}	else {
			vi = emMultiplePlink(Xtmp, P, Y, mu0, sigb0, Theta0, 
		                    Lambda0, Alpha0, fixlambda, 
		                    maxit, epsStopLogLik);
		} 
  	}
  	delete[] Xtmp;
	// save the snp names.
	t = clock() - t;
	List output = List::create( 
						       Rcpp::Named("Y")    = Y,
							   Rcpp::Named("snps") = snps,
							   Rcpp::Named("vi")   = vi,
							   Rcpp::Named("vi_Ind")   = vi_Ind,
							   Rcpp::Named("time")   = ((float)t) / CLOCKS_PER_SEC );

    
	return output;
}


