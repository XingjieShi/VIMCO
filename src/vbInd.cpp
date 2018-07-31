# include "vb.hpp"
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
using namespace Rcpp;
using namespace arma;

double logpexp(double x) {
    double y = x;
	if(x < 16) y = log(1 + exp(x));
    return y;
}


float VdotProd (fvec x, fvec y, uword n) {
    float z = 0;
    for (uword i = 0; i < n; i++)
        z += x[i] * y[i];
    return z;
}


fmat MdotProd(fmat X, mat Y)
{
    uword N = X.n_rows;
    uword p = X.n_cols;
    uword K = Y.n_cols;
    fmat XtY(p, K);
    for (uword j = 0; j < p; j++) {
        for (uword k = 0; k < K; k++) {
        	fvec a = conv_to<fvec>::from(X.col(j));
        	fvec b = conv_to<fvec>::from(Y.col(k));
            XtY(j, k) = VdotProd(a, b, N);
        }
    }
    return XtY;
}



double LogLikInd(fvec& xx, mat& E, 
           int p,      int    K,
           mat& mu,     mat&    s, 
           vec& sige,   vec&    sigb, 
		   mat& Alpha,  vec& alpha) {
	vec term(6, fill::zeros);
	
	int k,  j;
	
	// term 1  increasing
	double temp0 = 0;
    for (k=0; k < K; k++)
            temp0 = temp0 +  sum(E.col(k) % E.col(k))/sige(k);    			
    term(0) = -1.0/2*temp0;
	
	// term 2  decreasing
	double temp1 = 0;
	for (k=0; k < K; k++) 
		for (j=0; j < p; j++)
			temp1 = temp1 + xx(j)/sige(k)*
		                  (Alpha(j,k) * (mu(j,k)*mu(j,k) + s(j, k)) -
   						   Alpha(j,k)*Alpha(j,k)*mu(j,k)*mu(j,k));
	term(1) = -1.0/2*temp1;
	
	// term 3  decreasing
	mat temp2(p, K, fill::zeros);
	vec logodd = log(alpha/(1-alpha)); 
	for (k = 0; k < K; k++)
		for (j = 0; j < p; j++)
	        temp2(j, k) =  (Alpha(j, k) - 1)*logodd(k) - logpexp(-logodd(k)) ;
	
    term(2) = accu(temp2);
	
	//term 4 increasing
	mat temp3(p, K, fill::zeros);
	for (k = 0; k < K; k++)
		for (j = 0; j < p; j++)
	        temp3(j, k) =  Alpha(j, k)*log(Alpha(j, k)+ (Alpha(j, k)==0)) + 
							(1-Alpha(j, k))*log(1-Alpha(j, k) + (Alpha(j, k)==1));
  
    term(3) = - accu(temp3);
	
	
	//term5  increasing
	term(4) = -1.0/2 * E.n_rows*sum(log(sige));
	
	//term6  decreasing
	double temp5=0;
	for (j = 0; j < p; j++) {
		for (k = 0; k < K; k++)
			temp5 = temp5 + Alpha(j,k)*(1 + log(s(j,k)/sigb(k)) - (mu(j,k)*mu(j,k) + s(j, k))/sigb(k));
	}
	term(5) =  1.0/2*temp5;
	
	return sum(term);		   
}

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
	fmat XY = MdotProd(X, Y);
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
