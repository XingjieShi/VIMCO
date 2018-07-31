# include "vb.hpp"
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
using namespace Rcpp;
using namespace arma;


double uThree(mat& Theta, rowvec mu_j, rowvec Alpha_j, float XtX_jj){
	int K = Theta.n_cols;
	double ret = 0; 
	for (int s=0; s < K; s++) {
		double Alpha_js  = Alpha_j(s);
		double mu_js     = mu_j(s);
		for (int t=0; t < K; t++) {
			double theta_st = Theta(s, t);
			double Alpha_jt = Alpha_j(t);
			double mu_jt    = mu_j(t);
		    if (t != s) ret = ret + theta_st * Alpha_js * mu_js * Alpha_jt * mu_jt;
	    }
	}
	return ret * XtX_jj/2.0;
}




double mu_jkThree(rowvec theta_k, rowvec Alpha_j, rowvec mu_j, float XtX_jj, int k){
	int K = theta_k.n_elem;
	double ret = 0;
	for (int t=0; t < K; t++) {
		double theta_kt = theta_k(t);
		double Alpha_jt = Alpha_j(t);
		double mu_jt    = mu_j(t);
		if (t != k) ret = ret + theta_kt * Alpha_jt * mu_jt;
	}
	return ret * XtX_jj;
}



double LogLikMultiple(fvec& xx, mat& E, 
           int p,      int    K,
           mat& mu,     mat&    s, 
           mat& Theta,  vec&    sigb, 
		   vec& Lambda, double lambda,
		   mat& Alpha,  vec alpha) {
	vec term(7, fill::zeros);
	
	int k, t, j;
	// term 1
	double temp0 = 0;
    for (k=0; k < K; k++)
		for (int t=0; t < K; t ++)
            temp0 = temp0 + Theta(k, t) * sum(E.col(k) % E.col(t));    			
    term(0) = -1.0/2*temp0;
	
	// term 2
	double temp1 = 0;
	for (k=0; k < K; k++) 
		for (j=0; j < p; j++)
			temp1 = temp1 + Theta(k, k)*xx(j)*
		                  ( Lambda(j)*Alpha(j,k)*(mu(j,k)*mu(j,k) + s(j, k)) -
   						  Lambda(j)*Lambda(j)*Alpha(j,k)*Alpha(j,k)*mu(j,k)*mu(j,k) );
	term(1) = -1.0/2*temp1;
	
	// term 3
	double temp2 = 0;
	for (k = 0; k < K; k++) {
		for (t = 0; t < K; t++) {
            double temp2_st = 0;
			for (j = 0; j < p; j++) {
			    temp2_st = temp2_st + 
				    xx(j)*(Lambda(j) - Lambda(j)*Lambda(j))*Alpha(j, k)*Alpha(j, t)*mu(j, k)*mu(j, t);
			}
			if (t != k) temp2 = temp2 + Theta(k, t)*temp2_st;
		}		
	}	
	term(2) = -1.0/2*temp2;
	
	// term 4
	vec temp31(p, fill::zeros);
	temp31 =  Lambda%log(Lambda + (Lambda==0)) + 
				(1.0-Lambda)%log((1.0-Lambda) + (Lambda==1.0)) ;
	term(3) = - sum(temp31) + log(lambda/(1.0-lambda+(lambda==1.0)))*sum(Lambda-1.0) + p*log(lambda);
	
	//term 5
	mat temp41(p, K, fill::zeros);
	temp41 = Alpha%log(Alpha + (Alpha==0)) + 
	        (1.0-Alpha)%log((1.0-Alpha) + (Alpha==1.0)); 
    term(4) = - accu(temp41) + sum(log(alpha/(1.0-alpha +(alpha==1.0))) % vectorise(sum(Alpha-1.0, 0))) + p*sum(log(alpha));
	
	//term6
	term(5) = E.n_rows*log(det(Theta))/2.0;
	
	//term7
	double temp6=0;
	for (j = 0; j < p; j++) {
		for (k = 0; k < K; k++)
			temp6 = temp6 + Lambda(j)*Alpha(j,k)*(1.0 + log(s(j,k)/sigb(k)) - (mu(j,k)*mu(j,k) + s(j, k))/sigb(k));
		
	}
	term(6) = temp6/2.0;
	
	return sum(term);	
    	
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
