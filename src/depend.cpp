   
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <bitset>
#include <boost/algorithm/string.hpp>
#include "depend.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;



double logpexp(double x) {
    double y = x;
	if(x < 16) y = log(1 + exp(x));
    return y;
}

vec  xtx_diag(char* Xtmp, int n, int p, mat xmean, vec xsd) {
    vec y(p, fill::zeros);
    for (int j = 0; j < p; ++j)	{
    	char* col_j = Xtmp + j * n;
    	
    	for (int i = 0; i < n; ++i)
    	{
    		double tmp = (col_j[i] - xmean(0, j))/xsd[j];
    	  	y(j) += tmp * tmp;
    	}
    }
    return y;
}

vec VdotC (char* x, double c, int n, double xmean_j, double xsd_j) {
    vec z(n);
    for (int i = 0; i < n; i++)	{
        z[i] = (x[i] - xmean_j)/xsd_j * c;
    }
    return z;
}


double VdotProd (char* x, mat y, int n, double xmean_j, double xsd_j) {
    double z = 0;
    for (int i = 0; i < n; i++)	{
    	double f = y(i, 0);
        z += (x[i] - xmean_j)/xsd_j * f;
    }
    return z;
}


mat MdotProd(char* Xtmp, mat Y, int p, mat xmean, vec xsd)
{
	int n = Y.n_rows;
    int K = Y.n_cols;
    mat XtY(p, K);
    for (int j = 0; j < p; j++) {
    	char* col_j = Xtmp + j * n;
        for (int k = 0; k < K; k++) {
            XtY(j, k) = VdotProd(col_j, Y.col(k), n, xmean(0, j), xsd[j]);
        }
    }
    return XtY;
}


mat XdotB(char* Xtmp, mat B, int n, mat xmean, vec xsd)
{	
	int p = B.n_rows;
    int K = B.n_cols;
    double sum;
    //mat X = reshape(x, n, p)

    mat XB(n, K);

    for (int i=0; i < n; ++i)	{
    	for (int k = 0; k < K; ++k)	{
    		sum = 0;
    		for (int j=0; j < p; ++j)	{
    			char X_ij = *(Xtmp + j * n + i);
    			sum += (X_ij - xmean(0, j))/xsd[j] * B(j, k);
    		}
    		XB(i, k) = sum;
    	}
    }

    return XB;
}


mat Vouter(rowvec x, char* y, int n, double mean_j, double sd_j)	{
	int K = x.n_elem;
	mat outer(n, K, fill::zeros);
	for (int k = 0; k < K; ++k)
	{
		for (int i = 0; i < n; ++i)
		{
			outer(i, k) = x[k] * (y[i] - mean_j)/sd_j;
		}
	}
	return outer;
}

void Centering(char* Xtmp, double* meanOut, double* sdOut, int n, int p)	{
	double value;
	vec mean(p), sd(p);
	for (int j = 0; j < p; ++j)	{
		mean[j] = 0;
		sd[j] = 0;
    	char* col_j = Xtmp + j * n;
    	for (int i = 0; i < n; ++i)
    	{
    		value = col_j[i];
    		mean[j] += value;
    		sd[j] += value * value;
    	}
    	meanOut[j] = mean[j]/n;
    	sdOut[j] = sqrt(sd[j]/n - meanOut[j] * meanOut[j]);

    }	
}



double uThree(mat& Theta, rowvec mu_j, rowvec Alpha_j, double XtX_jj){
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




double mu_jkThree(rowvec theta_k, rowvec Alpha_j, rowvec mu_j, double XtX_jj, int k){
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


double LogLikInd(vec& xx, mat& E, 
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



double LogLikMultiple(vec& xx, mat& E, 
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

