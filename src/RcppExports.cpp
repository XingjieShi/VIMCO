// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// emInd
List emInd(arma::mat X, arma::mat Y, Rcpp::Nullable<Rcpp::NumericMatrix> mu0, Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0, Rcpp::Nullable<Rcpp::NumericVector> sigb0, Rcpp::Nullable<Rcpp::NumericVector> sige0, int maxit, double epsStopLogLik);
RcppExport SEXP _vimco_emInd(SEXP XSEXP, SEXP YSEXP, SEXP mu0SEXP, SEXP Alpha0SEXP, SEXP sigb0SEXP, SEXP sige0SEXP, SEXP maxitSEXP, SEXP epsStopLogLikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Alpha0(Alpha0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type sigb0(sigb0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type sige0(sige0SEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsStopLogLik(epsStopLogLikSEXP);
    rcpp_result_gen = Rcpp::wrap(emInd(X, Y, mu0, Alpha0, sigb0, sige0, maxit, epsStopLogLik));
    return rcpp_result_gen;
END_RCPP
}
// emMultiple
List emMultiple(arma::mat X, arma::mat Y, Rcpp::Nullable<Rcpp::NumericMatrix> mu0, Rcpp::Nullable<Rcpp::NumericVector> sigb0, Rcpp::Nullable<Rcpp::NumericMatrix> Theta0, Rcpp::Nullable<Rcpp::NumericVector> Lambda0, Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0, bool fixlambda, int maxit, double epsStopLogLik);
RcppExport SEXP _vimco_emMultiple(SEXP XSEXP, SEXP YSEXP, SEXP mu0SEXP, SEXP sigb0SEXP, SEXP Theta0SEXP, SEXP Lambda0SEXP, SEXP Alpha0SEXP, SEXP fixlambdaSEXP, SEXP maxitSEXP, SEXP epsStopLogLikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type sigb0(sigb0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Theta0(Theta0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type Lambda0(Lambda0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Alpha0(Alpha0SEXP);
    Rcpp::traits::input_parameter< bool >::type fixlambda(fixlambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsStopLogLik(epsStopLogLikSEXP);
    rcpp_result_gen = Rcpp::wrap(emMultiple(X, Y, mu0, sigb0, Theta0, Lambda0, Alpha0, fixlambda, maxit, epsStopLogLik));
    return rcpp_result_gen;
END_RCPP
}
// readPlink
Rcpp::List readPlink(std::string stringname, int nPheno);
RcppExport SEXP _vimco_readPlink(SEXP stringnameSEXP, SEXP nPhenoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type stringname(stringnameSEXP);
    Rcpp::traits::input_parameter< int >::type nPheno(nPhenoSEXP);
    rcpp_result_gen = Rcpp::wrap(readPlink(stringname, nPheno));
    return rcpp_result_gen;
END_RCPP
}
// vimco
Rcpp::List vimco(std::string stringname, int nPheno, bool fit, bool IndAsInit, Rcpp::Nullable<Rcpp::NumericMatrix> mu0, Rcpp::Nullable<Rcpp::NumericVector> sigb0, Rcpp::Nullable<Rcpp::NumericMatrix> Theta0, Rcpp::Nullable<Rcpp::NumericVector> Lambda0, Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0, bool fixlambda, int maxit, double epsStopLogLik);
RcppExport SEXP _vimco_vimco(SEXP stringnameSEXP, SEXP nPhenoSEXP, SEXP fitSEXP, SEXP IndAsInitSEXP, SEXP mu0SEXP, SEXP sigb0SEXP, SEXP Theta0SEXP, SEXP Lambda0SEXP, SEXP Alpha0SEXP, SEXP fixlambdaSEXP, SEXP maxitSEXP, SEXP epsStopLogLikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type stringname(stringnameSEXP);
    Rcpp::traits::input_parameter< int >::type nPheno(nPhenoSEXP);
    Rcpp::traits::input_parameter< bool >::type fit(fitSEXP);
    Rcpp::traits::input_parameter< bool >::type IndAsInit(IndAsInitSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type sigb0(sigb0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Theta0(Theta0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type Lambda0(Lambda0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type Alpha0(Alpha0SEXP);
    Rcpp::traits::input_parameter< bool >::type fixlambda(fixlambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsStopLogLik(epsStopLogLikSEXP);
    rcpp_result_gen = Rcpp::wrap(vimco(stringname, nPheno, fit, IndAsInit, mu0, sigb0, Theta0, Lambda0, Alpha0, fixlambda, maxit, epsStopLogLik));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_vimco_emInd", (DL_FUNC) &_vimco_emInd, 8},
    {"_vimco_emMultiple", (DL_FUNC) &_vimco_emMultiple, 10},
    {"_vimco_readPlink", (DL_FUNC) &_vimco_readPlink, 2},
    {"_vimco_vimco", (DL_FUNC) &_vimco_vimco, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_vimco(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
