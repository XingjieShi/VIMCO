#' @title
#' vimco
#' @description
#' Fit variational inference for mulitple correlated traits with variational Bayesion EM algorithm, all traits are jointly analyzed.
#'
#' @param X SNP genotypes recoded in terms of additive and dominant components, rows and columns are for n samples and p SNPs, respectively.
#' @param Y A phenotype matrix (n by K), each column is a  (quantitavie) phenotype.
#' @param mu0 Initial values for the variational mean of effect sizes (p by K).
#' @param sigb0 An initial vector (K by 1) for the standard deviation of effect size.
#' @param Theta0 An inital precision matrix (K by K) for errrors in multiple traits.
#' @param Lambda0 Inital values for variational means (p by 1) of nonzero indicators.
#' @param Alpha0 Initial values for the variational mean of nonzero indicators (p by K).
#' @param fixlambda a test parameter, always set to be true.
#' @param Maximum The maximal iteration in the VBEM algorithm.
#' @param epsStopLoglik The tolerance for convergence.
#' @return List of model parameters and local FDR estimates.
#' @examples
#'  data(datasim)
#'  fit_Mul = emMultiple(datasim$X, datasim$Y)
#' @details
#' \code{emMultiple} fits VIMCO for multiple correlated traits, all traits are modeled jointly.
#' @export

emMultiple <- function(X, Y, mu0 = NULL, sigb0 = NULL, Theta0 = NULL, Lambda0 = NULL, Alpha0 = NULL, fixlambda = TRUE, maxit = 10^3L, epsStopLogLik = 10^-5) {
  .Call(`_vimco_emMultiple`, X, Y, mu0, sigb0, Theta0, Lambda0, Alpha0, fixlambda, maxit, epsStopLogLik)
}
