#' @title
#' vimco
#' @description
#' Fit Bayesion variable selection regression with variational Bayesion EM algorithm, for each trait seperately.
#'
#' @param X SNP genotypes recoded in terms of additive and dominant components, rows and columns are for n samples and p SNPs, respectively.
#' @param Y A phenotype matrix (n by K), each column is a  (quantitavie) phenotype.
#' @param mu0 Initial values for the variational mean of effect sizes (p by K).
#' @param Alpha0 Initial values for the variational mean of nonzero indicators (p by K).
#' @param sigb0 An initial vector for the standard deviations of effect size.
#' @param sige0 An initial vector for the standard deviations of random errors in all trait.
#' @param maxit The maximal iteration in the VBEM algorithm.
#' @param epsStopLoglik A tolerance for convergence.
#' @return List of model parameters and local FDR estimates.
#' @examples
#'  data(datasim)
#'  fit_Ind = emInd(datasim$X, datasim$Y)
#' @details
#' \code{emInd} fits BVSR for multiple traits, but each trait is modeled seperately.
#' @export

emInd <- function(X, Y, mu0 = NULL, Alpha0 = NULL, sigb0 = NULL, sige0 = NULL, maxit = 10^3L, epsStopLogLik = 10^-5) {
  .Call(`_vimco_emInd`, X, Y, mu0, Alpha0, sigb0, sige0, maxit, epsStopLogLik)
}
