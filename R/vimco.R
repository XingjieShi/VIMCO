#' @title
#' vimco
#' @description Fit VIMCO for a GWAS data which is in a plink format (e.g. sim.bed, sim.bim, sim.fam).
#' @param stringname The plink file prefix.
#' @param nPheno The number of phenotypes.
#' @param fit whether fit (TRUE) model or not (FALSE).
#' @param IndAsInit whether take VBSR as initial values.
#' @param mu0 Initial values for the variational mean of effect sizes (p by K).
#' @param sigb0 An initial vector (K by 1) for the standard deviation of effect size.
#' @param Theta0 An inital precision matrix (K by K) for errrors in multiple traits.
#' @param Lambda0 Inital values for variational means (p by 1) of nonzero indicators.
#' @param Alpha0 Initial values for the variational mean of nonzero indicators (p by K).
#' @param fixlambda a test parameter, always set to be true.ı
#' @param maxit The maximal iteration in the VBEM algorithm.
#' @param epsStopLogLik The tolerance for convergence.
#'
#' @return List of model parameters and local FDR estimates.'
#' @examples
#' path <- system.file(package = "vimco")
#' setwd(path)
#' stringname <- "sim"
#' tmp <- vimco(stringname, nPheno = 4, fit = TRUE)
#' @details
#' \code{vimco} is the main function used for real GWAS analysis. If no intial values are provided and IndAsInit is set to be TRUE, it will fit BVSR first, and use its results as initial values for fitting VIMCO.
#' For a GWAS with n = 5000, p = 30000, it usually takes a few hours to finish the whole inference.
#' @references
#' VIMCO: Variational Inference for Multiple Correlated Outcomes in Genome-wide Association Studies[https://arxiv.org/abs/1807.10467]

vimco <- function(stringname, nPheno, fit = TRUE, IndAsInit = TRUE, mu0 = NULL, sigb0 = NULL, Theta0 = NULL, Lambda0 = NULL, Alpha0 = NULL, fixlambda = TRUE, maxit = 10^3L, epsStopLogLik = 10^-5) {
  .Call(`_vimco_vimco`, stringname, nPheno, fit, IndAsInit, mu0, sigb0, Theta0, Lambda0, Alpha0, fixlambda, maxit, epsStopLogLik)
}

