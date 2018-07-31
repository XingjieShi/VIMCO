# simulation for VIMCO
# AR  
# 


library(vimco)
library(Rcpp)
library(mvtnorm)

sourceCpp("vb_aux.cpp") 
source("functions.r")   # used for prerformance evaluation

assig = function(n_args){
    # example:  
	# n_args = c(2, 3, 4)
	cargs <- vector("list", length(n_args))
	for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
	t(expand.grid(cargs))
}
n_args = c(3, 3, 3, 50)
jobs = assig(n_args)



main <- function(number0) {

	type = c("AR", "FGE")[1]
	if(type == "AR")  rho <- c(0.2, .5, .8)
    if(type == "FGE") rho <- c(0.5, 0.9, 0.96)

	number <- as.numeric(number0)
	set.seed(number*100)	
	id = jobs[, number]
	g <- c(0, .15, .3)[id[1]]
	
	rhoX <- c(.2, .5, .8)[id[2]]
	
    rhoE <- rho[id[3]]

    r <- id[4]
	

	n = 5000
	p = 10000
	K = 4
	h = 0.3
	tau = 0.6
	alpha  = 0.01
	

	path = paste("h", h, "g", g, "rhoE", rhoE, "rhoX", rhoX, sep="_")
	ifelse(!dir.exists(path), dir.create(path), FALSE)
	setwd(paste("./", path, sep=""))

	sigx = AR(rhoX, p)



		#training data
		x    = rmvnorm(n, mean=rep(0, p), sigx)
		X <- x2snp(x)
		dat = rGwas(X, K, alpha, g, tau, type, rhoE, h) 
		x = dat$x
		group = dat$group
		y = dat$y
		bet = dat$bet	

		# fit Independent VB
		tic = proc.time()
		fit_Ind = emInd(x, y, maxit = 3000, epsStopLogLik = 10^-7)
		toc = proc.time()	
		print(toc - tic)
	
		# fit Multiple VB
		mu0     = fit_Ind$mu
		sigb0   = fit_Ind$sigb
		Theta0  = matrix(0, nrow=ncol(y), ncol=ncol(y))
		diag(Theta0)  =   1/c(fit_Ind$sige)
		Lambda0 = rep(1, p)
		Alpha0  = fit_Ind$Alpha 
		tic = proc.time()
		fit_Mul = emMultiple(x, y, mu0, 
							sigb0, Theta0, Lambda0, Alpha0,					 
							TRUE, 3000, 10^-7)
		toc = proc.time()
		print(toc - tic)
		

		# prediction
		yh_Mul <- cbind(1, x) %*% rbind(fit_Mul$Beta0, fit_Mul$Beta)	
		yh_Ind <- cbind(1, x) %*% rbind(fit_Ind$Beta0, fit_Ind$Beta)
		cor_Mul <- mean(diag(cor(dat$yt, yh_Mul)))
		cor_Ind <- mean(diag(cor(dat$yt, yh_Ind)))

    	# fit Gemma
  	    stringname <- paste("sim", r, sep="_")
  	    # plink and GEMMA path in your compupter or HPC.
	    plink_exe="/home/sysadmin/Software/plink/plink"
 	    tool_exe="~/gemma/bin/gemma"
 	    fit_Gemma <- fitGemma(x, y, stringname, plink_exe, tool_exe)
 	    fit_slm <- fit_Gemma$res_slm
   	    fit_mlm <- fit_Gemma$res_mlm

  		p_slm <- matrix(0, nrow = p, ncol = K)
		for (k in 1:K)	{
			p_slm[, k] <- fit_slm[[k]]$p_lrt 
		}
		p_mlm <- fit_mlm$p_lrt

  	    system(paste0("rm ", stringname, "*"))	
    	system(paste0("rm ", "output/", stringname, "*"))	



    	# output 
		res = list()
		res$bet     = bet
		res$group   = group
		res$sigx    = sigx
		res$fit_Mul = fit_Mul
		res$cor_Mul = cor_Mul
		res$fit_Ind = fit_Ind
		res$cor_Ind = cor_Ind
		res$p_mlm   = p_mlm
		res$p_slm 	= p_slm	
		saveRDS(res, 
				paste("h", h, "g", g, "rhoE", rhoE, "rhoX", rhoX, r, ".rds", sep="_"))
		setwd("../")

}
	

args <- commandArgs(TRUE)
main(args[1])

	


