library(Rcpp)
sourceCpp("vb_aux.cpp") 
source("functions.r")   # used for prerformance evaluation

assig = function(n_args){
    # example:  
	# n_args = c(2, 3, 4)
	cargs <- vector("list", length(n_args))
	for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
	t(expand.grid(cargs))
}
n_args = c(3, 3, 3)
jobs = assig(n_args)

	
gFDR <- 0.1

	type = c("AR", "FGE")[1]
	if(type == "AR")  rho <- c(0.2, .5, .8)
    if(type == "FGE") rho <- c(0.5, 0.9, 0.96)



main <- function(number0) 
{
	number <- as.numeric(number0)

	id = jobs[, number]
    #id <- c(1, 2, 3)
	
	g <- c(0, .15, .3)[id[1]]
	
	rhoX <- c(0)[id[2]]
	
    rhoE <- rho[id[3]]

	tau = 0.6
	alpha  = 0.01
	h <- 0.3

		path = paste("h", h, "g", g, "rhoE", rhoE, "rhoX", rhoX, sep="_")
		setwd(paste("./", path, sep=""))

		file_names = list.files()

    	case = file_names[grep(path, file_names)]

    	R <- length(case)
    	FDR_Ha = PWR_Ha = matrix(0, nrow = R, ncol=4)
    	FDR_Hb = PWR_Hb = matrix(0, nrow = R, ncol=4)
    	COR  = matrix(0, nrow = R, ncol=2)
    	AUC_Ha = AUC_Hb<- matrix(0, nrow = R, ncol=4)
    	for (r in 1:R)	
    	{
    		res = readRDS(case[r])
    		bet <- res$bet

    		p <- length(res$p_mlm)
    		K <- ncol(res$p_slm) 
			p_mlm <- matrix(0, nrow = p, ncol = K)
			for (k in 1:K)	{
				p_mlm[, k] <- res$p_mlm
			}

    		COR[r, ] = c(res$cor_Mul, res$cor_Ind)

    		AUC_Ha[r, ] = c(calAUC(abs(bet), 1 - c(res$fit_Mul$fdr)),
						 calAUC(abs(bet), 1 - c(res$fit_Ind$fdr)),
					 	 calAUC(abs(bet), 1 - p_mlm),
					 	 calAUC(abs(bet), 1 - res$p_slm)
						)	
    		
    		# Ha: \beta_jk = 0
			Ha_Mul <- calFDRsnp(fdr2Fdr(res$fit_Mul$fdr), gFDR,      res$bet, res$group) 	
			Ha_Ind <- calFDRsnp(fdr2Fdr(res$fit_Ind$fdr), gFDR,      res$bet, res$group) 	
			Ha_Slm <- calFDRsnp(res$p_slm,                0.05/p/K,  res$bet, res$group) 
			Ha_Mlm <- calFDRsnp(p_mlm,                    0.05/p,    res$bet, res$group) 

 			FDR_Ha[r, ] = c(Ha_Mul$FDR,   Ha_Ind$FDR,   Ha_Mlm$FDR,   Ha_Slm$FDR)
			PWR_Ha[r, ] = c(Ha_Mul$power, Ha_Ind$power, Ha_Mlm$power, Ha_Slm$power)


			# Hb: \beta_j = 0
			bet_Hb <- 1 * (rowSums(res$bet != 0) != 0)
   			fdr_Mul_Hb <- apply(res$fit_Mul$fdr, 1, "min")
			fdr_Ind_Hb <- apply(res$fit_Ind$fdr, 1, "min")
			p_slm_Hb <- apply(res$p_slm, 1, "min")

			Hb_Slm <- calFDRsnp(data.matrix(p_slm_Hb),  0.05/p, data.matrix(bet_Hb), res$group)
   			Hb_Mlm <- calFDRsnp(data.matrix(res$p_mlm), 0.05/p, data.matrix(bet_Hb), res$group)

 			Hb_Mul <- calFDRsnp(data.matrix(fdr_Mul_Hb), gFDR,  data.matrix(bet_Hb), res$group)
   			Hb_Ind <- calFDRsnp(data.matrix(fdr_Ind_Hb), gFDR,  data.matrix(bet_Hb), res$group)

 			FDR_Hb[r, ] = c(Hb_Mul$FDR,   Hb_Ind$FDR,   Hb_Mlm$FDR, Hb_Slm$FDR)
 			PWR_Hb[r, ] = c(Hb_Mul$power, Hb_Ind$power, Hb_Mlm$power, Hb_Slm$power)


 			AUC_Hb[r, ] = c(calAUC(abs(bet_Hb), 1 - fdr_Mul_Hb),
						 	calAUC(abs(bet_Hb), 1 - fdr_Ind_Hb),
					 	 	calAUC(abs(bet_Hb), 1 - res$p_mlm),
					 	 	calAUC(abs(bet_Hb), 1 - p_slm_Hb)
						   )

		}	

		measure = list()
		measure$PWR_Ha = colMeans(PWR_Ha)
		measure$PWR_Hb = colMeans(PWR_Hb)	
   		measure$FDR_Ha = colMeans(FDR_Ha)
		measure$FDR_Hb = colMeans(FDR_Hb)
		measure$COR    = colMeans(COR)
		measure$AUC_Ha = AUC_Ha
		measure$AUC_Hb = AUC_Hb 
 

		setwd("../")

	   	saveRDS(measure, 
				paste(path, ".rds", sep="_"))

	   	
}


args <- commandArgs(TRUE)
main(args[1])




