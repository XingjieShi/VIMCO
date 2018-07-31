# evaluate selection performance 
# similar to calFDR, but TP are counted differently.
# Taking the LD in genomes into account, when generate casual 
# SNPs, in each LD, no more than 1 snp can be casual.
# So numbers of groups can be counted, instead of individual SNPs.
#
calFDRsnp = function(Fdr, threshold, bet, group){
    p = nrow(bet)
    K = ncol(bet)
    act = 1*(bet != 0)
    sel = 1*(Fdr <= threshold)

    n_group = max(group)
    act_group = matrix(0, nrow=n_group, ncol=K)
    sel_group = matrix(0, nrow=n_group, ncol=K)

    for (j in 1:n_group) {
        for (k in 1:K) {
            act_group[j, k] = 1 * (sum(act[which(group == j), k]) != 0)
            sel_group[j, k] = 1 * (sum(sel[which(group == j), k]) != 0)
        }
    }

    TP  = sum(act_group * sel_group)
	pwr = TP/sum(act_group)
	
	#observed FDR.
	P = sum(sel_group)
	FDR = ifelse(P == 0, 0, (P - TP)/P)
	list(power = pwr,
	     FDR   = FDR)
}


# evaluate selection performance 
# with global Fdr control.
# threshold   --- expected glocal FDR.
calFDR = function(Fdr, threshold, bet){
    act = which(bet != 0)
    sel = which(Fdr <= threshold)
    TP  = length(intersect(act, sel))  


    pwr = TP/length(act)
    
    #observed FDR.
    P = length(sel)
    FDR = ifelse(P == 0, 0, (P - TP)/P)
    list(power = pwr,
         FDR   = FDR)
}



clockwise90 = function(a) t(a[nrow(a):1,])

# fractional Gussian Error in
# "Sparse Multivariate Regression With Covariance Estimation"
# JCGS 2010
FGE = function(H, p)
{
   sig = outer(1:p, 1:p, FUN = function(i,j) ((abs(i-j)+1)^(2*H) - 2*abs(i-j)^(2*H) + (abs(i-j)-1)^(2*H))/2)
   diag(sig)=1
   sig
}


AR = function(rho, p) {
    outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y))) 
}


## generate data with K > 1
library(mvtnorm)
rData = function(n, p, K, rhoX, lambda, alpha, ratio, type, rhoE, h)
{
    # correlation matrix of error E.
    if (type == "AR") { 
        sigma = AR(rhoE, K)
    }   else sigma = FGE(rhoE, K)

    sigx = AR(rhoX, p)
    x    = rmvnorm(n, mean=rep(0, p), sigx)
    X <- x2snp(x)
    #X <- x
    sigb = rep(sqrt(1), K)

    bet = matrix(0, nrow = p, ncol = K)
    eta = rbinom(p, 1, lambda) 
    gam = matrix(rbinom(p*K, 1, alpha), ncol=K)
    for (k in 1:K){
        bet[, k] = sigb[k] * rnorm(p) * gam[,k] * eta 
    }
    #image(clockwise90(bet), main = expression(beta))
    
	lp  = X %*% bet
	sd.lp = diag(sqrt(diag(var(lp)) * (1/h - 1)))
    sige = sd.lp %*% sigma %*% sd.lp
    err = rmvnorm(n, rep(0, K), sige)
    Y   = lp + err
	
	## heritability is  0.5.
    heritability = diag(var(lp)/(var(lp)+var(err)))
 
    ## testing data
	xt = rmvnorm(n, mean=rep(0, p), sigx)
    Xt <- x2snp(xt)
    #Xt <- xt
    Yt = Xt %*% bet + rmvnorm(n, rep(0, K), sige)
	
	list(x  = X, 
	     xt = Xt,
	     y  = Y,
		 yt = Yt,
		 bet = bet,
		 sigx = sigx,
		 heritability = heritability)
}
	


library(mvtnorm)
rData_nH = function(n, p, K, rhoX, lambda, alpha, ratio, type, rhoE, h)
{
    # correlation matrix of error E.
    if (type == "AR") { 
        sigma = AR(rhoE, K)
    }   else sigma = FGE(rhoE, K)

    sigx = AR(rhoX, p)
    x    = rmvnorm(n, mean=rep(0, p), sigx)
    X <- x2snp(x)
    #X <- x
    n1 <- round(p * K * lambda * alpha * (1 - ratio)/(1 + ratio))
    n2 <- round(ratio/(1 - ratio) * n1)

    gam <- matrix(0, nrow=p, ncol=K)

    index <- sample(1:p, n1 + n2, replace=FALSE)
    index_plei <- sample(index, n2, replace=FALSE)
    for (j in 1:p) {
        if (j %in% index) {
                if (j %in% index_plei) {
                    gam[j, sample(1:K, 2, replace=FALSE)] = 1
                }   else {
                    gam[j, sample(1:K, 1, replace=FALSE)] = 1
                }
        }
    }
    bet <- matrix(rnorm(p*K), ncol=K) * gam


    #image(clockwise90(bet), main = expression(beta))
    
    lp  = X %*% bet
    sd.lp = diag(sqrt(diag(var(lp)) * (1/h - 1)))
    sige = sd.lp %*% sigma %*% sd.lp
    err = rmvnorm(n, rep(0, K), sige)
    Y   = lp + err
    
    ## heritability is  0.5.
    heritability = diag(var(lp)/(var(lp)+var(err)))
 
    ## testing data
    xt = rmvnorm(n, mean=rep(0, p), sigx)
    Xt <- x2snp(xt)
    #Xt <- xt
    Yt = Xt %*% bet + rmvnorm(n, rep(0, K), sige)
    
    list(x  = X, 
         xt = Xt,
         y  = Y,
         yt = Yt,
         bet = bet,
         sigx = sigx,
         heritability = heritability)
}
    
# group all SNPs in disjoint blocks.
blockSNP = function(X, tau){
    p = ncol(X)
    index = 1:p
    cor_X <- cor(X)

    B <- 1:p
    S = numeric(0)
    tag <- vector("list", p)
    i = 1
    while (length(B) != 0)  {
            j <- min(B)
            tag[[i]] <- B[which(cor_X[j, B] > tau)]
            B = setdiff(B, tag[[i]])
            i = i + 1
    }
    g = tag[1:(i-1)]
    group <- numeric(p)
    for(k in 1:(i-1)) group[g[[k]]] = k
    group
}


# control the pleiotripy proportion in beta, based on binomial model.
# let x ~ B(p, K). 
# if Pr(x>1)/Pr(x>0) = g, find parameter p.
plei2pr = function(g, K) {
    if (g < 0.001) {
        stop("g should be larger than 0.001")
    }
    f = function(pr) {
        (1 - dbinom(0, K, pr) - dbinom(1, K, pr))/(1 - dbinom(0, K, pr)) - g
    }
    uniroot(f, c(0.001, 1))$root
}


rGwas = function(X, K, alpha, g, tau, type, rhoE, h)
{
    # correlation matrix of error E.
    if (type == "AR") { 
        sigma = AR(rhoE, K)
    }   else sigma = FGE(rhoE, K)


    n <- nrow(X)
    p <- ncol(X)
    group = blockSNP(X, tau)

    G <- max(group)
    n_nz <- p*alpha*K
    gam <- matrix(0, nrow=p, ncol=K)

    # g is proportian of pleiotripy, pr is the binomial parameters
    if (g != 0) {
        pr = plei2pr(g, K)
        # index of candidate groups
        ind_g <- sample(1:max(group), n_nz/K/pr, replace = FALSE)
        # index of candidate casual SNPs, only one per group.
        ind_snp <- numeric(length(ind_g))
        for (i in 1:length(ind_g)) ind_snp[i] <- which(group == ind_g[i]) [1]
     
        
        tmp = sample(1:(K*length(ind_g)), n_nz, replace = F)
        gam[ind_snp, ][tmp] = 1
    } else {
        ind = sample(p*K, n_nz, replace = FALSE)
        gam[ind] = 1
    }

    bet <- matrix(rnorm(p*K), ncol=K) * gam
    lp  = X %*% bet
    sd.lp = diag(sqrt(diag(var(lp)) * (1/h - 1)))
    sige = sd.lp %*% sigma %*% sd.lp
    err = rmvnorm(n, rep(0, K), sige)
    Y   = lp + err
    
    ## heritability is  0.5.
    heritability = diag(var(lp)/(var(lp)+var(err)))
    ## testing data
    #xt = X
    #Xt <- x2snp(xt)
    ##Xt <- xt
    Yt = X %*% bet + rmvnorm(n, rep(0, K), sige)
    
    list(x  = X, 
         group = group,
         #xt = Xt,
         y  = Y,
         yt = Yt,
         bet = bet,
         heritability = heritability)
}


# transform Normal x to SNPs.
x2snp = function (X)
{
    n <- nrow(X)
    p <- ncol(X)
    maf <- runif(p, 0.05, 0.5)
    AAprob = maf^2;
    Aaprob = 2 * maf * (1-maf);
    quanti = cbind(1-Aaprob-AAprob, 1- AAprob)  ## attention
    snp = matrix(0, n, p);
    for (j in 1:p){
        cutoff = qnorm(quanti[j,]);
        snp[X[,j] <  cutoff[1], j] = 0;
        snp[X[,j] >= cutoff[1]  & X[,j] < cutoff[2],j] = 1;  ## attention
        snp[X[,j] >= cutoff[2], j] = 2;
    }
    snp
}

# prepare binary ped files for GEMMA/Plink.
# imput:
#       x          --- simulated design matrix, enties are 
#                      counts of minor alleles.
#       y          --- mulvariate outcomes
#       stringname --- prefix for filenames
#       plink_exe  --- location of plink.exe
# output:
#       write files fro analysis in Gemma/Plink.
ad2bed = function(x, y, stringname, plink_exe="/home/sysadmin/Software/plink/plink")
{   
    n <- nrow(x)
    p <- ncol(x)
    snp <- ad2mm(x) # its body is in vb_aux.cpp.
    tped <- cbind(rep(1, p), 
                  paste0("snp", 1:p), 
                  rep(0, p),
                  rep(0, p),
                  snp
                 )
    write.table(tped, quote = FALSE, 
                row.names = FALSE, col.names = FALSE,
                paste0(stringname, ".tped"))
    tfam = data.frame(paste0("ID", 1:n), 
                      paste0("ID", 1:n), 
                      paste0("PD", 1:n),  
                      paste0("MD", 1:n),  
                      rep(1, n),
                      rep(-9,n))
    write.table(tfam, quote = FALSE, 
                row.names = FALSE, col.names = FALSE,
                paste0(stringname, ".tfam"))
    cmd = paste0(plink_exe, " --tfile ", stringname, " --noweb --make-bed --out ", stringname)
    system(cmd)

    # fill y in ".fam" files
    famfile <- read.table(paste0(stringname, ".fam"), 
                          stringsAsFactors=FALSE, header=FALSE)
    famfile$V6 <- NULL
    famfileK <- cbind(famfile, y)
    write.table(famfileK, quote = FALSE, 
                row.names = FALSE, col.names = FALSE,
                paste0(stringname, ".fam"))

}



# fit Gemma, both Multivariate and Single analysis.
# imput:
#       x          --- simulated normal design matrix
#       y          --- mulvariate quantitative outcomes
#       stringname --- prefix for filenames
# output:
#       a list contained multivariate and single analysis produced by Gemma.
fitGemma = function(x, y, stringname, plink_exe="/home/sysadmin/Software/plink/plink",
                    tool_exe="~/gemma/bin/gemma")
{
    K = ncol(y)
    # prepare data
    ad2bed(x, y, stringname, plink_exe)

    # generate relatedness matrix
    cmd0 <- paste0(tool_exe, " -bfile ", stringname, " -gk 1 -o ", stringname)
    system(cmd0)
   
    # multivariate LMM
    # bug: https://groups.google.com/forum/#!topic/gemma-discussion/uC1YciiSq8w
    # skip bugs.
    cmd1 <- paste0(tool_exe, " -bfile ", stringname, " -k ./output/", stringname, 
                    ".cXX.txt -lmm 2 -n ", paste(1:K, collapse=" "), " -o ", stringname)
    #err = tryCatch(system2(cmd1, stdout=TRUE, stderr=TRUE), error = function(e) NULL)
    #err = tryCatch(system(cmd1, ignore.stdout=F, ignore.stderr=F), error = function(e) NULL)
    
    res <- try(system(cmd1), silent = TRUE)
    
    #if (is.null(err)) {
    if (res == 22)  {
        #list(res_mlm = NA,
        #     res_slm = NA)
        stop("something wrong here")
    }   else    {
        res_mlm <- read.delim(paste0("./output/", stringname, ".assoc.txt"))
        # single LMM
        res_slm <- vector("list", K)
        for (k in 1:K) {
            cmd <- paste0(tool_exe, " -bfile ", stringname, " -k ./output/", stringname, 
                             ".cXX.txt -lmm 2 -n ", k, " -o ", paste0(stringname, k))
                system(cmd)
            res_slm[[k]] <- read.delim(paste0("./output/", stringname, k, ".assoc.txt"))
         }

        list(res_mlm = res_mlm,
             res_slm = res_slm)
    }
}



