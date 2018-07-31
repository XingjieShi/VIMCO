# simulation
# summary
# barplot
library(ggplot2)
library(cowplot)
fig_type = ".eps"

gFDR = 0.1
h <- 0.3

type = c("AR", "FGE")[1]
rho <- c(.2, .5, .8)
 
n_rho <- length(rho)    

for (i in 1:3) 
{
    for (j in 1:3) 
    {

        g <- c(0, 0.15, 0.3)[i]
	
        rhoX <- c(.2, .5, .8)[j]

        rhoE <- method <- character(0)
        AUC <- numeric(0)
        for (k in 1:n_rho) 
        {
            path = paste("h", h, "g", g, "rhoE", rho[k], "rhoX", rhoX, sep="_")
            res <- readRDS(paste0(path, "_.rds"))
            R <- nrow(res$AUC)
            method <- c(method, c("VIMCO", "BVSR", "sLMM"))
            rhoE <- c(rhoE, rep(rho[k], 3))
            AUC <- c(AUC, c(apply(res$AUC_Ha[,-3], 2, "mean")))
        }
        result =  data.frame(AUC = AUC, 
                             rhoE    = factor(rhoE, levels=as.character(rho)),
                             method = factor(method, levels=c("VIMCO", "BVSR", "sLMM")), 
                             stringsAsFactors = FALSE)



        fig_AUC <- ggplot(result, aes(x=rhoE, y=AUC, fill=method)) +

           geom_bar(position = "dodge", stat="identity") + 
            xlab(expression(rho[E])) + ylab("AUC") + theme_bw() +
              coord_cartesian(ylim=c(.7, 1)) + 

            scale_fill_manual(values=c("#1F78B4", "#FF7F00", "#6A3D9A"))+
            labs(fill='Method')  + theme_classic() + 
             theme(axis.title.x = element_text(size=20, face = "bold"),
                    axis.text.x = element_text(size=20, face = "bold"),
                    axis.title.y = element_text(size=20, face = "bold"),
                    axis.text.y = element_text(size=20,face = "bold"),
                    legend.position = c(0.88, 0.9),
                     legend.title = element_text(size=20,face = "bold"),
                     legend.text  = element_text(size=20,face = "bold"))
     

        
        method <- rhoE <- character(0)
        PWR_Ha <- numeric(0)
        FDR_Ha <- numeric(0)
        for (k in 1:n_rho) 
        {
            path = paste("h", h, "g", g, "rhoE", rho[k], "rhoX", rhoX, sep="_")
            res <- readRDS(paste0(path, "_.rds"))

            method <- c(method, "VIMCO", "BVSR")
            rhoE <- c(rhoE, rep(rho[k], 2))
            PWR_Ha <- c(PWR_Ha, res$PWR_Ha[1:2])
            FDR_Ha <- c(FDR_Ha, res$FDR_Ha[1:2])
        }

        # PWR
        result <- data.frame(method = factor(method, levels = c("VIMCO", "BVSR")),
                                    rhoE   = factor(rhoE,   levels=as.character(rho)),
                                       PWR_Ha = PWR_Ha,
                                        stringsAsFactors = FALSE)
        

        fig_PWR <- ggplot(result, aes(x=result$rho, y=PWR_Ha, fill=method)) +

             geom_bar(position = "dodge", stat="identity") + 

            xlab(expression(rho[E])) + ylab("Power") +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 0.832)) +
            scale_fill_manual(values=c("#1F78B4", "#FF7F00"))+
            labs(fill='Method')  + theme_classic() + 
               theme(axis.title.x = element_text(size=20, face = "bold"),
                    axis.text.x = element_text(size=20, face = "bold"),
                    axis.title.y = element_text(size=20, face = "bold"),
                    axis.text.y = element_text(size=20,face = "bold"),
                    panel.border = element_blank(),
                    legend.position = "none",
                     legend.title = element_text(size=20,face = "bold"),
                     legend.text  = element_text(size=20,face = "bold")) 



        # FDR
        result_FDR <- data.frame(method = factor(method, levels = c("VIMCO", "BVSR")),
                                 rhoE   = factor(rhoE,   levels=as.character(rho)),
                                 FDR_Ha = FDR_Ha,
                                 stringsAsFactors = FALSE)
    

        fig_FDR <- ggplot(result_FDR, aes(x=result$rho, y=FDR_Ha, fill=method)) +

           geom_bar(position = "dodge", stat="identity") + 

           geom_hline(aes(yintercept = gFDR), colour = 'black', lwd=1.5, lty=2) + 


            xlab(expression(rho[E])) + ylab("FDR") + theme_bw() +
               scale_y_continuous(expand = c(0, 0), limits =  c(0, 0.52)) +

            scale_fill_manual(values=c("#1F78B4", "#FF7F00"))+
            labs(fill='Method')  + theme_classic() + 
             theme(axis.title.x = element_text(size=20, face = "bold"),
                    axis.text.x = element_text(size=20, face = "bold"),
                    axis.title.y = element_text(size=20, face = "bold"),
                    axis.text.y = element_text(size=20,face = "bold"),
                    legend.position = "none",
                     legend.title = element_text(size=20,face = "bold"),
                     legend.text  = element_text(size=20,face = "bold")) 


        Ha = plot_grid(fig_PWR, fig_FDR, fig_AUC, nrow=1)
        ggsave(paste0("Ha", "g", g*100, "rhoX", rhoX*10, fig_type), Ha, width=64, height=20, units='cm')

       #########################################################################################################
       #########################################################################################################
       #########################################################################################################
       #########################################################################################################
       
        rhoE <- method <- character(0)
        AUC <- numeric(0)
        for (k in 1:n_rho) 
        {
            path = paste("h", h, "g", g, "rhoE", rho[k], "rhoX", rhoX, sep="_")
            res <- readRDS(paste0(path, "_.rds"))
            R <- nrow(res$AUC)
            method <- c(method, c("VIMCO", "BVSR", "sLMM", "mvLMM"))
            rhoE <- c(rhoE, rep(rho[k], 4))
            AUC <- c(AUC, c(apply(res$AUC_Hb[, c(1,2,4,3)], 2, "mean")))
        }
        result =  data.frame(AUC = AUC, 
                             rhoE    = factor(rhoE, levels=as.character(rho)),
                             method = factor(method, levels=c("VIMCO", "BVSR","sLMM", "mvLMM")), 
                             stringsAsFactors = FALSE)



        fig_AUC <- ggplot(result, aes(x=rhoE, y=AUC, fill=method)) +

           geom_bar(position = "dodge", stat="identity") + 
            xlab(expression(rho[E])) + ylab("AUC") + theme_bw() +
              coord_cartesian(ylim=c(.7, 1)) + 

            scale_fill_manual(values=c("#A6CEE3", "#FDBF6F", "#CAB2D6", "#B2DF8A"))+
            labs(fill='Method')  + theme_classic() + 
             theme(axis.title.x = element_text(size=20, face = "bold"),
                    axis.text.x = element_text(size=20, face = "bold"),
                    axis.title.y = element_text(size=20, face = "bold"),
                    axis.text.y = element_text(size=20,face = "bold"),
                    legend.position = c(0.88, 0.9),
                     legend.title = element_text(size=20,face = "bold"),
                     legend.text  = element_text(size=20,face = "bold"))
     
        method <- rhoE <- character(0)
        PWR_Hb <- numeric(0)
        FDR_Hb <- numeric(0)
        for (k in 1:n_rho) 
        {
            path = paste("h", h, "g", g, "rhoE", rho[k], "rhoX", rhoX, sep="_")
            res <- readRDS(paste0(path, "_.rds"))

            method <- c(method, "VIMCO", "BVSR")
            rhoE <- c(rhoE, rep(rho[k], 2))
            PWR_Hb <- c(PWR_Hb, res$PWR_Hb[1:2])
            FDR_Hb <- c(FDR_Hb, res$FDR_Hb[1:2])
        }

        # PWR
        result <- data.frame(method = factor(method, levels = c("VIMCO", "BVSR")),
                                    rhoE   = factor(rhoE,   levels=as.character(rho)),
                                       PWR_Hb = PWR_Hb,
                                        stringsAsFactors = FALSE)
        

        fig_PWR <- ggplot(result, aes(x=result$rho, y=PWR_Hb, fill=method)) +

             geom_bar(position = "dodge", stat="identity") + 

            xlab(expression(rho[E])) + ylab("Power") +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 0.832)) +
            scale_fill_manual(values=c("#A6CEE3", "#FDBF6F"))+
            labs(fill='Method')  + theme_classic() + 
               theme(axis.title.x = element_text(size=20, face = "bold"),
                    axis.text.x = element_text(size=20, face = "bold"),
                    axis.title.y = element_text(size=20, face = "bold"),
                    axis.text.y = element_text(size=20,face = "bold"),
                    panel.border = element_blank(),
                    legend.position = "none",
                     legend.title = element_text(size=20,face = "bold"),
                     legend.text  = element_text(size=20,face = "bold")) 



        # FDR
        result_FDR <- data.frame(method = factor(method, levels = c("VIMCO", "BVSR")),
                                 rhoE   = factor(rhoE,   levels=as.character(rho)),
                                 FDR_Hb = FDR_Hb,
                                 stringsAsFactors = FALSE)
    

        fig_FDR <- ggplot(result_FDR, aes(x=result$rho, y=FDR_Hb, fill=method)) +

           geom_bar(position = "dodge", stat="identity") + 

           geom_hline(aes(yintercept = gFDR), colour = 'black', lwd=1.5, lty=2) + 


            xlab(expression(rho[E])) + ylab("FDR") + theme_bw() +
               scale_y_continuous(expand = c(0, 0), limits =  c(0, 0.52)) +

            scale_fill_manual(values=c("#A6CEE3", "#FDBF6F"))+
            labs(fill='Method')  + theme_classic() + 
             theme(axis.title.x = element_text(size=20, face = "bold"),
                    axis.text.x = element_text(size=20, face = "bold"),
                    axis.title.y = element_text(size=20, face = "bold"),
                    axis.text.y = element_text(size=20,face = "bold"),
                    legend.position = "none",
                     legend.title = element_text(size=20,face = "bold"),
                     legend.text  = element_text(size=20,face = "bold")) 


        Hb = plot_grid(fig_PWR, fig_FDR, fig_AUC, nrow=1)
        ggsave(paste0("Hb", "g", g*100, "rhoX", rhoX*10, fig_type), Hb, width=64, height=20, units='cm')
      


       
    }
}


