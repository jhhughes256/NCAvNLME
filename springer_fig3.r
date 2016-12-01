###3dplot_final3.r
##Goal: Create facet plots that represent data from NCAvNLME study
#Changes: No all scenario plots, more looking into how different limited scenarios can show different data

# Remove any previous objects in the workspace
   rm(list=ls(all=TRUE))
   graphics.off()
   
# Load libraries
   library(ggplot2)
   library(grid)
   library(plyr)
   
# Source functions file and set directory
   master.dir <- "E:/hscpw-df1/Data1/Jim Hughes/Interface"      ### Directory containing source files
   setwd(paste(master.dir,"3D Graphs",sep="/"))
   source("functions_NCAvNLME_2016.r")
   work.dir <- paste(master.dir,"Result History (2016)",sep="/")
   numscen <- 15												#Number of scenarios per run (number of runs should be an even number, 9 scen in first run + 6 in second)

# Read in results
   filelist <- list.files(work.dir)								#Record file names
   filenum <- length(filelist)									#Record number of files
   balfile <- filelist[(1:(filenum/2))*2]						#Isolate balfazar and ANOVA file names
   aovfile <- filelist[(1:(filenum/2))*2-1]
   
   baldf <- data.frame(matrix(NA, nrow = filenum*numscen/2, ncol = 80))		#Create empty df for results to be placed in
   colnames(baldf) <- c("FileList","RUVprop","RUVadd","X.RUV.LLOQ","LLOQ","X.BLOQ","F1pop","F1bsv",
					    "IPRED_PBIOQ","NCA_PBIOQ","M1F1_PBIOQ","M1PH_PBIOQ","M3F1_PBIOQ","M3PH_PBIOQ",
					    "IPRED_PCRAT","NCA_PCRAT","M1_PCRAT","M3_PCRAT","NCA_FT1","M1F1_FT1","M3F1_FT1",
					    "M1PH_FT1","M3PH_FT1","NCA_FT2","M1F1_FT2","M3F1_FT2","M1PH_FT2","M3PH_FT2",
					    "NCA_CT1","M1_CT1","M3_CT1","NCA_CT2","M1_CT2","M3_CT2","F_NCA_Sens","F_M1F1_Sens",
					    "F_M3F1_Sens","F_M1PH_Sens","F_M3PH_Sens","F_NCA_Spec","F_M1F1_Spec","F_M3F1_Spec",
					    "F_M1PH_Spec","F_M3PH_Spec","F_NCA_Acc","F_M1F1_Acc","F_M3F1_Acc","F_M1PH_Acc",
					    "F_M3PH_Acc","F_NCA_PPV","F_M1F1_PPV","F_M3F1_PPV","M1PH_PPV","F_M3PH_PPV",
					    "F_NCA_NPV","F_M1F1_NPV","F_M3F1_NPV","F_M1PH_NPV","F_M3PH_NPV","C_NCA_Sens",
					    "C_M1_Sens","C_M3_Sens","C_NCA_Spec","C_M1_Spec","C_M3_Spec","C_NCA_Acc",
					    "C_M1_Acc","C_M3_Acc","C_NCA_PPV","C_M1_PPV","C_M3_PPV","C_NCA_NPV","C_M1_NPV",
					    "C_M3_NPV","X.M1msuc","X.M3msuc","X.M1csuc","X.M3csuc","m1nsim","m3nsim")
   aovdf <- baldf
   
   for(i in 1:(filenum/2))										#Copying across values to placeholder df
   {
      tempbal <- read.csv(paste(work.dir,balfile[i],sep="/"))
	  baldf[(1:15)+15*(i-1),] <- tempbal[-1]
	  tempaov <- read.csv(paste(work.dir,aovfile[i],sep="/"))
	  aovdf[(1:15)+15*(i-1),] <- tempaov[-1]
   }

#Creating separate df sets to be used to highlight effects of BOV   
   bovset <- rep(c("noBOV","BOVonCL","BOV"),each=numscen)
   baldf.bov <- baldf[c(76:105,1:15),]
   baldf.bov <- cbind(baldf.bov,bovset)
   aovdf.bov <- aovdf[c(76:105,1:15),]
   aovdf.bov <- cbind(aovdf.bov,bovset)

   df01 <- data.frame(rep(baldf$FileList,times=5),rep(c("NCA","M1Frel","M1PH","M3Frel","M3PH"),each=filenum*numscen/2),
 	 			      rep(baldf$X.BLOQ,times=5),rep(baldf$RUVadd,times=5),c(baldf$F_NCA_Acc,baldf$F_M1F1_Acc,baldf$F_M1PH_Acc,baldf$F_M3F1_Acc,baldf$F_M3PH_Acc),
					  c(baldf$F_NCA_Spec,baldf$F_M1F1_Spec,baldf$F_M1PH_Spec,baldf$F_M3F1_Spec,baldf$F_M3PH_Spec),
					  c(baldf$F_NCA_Sens,baldf$F_M1F1_Sens,baldf$F_M1PH_Sens,baldf$F_M3F1_Sens,baldf$F_M3PH_Sens),
					  c(rep("BOV",times=numscen*(filenum-4)/2),rep("noBOV",times=numscen),rep("BOVonCL",times=numscen)))
   colnames(df01) <- c("Run Name","Method","BLOQ","RUVadd","Accuracy","Specificity","Sensitivity","BOV") #All BALTH scenarios df 
  
   df02 <- df01[c(1:30,(1:30)+(filenum*numscen/2),(1:30)+(filenum*numscen),(1:30)+(filenum*numscen*3/2),(1:30)+(filenum*numscen*2)),] #Limited BALTH scenarios df 
   
   df03 <- data.frame(rep(baldf.bov$FileList,times=5),rep(c("NCA","M1Frel","M1PH","M3Frel","M3PH"),each=numscen*3),
 	 			      rep(baldf.bov$X.BLOQ,times=5),rep(baldf.bov$RUVadd,times=5),c(baldf.bov$F_NCA_Acc,baldf.bov$F_M1F1_Acc,baldf.bov$F_M1PH_Acc,baldf.bov$F_M3F1_Acc,baldf.bov$F_M3PH_Acc),
					  c(baldf.bov$F_NCA_Spec,baldf.bov$F_M1F1_Spec,baldf.bov$F_M1PH_Spec,baldf.bov$F_M3F1_Spec,baldf.bov$F_M3PH_Spec),
					  c(baldf.bov$F_NCA_Sens,baldf.bov$F_M1F1_Sens,baldf.bov$F_M1PH_Sens,baldf.bov$F_M3F1_Sens,baldf.bov$F_M3PH_Sens),
					  rep(baldf.bov$bovset,times=5))
   colnames(df03) <- colnames(df01) #BOV BALTH scenarios df 
   
   df04 <- data.frame(rep(aovdf$FileList,times=5),rep(c("NCA","M1Frel","M1PH","M3Frel","M3PH"),each=filenum*numscen/2),
 	 			      rep(aovdf$X.BLOQ,times=5),rep(aovdf$RUVadd,times=5),c(aovdf$F_NCA_Acc,aovdf$F_M1F1_Acc,aovdf$F_M1PH_Acc,aovdf$F_M3F1_Acc,aovdf$F_M3PH_Acc),
					  c(aovdf$F_NCA_Spec,aovdf$F_M1F1_Spec,aovdf$F_M1PH_Spec,aovdf$F_M3F1_Spec,aovdf$F_M3PH_Spec),
					  c(aovdf$F_NCA_Sens,aovdf$F_M1F1_Sens,aovdf$F_M1PH_Sens,aovdf$F_M3F1_Sens,aovdf$F_M3PH_Sens),
					  c(rep("BOV",times=numscen*(filenum-4)/2),rep("noBOV",times=numscen),rep("BOVonCL",times=numscen)))
   colnames(df04) <- colnames(df01) #All ANOVA scenarios df 
   
   df05 <- df01[c(1:30,(1:30)+(filenum*numscen/2),(1:30)+(filenum*numscen),(1:30)+(filenum*numscen*3/2),(1:30)+(filenum*numscen*2)),] #Limited ANOVA scenarios df 
   
   df06 <- data.frame(rep(aovdf.bov$FileList,times=5),rep(c("NCA","M1Frel","M1PH","M3Frel","M3PH"),each=numscen*3),
 	 			      rep(aovdf.bov$X.BLOQ,times=5),rep(aovdf.bov$RUVadd,times=5),c(aovdf.bov$F_NCA_Acc,aovdf.bov$F_M1F1_Acc,aovdf.bov$F_M1PH_Acc,aovdf.bov$F_M3F1_Acc,aovdf.bov$F_M3PH_Acc),
					  c(aovdf.bov$F_NCA_Spec,aovdf.bov$F_M1F1_Spec,aovdf.bov$F_M1PH_Spec,aovdf.bov$F_M3F1_Spec,aovdf.bov$F_M3PH_Spec),
					  c(aovdf.bov$F_NCA_Sens,aovdf.bov$F_M1F1_Sens,aovdf.bov$F_M1PH_Sens,aovdf.bov$F_M3F1_Sens,aovdf.bov$F_M3PH_Sens),
					  rep(aovdf.bov$bovset,times=5)) 
   colnames(df06) <- colnames(df01) #BOV ANOVA df
   
   ## New set of limited scenario df
   # Odd numbers BALTH, even numbers ANOVA
   # Choosing Sim 1, 4, 5 and half of Sim 3
   
   df07 <- df01[c(c(1:9,46:54),c(1:9,46:54)+(filenum*numscen/2),c(1:9,46:54)+(filenum*numscen),c(1:9,46:54)+(filenum*numscen*3/2),c(1:9,46:54)+(filenum*numscen*2)),]
   df07$RUVaddf <- as.factor(df07$RUVadd)
   levels(df07$RUVaddf) <- c("0.0005","0.001","0.0015")
   
   df08 <- df04[c(c(1:9,46:54),c(1:9,46:54)+(filenum*numscen/2),c(1:9,46:54)+(filenum*numscen),c(1:9,46:54)+(filenum*numscen*3/2),c(1:9,46:54)+(filenum*numscen*2)),]
   df08$RUVaddf <- as.factor(df08$RUVadd)
   levels(df08$RUVaddf) <- c("0.0005","0.001","0.0015")
   
   # Try adding more plots to BOV limited scenarios
   # Choosing Sim 1, 4
   
   bovset2 <- rep(c("noBOV","BOVonCL","BOV","BOV"),each=numscen)
   baldf.bov2 <- baldf[c(76:105,1:15,46:60),]
   baldf.bov2 <- cbind(baldf.bov2,bovset2)
   aovdf.bov2 <- aovdf[c(76:105,1:15,46:60),]
   aovdf.bov2 <- cbind(aovdf.bov2,bovset2)
   
   df09 <- data.frame(rep(baldf.bov2$FileList,times=5),rep(c("NCA","M1Frel","M1PH","M3Frel","M3PH"),each=numscen*4),
 	 			      rep(baldf.bov2$X.BLOQ,times=5),rep(baldf.bov2$RUVadd,times=5),c(baldf.bov2$F_NCA_Acc,baldf.bov2$F_M1F1_Acc,baldf.bov2$F_M1PH_Acc,baldf.bov2$F_M3F1_Acc,baldf.bov2$F_M3PH_Acc),
					  c(baldf.bov2$F_NCA_Spec,baldf.bov2$F_M1F1_Spec,baldf.bov2$F_M1PH_Spec,baldf.bov2$F_M3F1_Spec,baldf.bov2$F_M3PH_Spec),
					  c(baldf.bov2$F_NCA_Sens,baldf.bov2$F_M1F1_Sens,baldf.bov2$F_M1PH_Sens,baldf.bov2$F_M3F1_Sens,baldf.bov2$F_M3PH_Sens),
					  rep(baldf.bov2$bovset,times=5))
   colnames(df09) <- colnames(df01) #BOV BALTH df 
   
   df10 <- data.frame(rep(aovdf.bov2$FileList,times=5),rep(c("NCA","M1Frel","M1PH","M3Frel","M3PH"),each=numscen*4),
 	 			      rep(aovdf.bov2$X.BLOQ,times=5),rep(aovdf.bov2$RUVadd,times=5),c(aovdf.bov2$F_NCA_Acc,aovdf.bov2$F_M1F1_Acc,aovdf.bov2$F_M1PH_Acc,aovdf.bov2$F_M3F1_Acc,aovdf.bov2$F_M3PH_Acc),
					  c(aovdf.bov2$F_NCA_Spec,aovdf.bov2$F_M1F1_Spec,aovdf.bov2$F_M1PH_Spec,aovdf.bov2$F_M3F1_Spec,aovdf.bov2$F_M3PH_Spec),
					  c(aovdf.bov2$F_NCA_Sens,aovdf.bov2$F_M1F1_Sens,aovdf.bov2$F_M1PH_Sens,aovdf.bov2$F_M3F1_Sens,aovdf.bov2$F_M3PH_Sens),
					  rep(aovdf.bov2$bovset,times=5)) 
   colnames(df10) <- colnames(df01) #BOV ANOVA df
   
   ## New New set of limited scenario df
   # Odd numbers BALTH, even numbers ANOVA
   # Choosing half of each sim
   vec01 <- c(1:9,1:9+45,1:9+60)
   df11 <- df01[c(vec01,vec01+(filenum*numscen/2),vec01+(filenum*numscen),vec01+(filenum*numscen*3/2),vec01+(filenum*numscen*2)),]
   df11$RUVaddf <- as.factor(df11$RUVadd)
   levels(df11$RUVaddf) <- c("0.0005","0.001","0.0015")
   df12 <- df04[c(vec01,vec01+(filenum*numscen/2),vec01+(filenum*numscen),vec01+(filenum*numscen*3/2),vec01+(filenum*numscen*2)),]
   df12$RUVaddf <- as.factor(df12$RUVadd)
   levels(df12$RUVaddf) <- c("0.0005","0.001","0.0015")
   
   vec02 <- c(1:15,1:15+30,1:15+45,1:15+60)
   df13 <- df01[c(vec02,vec02+(filenum*numscen/2),vec02+(filenum*numscen),vec02+(filenum*numscen*3/2),vec02+(filenum*numscen*2)),]
   df14 <- df04[c(vec02,vec02+(filenum*numscen/2),vec02+(filenum*numscen),vec02+(filenum*numscen*3/2),vec02+(filenum*numscen*2)),]

   ## DF's in use  3  6  7  8
     
	 theme_set(theme_bw()) 
   
   df07$RUVadd <- "RUV"
#### RUV LOQ RES
  ##		 HALF SIMS 1,4
   plotobj01 <- NULL
   plotobj01 <- ggplot(df07, aes(x=BLOQ,y=Accuracy*100,colour=Method,shape=Method)) 
   plotobj01 <- plotobj01 + geom_line(size=1)
   plotobj01 <- plotobj01 + geom_point(size=2,fill="white")
   plotobj01 <- plotobj01 + scale_shape_manual(values = c(21, 22, 23, 24, 25))
   plotobj01 <- plotobj01 + scale_y_continuous("Percent Accuracy\n",limits=c(30,100))
   plotobj01 <- plotobj01 + scale_x_continuous("\nPercent Below LOQ")
   plotobj01 <- plotobj01 + facet_wrap(~RUVadd+RUVaddf)
   ggsave("NCA_Fig3.eps")
   ggsave("NCA_Fig3.pdf")
  
   plotobj02 <- NULL
   plotobj02 <- ggplot(df07, aes(x=BLOQ,y=Specificity*100,colour=Method,shape=Method)) 
   plotobj02 <- plotobj02 + geom_line(size=1)
   plotobj02 <- plotobj02 + geom_point(size=2,fill="white")
   plotobj02 <- plotobj02 + scale_shape_manual(values = c(21, 22, 23, 24, 25))
   plotobj02 <- plotobj02 + scale_y_continuous("Percent Specificity\n",limits=c(30,100))
   plotobj02 <- plotobj02 + scale_x_continuous("\nPercent Below LOQ")
   plotobj02 <- plotobj02 + facet_wrap(~RUVadd+RUVaddf)
   ggsave("NCA_Fig5.eps")
   ggsave("NCA_Fig5.pdf")

   plotobj03 <- NULL
   plotobj03 <- ggplot(df07, aes(x=BLOQ,y=Sensitivity*100,colour=Method,shape=Method))
   plotobj03 <- plotobj03 + geom_line(size=1)
   plotobj03 <- plotobj03 + geom_point(size=2,fill="white")
   plotobj03 <- plotobj03 + scale_shape_manual(values = c(21, 22, 23, 24, 25))
   plotobj03 <- plotobj03 + scale_y_continuous("Percent Sensitivity\n",limits=c(30,100))
   plotobj03 <- plotobj03 + scale_x_continuous("\nPercent Below LOQ")
   plotobj03 <- plotobj03 + facet_wrap(~RUVadd+RUVaddf)
   ggsave("NCA_Fig4.eps")
   ggsave("NCA_Fig4.pdf")
   
#### BOV LOQ RES
  ##		 
  ##		 
   
   plotobj04 <- NULL
   plotobj04 <- ggplot(df03, aes(x=BLOQ,y=Accuracy*100,colour=Method))
   plotobj04 <- plotobj04 + geom_line(size=1)
   plotobj04 <- plotobj04 + scale_y_continuous("Percent Accuracy\n",limits=c(30,100))
   plotobj04 <- plotobj04 + scale_x_continuous("\nPercent Below LOQ")
   plotobj04 <- plotobj04 + facet_wrap(~BOV)
   ggsave("NCA_Fig6.eps")
   ggsave("NCA_Fig6.pdf")