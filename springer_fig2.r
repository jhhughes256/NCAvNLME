### NONMEM vs. NCA in Bioequivalence Studies
  # All values needing definition occur before simulation code

# Remove any previous objects in the workspace
   rm(list=ls(all=TRUE))
   graphics.off()

# Define dose, blq and type of estimation code to be used
   blq <- 0.01
   ruv.blq <- 0.2
   amt <- 125
   #dif.ka <- 0       # Different ka for formulations;   0 <- OFF   1 <- ON   WIP

# Source functions file and NONMEM .ctl reference txt
   master.dir <- "E:/hscpw-df1/Data1/Jim Hughes/2016/01_MDV"      ### Directory containing source files
   setwd(master.dir)
   source("functions_NCAvNLME_2016.r")
   ctlm1 <- readLines("NONMEM_ref_M1_noqv3bsv+bov16.ctl")          # Suffix to file determines how ambitious your model will be
   ctlm3 <- readLines("NONMEM_ref_M3_noqv3bsv+bov16.ctl")          # noqbsv+bov&noclbov          (example)
   wfn.dir <- "c:/nm72/wfn7/bin/wfn.bat"

# Clear and set working directory
   SIM.name.out <- "Run1_Scen1"            ### Change for each new study design to avoid rewriting data
   SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
   SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
   EST.dir <- paste(SIM.dir,"ctl",sep="/")
   EST.file <- paste(EST.dir,SIM.name.out,sep="/")
   FIT.dir <- paste(SIM.dir,"fit",sep="/")

# Load libraries
   library(ggplot2)
   library(doBy)
   library(plyr)
   library(grid)
   library(MASS)
   library(MBESS)
   library(Cairo)
   set.seed(1234)

# Limit sampling for NCA and NONMEM analysis
   sslimit <- c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96)
   clastlim <- tail(sslimit,1)

# Set up rich sampling times
   T.1 <- seq(from=0, to=3, by=0.05)
   T.25 <- seq(from=3.25, to=6, by=0.25)
   T.5 <- seq(from=6.5, to=12, by=0.5)
   T1 <- seq(from=13, to=clastlim, by=1)
   TIME <- c(T.1,T.25,T.5,T1)

# Limits for bioequivalence
   limitlo <- 0.8
   limithi <- 1.25

# State number of cores to be used while using the server for NONMEM
   ncore.m1 <- 5
   ncore.m3 <- 20

# Define study design
   nid  <- 24            # Number of subjects per study
   nsim <- 500           # Number of studies; must be a multiple of ncore
   nsub <- nsim*nid      # Number of subjects
   nobs <- length(TIME)  # Number of observations per subject

# Simulation Parameters
# Define population parameter (THETA) values
   CLpop   <- 20
   V2pop   <- 100
   Qpop    <- 35
   V3pop   <- 400
   KApop   <- 0.5
   F1pop   <- 1

# Define population variability (OMEGA) values
# Written in variances, converted to std dev in ETA section
# Between subject
   CLbsv   <- 0.045
   V2bsv   <- 0.045
   Qbsv    <- 0.045
   V3bsv   <- 0.045
   KAbsv   <- 0.01
   F1bsv   <- 0.1225
# Between occasion
   CLbov   <- 0.045
   V2bov   <- 0.045
   Qbov    <- 0.045
   V3bov   <- 0.045

# Define random unexplained variability (SIGMA) values; can't both equal zero
# RUV norm <- std dev
# RUV.nm   <- variance
   RUVprop <- 0.05              # 5 percent proportional
   RUVadd  <- 0.0015   # Ensures 20 percent proportional at lloq
   ifelse(blq==0 || RUVadd==0,
     RUVadd.nm <- paste(RUVadd,"FIX",sep=" "),
     RUVadd.nm <- RUVadd)
   trunc.blq <- 0.01

# Correlation between parameters
#               CL  V2  V3   Q
   CORvec <- c(1,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,1)

   CORmat <- matrix(CORvec,nrow=4,ncol=4)

# Convert BSV to SD and combine with correlation matrix to form covariance matrix
   STD <- c(sqrt(CLbsv),sqrt(V2bsv),sqrt(V3bsv),sqrt(Qbsv))
   OMEGA <- cor2cov(CORmat,STD)

   parameters <- read.csv(paste(SIM.file,"PARAMETERS.csv",sep="_"))
   #simdata <- read.csv(paste(SIM.file,"_RAW.csv", sep=""))

   ipredfresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_IPREDFresult.csv",sep=""))
   ipredcresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_IPREDCresult.csv",sep=""))

   ncafresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_NCAFresult.csv",sep=""))
   ncacresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_NCACresult.csv",sep=""))

   trunc.file <- paste(SIM.name.out,"TRUNCATED.csv",sep="_")
   limdata <- read.csv(paste(SIM.dir,trunc.file,sep="/"))
   limobs <- length(sslimit)
   per.bloq <- percent.blq(limdata$DV,limdata$TIME,trunc.blq)

   m1sim.in <- read.csv(paste(SIM.file,"M1_NMTHETAS.csv", sep="_"))
   m3sim.in <- read.csv(paste(SIM.file,"M3_NMTHETAS.csv", sep="_"))
   m1fresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M1Fresult.csv",sep=""))
   m1cresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M1Cresult.csv",sep=""))
   m1phfresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M1PHFresult.csv",sep=""))
   m3fresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3Fresult.csv",sep=""))
   m3cresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3Cresult.csv",sep=""))
   m3phfresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3PHFresult.csv",sep=""))

   conf.int <- read.csv(paste(SIM.file,"CONF_INTtable.csv",sep="_"))
   bioq.table <- read.csv(paste(SIM.file,"BIOQtable.csv",sep="_"))
   final.table <- read.csv(paste(SIM.file,"FINALtable.csv",sep="_"))
   aovbioqtable <- read.csv(paste(SIM.file,"AOVbioq.csv",sep="_"))

   ipredftable <- glohipercent.func(ipredfresult,"Simulation","FREL")
   ncaftable <- glohipercent.func(ncafresult,"Non-Compartmental Analysis","FREL")
   m1ftable <- glohipercent.func(m1fresult,"Non-Linear Mixed Effects M1 Frel","FREL")
   m1phfresult.in <- m1phfresult
   colnames(m1phfresult.in) <- c("SIM_ID","FREL_GMEAN","FREL_GLO95","FREL_GHI95")
   m1phftable <- glohipercent.func(m1phfresult.in,"Non-Linear Mixed Effects M1 Post-Hoc","FREL")
   m3ftable <- glohipercent.func(m3fresult,"Non-Linear Mixed Effects M3 Frel","FREL")
   m3phfresult.in <- m3phfresult
   colnames(m3phfresult.in) <- c("SIM_ID","FREL_GMEAN","FREL_GLO95","FREL_GHI95")
   m3phftable <- glohipercent.func(m3phfresult.in,"Non-Linear Mixed Effects M3 Post-Hoc","FREL")
   allftable <- rbind(ipredftable,ncaftable,m1ftable,m1phftable,m3ftable,m3phftable)

   m1.nsim <- final.table[1, 35]
   m3.nsim <- final.table[1, 36]

   ipredbioq <- bioq.table[bioq.table$Data == "IPRED", ]
   ncabioq <- bioq.table[bioq.table$Data == "NCA", ]
   m1f1bioq <- bioq.table[bioq.table$Data == "M1F1", ]
   m3f1bioq <- bioq.table[bioq.table$Data == "M3F1", ]
   m1phbioq <- bioq.table[bioq.table$Data == "M1PH", ]
   m3phbioq <- bioq.table[bioq.table$Data == "M3PH", ]

   theme_set(theme_bw())

   outputdf3 <- allftable[1:5]
   outputdf3$ORDER <- outputdf3$FREL_GHI95-outputdf3$FREL_GLO95
   outputdf3 <- orderBy(~Method+ORDER,outputdf3)
   outputdf3$SIMorder <- c(rep(1:nsim,times=2),rep(1:m1.nsim,times=2),rep(1:m3.nsim,times=2))
   outputdf3$LORED <- ifelse(outputdf3$FREL_GLO95<0.8|outputdf3$FREL_GLO95>1.25,outputdf3$FREL_GLO95,NA)
   outputdf3$HIRED <- ifelse(outputdf3$FREL_GHI95<0.8|outputdf3$FREL_GHI95>1.25,outputdf3$FREL_GHI95,NA)
   outputdf3$LOBLU <- ifelse(outputdf3$FREL_GLO95>=0.8&outputdf3$FREL_GLO95<=1.25,outputdf3$FREL_GLO95,NA)
   outputdf3$HIBLU <- ifelse(outputdf3$FREL_GHI95>=0.8&outputdf3$FREL_GHI95<=1.25,outputdf3$FREL_GHI95,NA)

   annote.temp <- c(paste(signif(mean(ipredbioq$p)*100,3),"%",sep=""),paste(signif(mean(ncabioq$p)*100,3),"%",sep=""),paste(signif(mean(m1f1bioq$p)*100,3),"%",sep=""),
			        paste(signif(mean(m1phbioq$p)*100,3),"%",sep=""),paste(signif(mean(m3f1bioq$p)*100,3),"%",sep=""),paste(signif(mean(m3phbioq$p)*100,3),"%",sep=""))
   annote.text1	<- data.frame(c("Simulation","Non-Compartmental Analysis","Non-Linear Mixed Effects M1 Frel","Non-Linear Mixed Effects M1 Post-Hoc",
                                "Non-Linear Mixed Effects M3 Frel","Non-Linear Mixed Effects M3 Post-Hoc"),annote.temp)
   colnames(annote.text1) <- c("Method","Bioq")

   plotobj3 <- NULL
   plotobj3 <- ggplot(outputdf3)
   plotobj3 <- plotobj3 + geom_point(aes(x=SIMorder, y=HIBLU*100), size=2, colour="blue",alpha=0.5, shape=20)
   plotobj3 <- plotobj3 + geom_point(aes(x=SIMorder, y=LOBLU*100), size=2, colour="blue",alpha=0.5, shape=20)
   plotobj3 <- plotobj3 + geom_point(aes(x=SIMorder, y=HIRED*100), size=2, colour="red",alpha=0.5, shape=18)
   plotobj3 <- plotobj3 + geom_point(aes(x=SIMorder, y=LORED*100), size=2, colour="red",alpha=0.5, shape=18)
   plotobj3 <- plotobj3 + geom_text(data=annote.text1, aes(x = (0.1*nsim), y = 140, label=Bioq), size = 4, colour="black")
   plotobj3 <- plotobj3 + geom_hline(yintercept=limitlo*100, colour="black", linetype="dashed")
   plotobj3 <- plotobj3 + geom_hline(yintercept=limithi*100, colour="black", linetype="dashed")
   plotobj3 <- plotobj3 + scale_x_continuous("Study ordered by CI width",lim=c(0,nsim))
   plotobj3 <- plotobj3 + scale_y_continuous("Lower and Upper CI", lim=c(65,150))
   plotobj3 <- plotobj3 + facet_wrap(~Method,ncol=2)
   ggsave("NCA_Fig2b.pdf")

   #outputdf3$Method <- factor(outputdf3$Method, levels(outputdf3$Method)[c(1,3,4,2,5,6)])
   plotobj4 <- NULL
   plotobj4 <- plotobj3 + facet_wrap(~Method,ncol=3)
