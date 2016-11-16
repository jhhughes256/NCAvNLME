### NONMEM vs. NCA in Bioequivalence Studies
  # All values needing definition occur before simulation code

# Remove any previous objects in the workspace
   rm(list=ls(all=TRUE))
   graphics.off()

# Define dose, blq and type of estimation code to be used
   blq <- SUB01
   ruv.blq <- SUB02 
   amt <- SUB03
   #dif.ka <- 0       # Different ka for formulations;   0 <- OFF   1 <- ON   WIP
   
# Source functions file and NONMEM .ctl reference txt
   master.dir <- "SUB04"      ### Directory containing source files
   setwd(master.dir)
   source("SUB29")
   ctlm1 <- readLines("SUB30")          # Suffix to file determines how ambitious your model will be
   ctlm3 <- readLines("SUB31")          # noqbsv+bov&noclbov          (example)
   wfn.dir <- "c:/nm72/wfn7/bin/wfn.bat"
   
# Clear and set working directory
   SIM.name.out <- "SUB05"            ### Change for each new study design to avoid rewriting data
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
   set.seed(1234)

# Limit sampling for NCA and NONMEM analysis
   sslimit <- SUB36
   clastlim <- tail(sslimit,1)
   
# Set up rich sampling times
   T.1 <- SUB37
   T.25 <- SUB38
   T.5 <- SUB39
   T1 <- SUB40
   TIME <- c(T.1,T.25,T.5,T1)
  
# Limits for bioequivalence
   limitlo <- SUB34
   limithi <- SUB35
   
# State number of cores to be used while using the server for NONMEM
   ncore.m1 <- SUB07
   ncore.m3 <- SUB08

# Define study design
   nid  <- SUB33            # Number of subjects per study
   nsim <- SUB09           # Number of studies; must be a multiple of ncore
   nsub <- nsim*nid      # Number of subjects 
   nobs <- length(TIME)  # Number of observations per subject

# Simulation Parameters 
# Define population parameter (THETA) values
   CLpop   <- SUB10
   V2pop   <- SUB11
   Qpop    <- SUB12
   V3pop   <- SUB13
   KApop   <- SUB14
   F1pop   <- SUB15

# Define population variability (OMEGA) values
# Written in variances, converted to std dev in ETA section
# Between subject
   CLbsv   <- SUB16
   V2bsv   <- SUB17
   Qbsv    <- SUB18
   V3bsv   <- SUB19
   KAbsv   <- SUB20
   F1bsv   <- SUB21
# Between occasion
   CLbov   <- SUB22
   V2bov   <- SUB23
   Qbov    <- SUB24
   V3bov   <- SUB25

# Define random unexplained variability (SIGMA) values; can't both equal zero
# RUV norm <- std dev
# RUV.nm   <- variance
   RUVprop <- SUB26              # 5 percent proportional                       
   RUVadd  <- SUB27   # Ensures 20 percent proportional at lloq
   ifelse(blq==0 || RUVadd==0,
     RUVadd.nm <- paste(RUVadd,"FIX",sep=" "),
     RUVadd.nm <- RUVadd)        
   trunc.blq <- SUB28
  
# Correlation between parameters             
#               CL  V2  V3   Q
   CORvec <- SUB32
   
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
   m3fresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3Fresult.csv",sep=""))
   m3cresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3Cresult.csv",sep=""))
   
   conf.int <- read.csv(paste(SIM.file,"CONF_INTtable.csv",sep="_"))
   bioq.table <- read.csv(paste(SIM.file,"BIOQtable.csv",sep="_"))
   final.table <- read.csv(paste(SIM.file,"FINALtable.csv",sep="_"))
   aovbioqtable <- read.csv(paste(SIM.file,"AOVbioq.csv",sep="_"))
