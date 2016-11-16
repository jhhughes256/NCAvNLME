### NONMEM vs. NCA in Bioequivalence Studies
  # All values needing definition occur before simulation code
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Source functions file and NONMEM .ctl reference txt
  master.dir <- "E:/hscpw-df1/Data1/Jim Hughes/DDPLY"      ### Directory containing source files
  setwd(master.dir)
  source("functions_NCAvNLME_2016.r")
  ctlm1 <- readLines("NONMEM_ref_M1.ctl")  #No BSV or BOV on Q & V3
  ctlm3 <- readLines("NONMEM_ref_M3.ctl")
  rscript <- readLines("Rscript.r")
  wfn.dir <- "c:/nm72/wfn7/bin/wfn.bat"

# Load libraries
  library(ggplot2)
  library(doBy)
  library(plyr)
  library(grid)
  library(MASS)
  library(MBESS)
  library(parallel)
  set.seed(1234)

# Specify run characteristics
  #set up values that are variable between runs
  rundf <- data.frame(
    RUN = rep(1:7, each = 15), #rep(number.of.runs, each=number.of.scenarios)
    SCEN = rep(1:15, times = 7), #opposite of above
    RUV.TYPE = rep(c(rep(1, 3),rep(2, 2)),21), #changing calculation of RUV for scenarios 10:15
    RUV.BLQ = rep(c(0.2, 0.15, 0.1, 0.1, 0.5), 21),
    F1.POP = rep(c(rep(1.0, 5), rep(0.9, 5), rep(1.11, 5)), 7), #changing Frel of generic
    F1.BSV = rep(c(rep(0.1225, 5), rep(0.0484, 5), rep(0.0529, 5)), 7), #changing BSV on frel
    BLQ = c(rep(0.01, 15), rep(0.1, 15),rep(0.01, 75)),  #Run 2 - raised LLOQ
    RUV.PROP = c(rep(0.05, 30), rep(0.09, 15), rep(0.05, 60)),  #Run 3 - increased proportional RUV
    SS.TYPE = c(rep(1, 45), rep(2, 15), rep(1, 45)),  #Run 4 - reduced sampling schedule
    KA.TYPE = c(rep(1, 60), rep(2, 15), rep(1, 30)),  #Run 5 - 20% lower generic KA
    BOV.TYPE = c(rep(1, 75), rep(2, 15), rep(3, 15)))  #Run 6 & 7 - BOV testing

  #set up values that are constant between runs
  runvec <- c(
    NID = 24,
    NSIM = 80,
    LIMITLO = 0.8,
    LIMITHI = 1.25,
    NCORE.M1 = 5,
    NCORE.M3 = 20,
    AMT = 125,
    CL.POP = 20,
    V2.POP = 100,
    Q.POP = 35,
    V3.POP = 400,
    KA.POP = 0.5,
    KA.GEN = 0.3,
    CL.BSV = 0.045,
    V2.BSV = 0.045,
    V3.BSV = 0.045,
    Q.BSV = 0.045,
    KA.BSV = 0.01,
    CL.BOV = 0.045,
    V2.BOV = 0.045,
    V3.BOV = 0.045,
    Q.BOV = 0.045)

  #set up time vector for simulation
  timevec <- c(
    seq(0, 3, 0.05),
    seq(3.25, 6, 0.25),
    seq(6.5, 12, 0.5),
    seq(13, 96, 1))
  TIME <- timevec

  #set up correlation vector for simulation
  corvec <- c(
    1, 0.3, 0.3, 0.3,
    0.3, 1, 0.3, 0.3,
    0.3, 0.3, 1, 0.3,
    0.3, 0.3 ,0.3 , 1)

### FIRST HALF -----------------------------------------------------------------
  # First loop includes
  # 1. Simulation of data
  # 2. Initial analysis of simulated data
  # 3. Non-Compartmental Analysis
  # 4. Initial analysis of NCA data
  # 5. Non-Linear Mixed Effects
  # 6. Create .r script for loading

  ddply(rundf[1,], .(RUN, SCEN), function(df, vec, time, cor, rtemplate) {
### 1. Simulation of data
  #Set working directory
    SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
    SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
    SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
    EST.dir <- paste(SIM.dir,"ctl",sep="/")
    EST.file <- paste(EST.dir,SIM.name.out,sep="/")
    FIT.dir <- paste(SIM.dir,"fit",sep="/")

  # Define study design
    nid <- vec["NID"]
    nsim <- vec["NSIM"]
    nsub <- nid*nsim
    nobs <- length(time)
    if(df$SS.TYPE == 1) {
      sstime <- c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96)
    }else{
      sstime <- c(0,0.25,0.5,1,2,4,8,16,36,96)
    }

  # Define random unexplained variability (SIGMA) values
    ruv.prop <- df$RUV.PROP
    ruv.blq <- df$RUV.BLQ
    blq <- df$BLQ
    if(df$RUV.TYPE == 1) {  #Scenarios 1-9
      ruv.add <- (ruv.blq - ruv.prop) * blq
      trunc.blq <- blq
    }
    if(df$RUV.TYPE == 2) {  #Scenarios 10-15
      ruv.add <- (0.2 - ruv.prop) * blq
      trunc.blq <- ruv.add/(ruv.blq - ruv.prop)
    }
    ruv.add.nm <- ifelse(blq == 0 || ruv.add == 0,
      paste(ruv.add, "FIX"),
      ruv.add
    )

  # Set cores
    ncore.m1 <- vec["NCORE.M1"]
    ncore.m3 <- vec["NCORE.M3"]

  # Correlation between parameters
    cormat <- matrix(cor, nrow = 4, ncol = 4)
    std.bsv <- c(sqrt(vec["CL.BSV"]),sqrt(vec["V2.BSV"]),sqrt(vec["Q.BSV"]),sqrt(vec["V3.BSV"]))
    omega <- cor2cov(cormat, std.bsv)

  # Account for random BSV
    etamat <- mvrnorm(n = nsub, mu = c(0, 0, 0, 0), omega)
    ETA1 <- rep(etamat[,1],each=2)  # CLbsv
    ETA2 <- rep(etamat[,2],each=2)  # V2bsv
    ETA3 <- rep(etamat[,3],each=2)  # Qbsv
    ETA4 <- rep(etamat[,4],each=2)  # V3bsv
    if(vec["KA.BSV"]==0){
      ETA5 <- rep(0, times=nsub)  # KA (no bsv)
    }else{
      ETA5 <- rnorm(nsub,mean=0,sd=sqrt(vec["KA.BSV"]))  # KA (bsv)
    }
    if(df$F1.BSV==0){
      ETA6 <- rep(0, times=nsub)  # F1 (no bsv)
    }else{
      ETA6 <- rnorm(nsub,mean=0,sd=sqrt(df$F1.BSV))  # F1 (bsv)
    }
    ETA7 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["CL.BOV"]))  # CLbov
    ETA8 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["V2.BOV"]))  # V2bov
    ETA9 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["Q.BOV"]))   # Qbov
    ETA10 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["V3.BOV"]))  # V3bov
    if(vec["KA.BSV"]==0){
      ETA11 <- rep(0, times=nsub)  # KA (no bsv)
    }else{
      ETA11 <- rnorm(nsub,mean=0,sd=sqrt(vec["KA.BSV"]))  # KA (bsv)
    }

    if(df$KA.TYPE == 1) {
      ka.val <- rep(vec["KA.POP"] * exp(ETA5), each = 2)
    }else{
      ka.val <- as.vector(rbind(vec["KA.POP"] * exp(ETA5), vec["KA.GEN"] * exp(ETA11)))
    }

  # Create/clear directories before run
    dir.setup(SIM.dir)
    dir.setup(EST.dir)
    dir.setup(FIT.dir)
    nm.clear(EST.dir, nsim*2)

  # Define theta table
    thetadf <- data.frame(
      ID = rep(1:nid, times = nsim, each = 2),
      AMT = rep(vec["AMT"], times = nsub * 2),
      FORM = rep(1:2, times = nsub),
      CL = vec["CL.POP"] * exp(ETA1 + ETA7),
      V2 = vec["V2.POP"] * exp(ETA2 + ETA8),
      Q = vec["Q.POP"] * exp(ETA3 + ETA9),
      V3 = vec["V3.POP"] * exp(ETA4 + ETA10),
      KA = ka.val,
      F1 = as.vector(rbind(df$F1.POP * exp(ETA6), rep(1, nsim))))
    write.csv(thetadf, file = paste(SIM.file, "THETAS.csv", sep = "_"), row.names = F)
    thetadf$FORM <- NULL

  # Population simulation
    simdataPRED1 <- simulate.2comp.abs(
      ID = 0,
      AMT = vec["AMT"],
      CL = vec["CL.POP"],
      Q = vec["Q.POP"],
      V2 = vec["V2.POP"],
      V3 = vec["V3.POP"],
      KA = vec["KA.POP"],
      F1 = 1)
    names(simdataPRED1) <- c("TIME", "PRED")
    simdataPRED1$FORM <- rep(1, times = nobs)

    simdataPRED2 <- simulate.2comp.abs(
      ID = 0,
      AMT = vec["AMT"],
      CL = vec["CL.POP"],
      Q = vec["Q.POP"],
      V2 = vec["V2.POP"],
      V3 = vec["V3.POP"],
      KA = vec["KA.POP"],
      F1 = df$F1.POP)
    names(simdataPRED2) <- c("TIME", "PRED")
    simdataPRED2$FORM <- rep(2, times = nobs)

    simdataPRED <- rbind(simdataPRED1, simdataPRED2)

  # Individual simulation
    simdataIPRED <- mdply(thetadf,simulate.2comp.abs)
    names(simdataIPRED)[10] <- "IPRED"
    simdataIPRED$UID  <- rep(1:nsub, each = nobs * 2)
    simdataIPRED$SIM  <- rep(1:nsim, each = nobs * nid * 2)
    simdataIPRED$FORM <- ifelse(simdataIPRED$F1 == 1, 1, 2)

  # Combine PRED and IPRED
    simdata <- orderBy(~UID + FORM + TIME,
      merge(simdataIPRED, simdataPRED, all = T)
    )[, c(12 ,3, 13, 1, 4, 2, 5:10, 14, 11)]

  # Add residual error
    CP <- simdata$IPRED
    tnobs <- length(CP)
    EPS1 <- rnorm(tnobs, mean = 0, sd = ruv.prop)
    EPS2 <- rnorm(tnobs, mean = 0, sd = ruv.add)
    simdata$DV <- CP * (1 + EPS1) + EPS2
    simdata$DV <- ifelse(simdata$TIME == 0, 0, simdata$DV)

  # Save simulation file
    write.csv(simdata, file = paste(SIM.file, "_RAW.csv", sep = ""), row.names = F)

### 2. Initial analysis of simulated data
  # Process results for IPRED (see functions utility)
    data.process(simdata,TIME,nid,nsim,SIM.file,trunc.blq,mode=1)

### 3. Non-Compartmental Analysis
### 4. Initial analysis of NCA data
  # Process results for NCA with truncation (see functions utility)
    data.process(simdata,sstime,nid,nsim,SIM.file,trunc.blq,mode=2)

### 5. Non-Linear Mixed Effects
  # Source truncated dataset from NCA data.process function
    trunc.file <- paste(SIM.name.out,"TRUNCATED.csv",sep="_")
    limdata <- read.csv(paste(SIM.dir,trunc.file,sep="/"))
    limobs <- length(sstime)
    per.bloq <- percent.blq(limdata$DV,limdata$TIME,trunc.blq)

  # Create .ctl for NONMEM to run for each SIM (and a .bat placeholder)
  # Replaces phrases in the .ctl reference file to values stated at the beginning of the script
    input.vector <- c(vec["CL.POP"], vec["V2.POP"], vec["V3.POP"], vec["Q.POP"], df$F1.POP, vec["KA.POP"],
      vec["CL.BSV"], vec["V2.BSV"], vec["V3.BSV"], vec["Q.BSV"], vec["KA.BSV"], df$F1.BSV,
      trunc.blq, vec["CL.BOV"], vec["V2.BOV"], vec["V3.BOV"], vec["Q.BOV"],
      omega[2,1], omega[3,1], omega[3,2], omega[4,1], omega[4,2], omega[4,3],
      ruv.prop, ruv.add.nm)
    phrase.vector <- c("clpop", "v2pop", "v3pop", "qpop", "f1pop", "kapop",
      "clbsv", "v2bsv", "v3bsv", "qbsv", "kabsv", "f1bsv",
      "blq", "clbov", "v2bov", "v3bov", "qbov",
      "mat21", "mat31", "mat32", "mat41", "mat42", "mat43",
      "ruvprop", "ruvadd")
    ctlm1 <- ctl.rep(phrase.vector,input.vector,ctlm1)
    ctlm3 <- ctl.rep(phrase.vector,input.vector,ctlm3)

  # Create .bat files for NONMEM shell and run in NONMEM (number of .bats determined by ncore value)
    nmbat1 <- rep("",times=nsim)
    nmbat1 <- nm.prep(nsim,ncore.m1,nid,EST.file,1,SIM.name.out,limdata,limobs,trunc.blq,ctlm1,nmbat1,vec["AMT"])
    nmbat2 <- rep("",times=nsim)
    nmbat2 <- nm.prep(nsim,ncore.m3,nid,EST.file,3,SIM.name.out,limdata,limobs,trunc.blq,ctlm3,nmbat2,vec["AMT"])
    setwd(EST.dir)
    cd.EST.dir <- paste("cd",EST.dir,sep=" ")

    for(i in 1:ncore.m1) {
      nmbat.split1 <- nmbat1[1:(nsim/ncore.m1)+(nsim/ncore.m1)*(i-1)]
  	  nmbat.split1 <- c(paste("call ",wfn.dir,sep=""),"E:",cd.EST.dir,nmbat.split1)
  	  bat.name1 <- paste(EST.file,"_M1",i,".bat",sep="")
  	  cmd <- paste(SIM.name.out,"_M1",i,".bat",sep="")
  	  writeLines(nmbat.split1,bat.name1)
  	  system(cmd, input=nmbat1, invisible=F, show.output.on.console=F, wait=F)
    }

    for (i in 1:ncore.m3) {
      nmbat.split2 <- nmbat2[1:(nsim/ncore.m3)+(nsim/ncore.m3)*(i-1)]
  	  nmbat.split2 <- c(paste("call ",wfn.dir,sep=""),"E:",cd.EST.dir,nmbat.split2)
  	  bat.name2 <- paste(EST.file,"_M3",i,".bat",sep="")
  	  cmd <- paste(SIM.name.out,"_M3",i,".bat",sep="")
      writeLines(nmbat.split2,bat.name2)
  	  if(i!=ncore.m3) {
        system(cmd, input=nmbat2, invisible=F, show.output.on.console=F, wait=F)
      }else{
        Sys.sleep(60)
        system(cmd, input=nmbat2, invisible=F, show.output.on.console=F, wait=F)
        wait.file <- paste0(SIM.name.out,"3_model",nsim)
        start.time <- Sys.time()
        while(!file.exists(paste0(EST.dir,"/",wait.file,".nm7/",tolower(wait.file),".fit"))) {
          Sys.sleep(60)
          print(Sys.time() - start.time)
        }
  	  }
    }
    setwd(master.dir)

  # Create r.script that will bring up all those lovely numbers you want
    input.vector1 <- c(blq,ruv.blq,vec["amt"],master.dir,SIM.name.out,ncore.m1,ncore.m3,vec["nsim"],
      vec["CL.POP"],vec["V2.POP"],vec["V3.POP"],vec["Q.POP"],vec["KA.POP"],vec["F1.POP"],
      vec["CL.BSV"],vec["V2.BSV"],vec["V3.BSV"],vec["Q.BSV"],vec["KA.BSV"],vec["F1.BSV"],
      vec["CL.BOV"],vec["V2.BOV"],vec["V3.BOV"],vec["Q.BOV"],ruv.prop,ruv.add,trunc.blq,
  		"functions_NCAvNLME_2016.r", "NONMEM_ref_M1.ctl", "NONMEM_ref_M3.ctl",
      "c(1,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,1)",
      nid,vec["LIMITLO"],vec["LIMITHI"],"c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96)",
      "seq(from=0, to=3, by=0.05)", "seq(from=3.25, to=6, by=0.25)",
      "seq(from=6.5, to=12, by=0.5)", "seq(from=13, to=clastlim, by=1)")
    phrase.vector1 <- c("SUB01", "SUB02", "SUB03", "SUB04", "SUB05", "SUB07", "SUB08",
      "SUB09", "SUB10", "SUB11", "SUB12", "SUB13", "SUB14", "SUB15", "SUB16", "SUB17",
      "SUB18", "SUB19", "SUB20", "SUB21", "SUB22", "SUB23", "SUB24", "SUB25", "SUB26",
      "SUB27", "SUB28", "SUB29", "SUB30", "SUB31", "SUB32", "SUB33", "SUB34", "SUB35",
      "SUB36", "SUB37", "SUB38", "SUB39", "SUB40")
    complete.r <- ctl.rep(phrase.vector1,input.vector1,rtemplate)
    writeLines(complete.r,paste(SIM.file,".R",sep=""))
    print(paste(SIM.name.out,"complete"))
  }, vec = runvec, time = timevec, cor = corvec, rtemplate = rscript)

### SECOND HALF ----------------------------------------------------------------

  ddply(rundf[1,], .(RUN, SCEN), function(df, vec) {
  #Set working directory
    SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
    SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
    SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
    EST.dir <- paste(SIM.dir,"ctl",sep="/")
    EST.file <- paste(EST.dir,SIM.name.out,sep="/")
    FIT.dir <- paste(SIM.dir,"fit",sep="/")

  # Define study design
    nid <- vec["NID"]
    nsim <- vec["NSIM"]
    nsub <- nid*nsim
    nobs <- length(time)
    sstime <- ifelse(df$SS.TYPE == 1,
      c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96),
      c(0,0.25,0.5,1,2,4,8,16,36,96))

  # Define random unexplained variability (SIGMA) values
    ruv.prop <- df$RUV.PROP
    ruv.blq <- df$RUV.BLQ
    blq <- df$BLQ
    if(df$RUV.TYPE == 1) {
      ruv.add <- (ruv.blq - ruv.prop) * blq
      trunc.blq <- blq
    }
    if(df$RUV.TYPE == 2) {
      ruv.add <- (0.2 - ruv.prop) * blq
      trunc.blq <- ruv.add/(ruv.blq - ruv.prop)
    }
    ruv.add.nm <- ifelse(blq == 0 || ruv.add == 0,
      paste(ruv.add, "FIX"),
      ruv.add
    )

  # Load data files for processing
    #simdata <- read.csv(paste(SIM.file,"_RAW.csv", sep="")) #Only do this if you need it for troubleshooting, RAW.csv can be 500MB+

    ipredresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_IPREDresult.csv",sep=""))
    ipredfresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_IPREDFresult.csv",sep=""))
    ipredcresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_IPREDCresult.csv",sep=""))

    ncaresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_NCAresult.csv",sep=""))
    ncafresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_NCAFresult.csv",sep=""))
    ncacresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_NCACresult.csv",sep=""))

    trunc.file <- paste(SIM.name.out,"TRUNCATED.csv",sep="_")
    limdata <- read.csv(paste(SIM.dir,trunc.file,sep="/"))
    limobs <- length(sstime)
    per.bloq <- percent.blq(limdata$DV,limdata$TIME,trunc.blq)

  # Setup mbt .bat file (is run manually)
    setwd(EST.dir)
    cd.EST.dir <- paste("cd",EST.dir,sep=" ")
    mbtcall <- c(paste("call ",wfn.dir,sep=""),"E:",cd.EST.dir,"call nmmbt")
    mbtbat <- "nmmbtrun.bat"
    mbtbat.dir <- paste(EST.dir,mbtbat,sep="/")
    writeLines(mbtcall,mbtbat.dir)
    system(mbtbat,invisible=FALSE)
    setwd(master.dir)

   # Process fit files into results table (see functions utility)
    m1.nlme.fit.out <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M1",nsim,nid,limobs)
    m1.nsim <- m1.nlme.fit.out[1]
    m1.fitfail <- m1.nlme.fit.out[-1]
    m1sim.in <- read.csv(paste(SIM.file,"M1_NMTHETAS.csv", sep="_"))
    m1sim.time <- nlme.simtime(m1sim.in)
    m1sim.out <- nlme.sim(m1sim.in,m1sim.time,nid,m1.nsim)
    data.process(m1sim.out,m1sim.time,nid,m1.nsim,SIM.file,trunc.blq,mode=3)
    m1result <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M1result.csv",sep=""))
    m1fresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M1Fresult.csv",sep=""))
    m1cresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M1Cresult.csv",sep=""))
    m1phfresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M1PHFresult.csv",sep=""))

    m3.nlme.fit.out <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M3",nsim,nid,limobs)
    m3.nsim <- m3.nlme.fit.out[1]
    m3.fitfail <- m3.nlme.fit.out[-1]
    m3sim.in <- read.csv(paste(SIM.file,"M3_NMTHETAS.csv", sep="_"))
    m3sim.time <- nlme.simtime(m3sim.in)
    m3sim.out <- nlme.sim(m3sim.in,m3sim.time,nid,m3.nsim)
    data.process(m3sim.out,m3sim.time,nid,m3.nsim,SIM.file,trunc.blq,mode=4)
    m3result <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3result.csv",sep=""))
    m3fresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3Fresult.csv",sep=""))
    m3cresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3Cresult.csv",sep=""))
    m3phfresult <- read.csv(paste(SIM.dir,"/",SIM.name.out,"_M3PHFresult.csv",sep=""))

  ### Determine for bioequvalence
  # FREL
    ipredftable <- glohipercent.func(ipredfresult,"Simulation", "FREL")
    ncaftable <- glohipercent.func(ncafresult,"Non-Compartmental Analysis", "FREL")
    m1ftable <- glohipercent.func(m1fresult,"Non-Linear Mixed Effects M1 Frel", "FREL")
    m1phfresult.in <- m1phfresult
    colnames(m1phfresult.in) <- c("SIM_ID", "FREL_GMEAN", "FREL_GLO95", "FREL_GHI95")
    m1phftable <- glohipercent.func(m1phfresult.in,"Non-Linear Mixed Effects M1 Post-Hoc", "FREL")
    m3ftable <- glohipercent.func(m3fresult,"Non-Linear Mixed Effects M3 Frel", "FREL")
    m3phfresult.in <- m3phfresult
    colnames(m3phfresult.in) <- c("SIM_ID", "FREL_GMEAN", "FREL_GLO95", "FREL_GHI95")
    m3phftable <- glohipercent.func(m3phfresult.in,"Non-Linear Mixed Effects M3 Post-Hoc", "FREL")
    allftable <- rbind(ipredftable,ncaftable,m1ftable,m1phftable,m3ftable,m3phftable)

  # CMAX
    ipredctable <- glohipercent.func(ipredcresult,"Simulation", "CMAX")
    ncactable <- glohipercent.func(ncacresult,"Non-Compartmental Analysis", "CMAX")
    m1ctable <- glohipercent.func(m1cresult,"Non-Linear Mixed Effects M1", "CMAX")
    m3ctable <- glohipercent.func(m3cresult,"Non-Linear Mixed Effects M3", "CMAX")
    allctable <- rbind(ipredctable,ncactable,m1ctable,m3ctable)

  ### Confidence Intervals of Means and Confidence Intervals (ConfidenceInterval-ception)
  # FREL
    ipredfci <- confint.func(ipredfresult,"IPRED", "FREL")
    ncafci   <- confint.func(ncafresult,"NCA", "FREL")
    m1fci  <- confint.func(m1fresult,"M1 F1", "FREL")
    m1phfci <- confint.func(m1phfresult,"M1 PH", "FREL")
    m3fci  <- confint.func(m3fresult,"M3 F1", "FREL")
    m3phfci <- confint.func(m3phfresult,"M3 PH", "FREL")
  # CMAX
    ipredcci <- confint.func(ipredcresult,"IPRED", "CMAX")
    ncacci   <- confint.func(ncacresult,"NCA", "CMAX")
    m1cci  <- confint.func(m1cresult,"M1", "CMAX")
    m3cci  <- confint.func(m3cresult,"M3", "CMAX")
  # And combine
    confint.table <- as.data.frame(rbind(ipredfci,ncafci,m1fci,m1phfci,m3fci,m3phfci,ipredcci,ncacci,m1cci,m3cci))
    confint.table <- orderBy(~VARIABLE+ANALYSIS ,confint.table)
    write.csv(confint.table, file=paste(SIM.file,"CONF_INTtable.csv",sep="_"),row.names=FALSE)
    confint.table2 <- confint.table
    confint.table2$VARIABLE <- c(rep((1:3),each=6),rep((4:6),each=4))
    confint.table2$ORDER <- confint.table$CIHI95-confint.table$CILO95

  # Find percentage of
    mbt <- read.table(file=paste(EST.dir,"nmmbt.nm7.txt",sep="/"),header=TRUE)
    mbt <- orderBy(~Run,mbt)
    mbt$Run <- rep(1:2,each=nsim)
    m1mbt <- subset(mbt[c(1,5,6)],Run==1) #m1
    m1mbt$Min <- gsub("SUCCESSFUL",1,m1mbt$Min)
    m1mbt$Min <- gsub("TERMINATED",0,m1mbt$Min)
    m1mbt$Cov <- gsub("NONE",0,m1mbt$Cov)
    m1mbt$Cov <- gsub("OK",1,m1mbt$Cov)
    m1mbt$Cov <- gsub("ABORTED",0,m1mbt$Cov)
    m1mbt$Cov <- gsub("UNOBTAINABLE",0,m1mbt$Cov)
    m3mbt <- subset(mbt[c(1,5,6)],Run==2) #m3
    m3mbt$Min <- gsub("SUCCESSFUL",1,m3mbt$Min)
    m3mbt$Min <- gsub("TERMINATED",0,m3mbt$Min)
    m3mbt$Cov <- gsub("NONE",0,m3mbt$Cov)
    m3mbt$Cov <- gsub("OK",1,m3mbt$Cov)
    m3mbt$Cov <- gsub("ABORTED",0,m3mbt$Cov)
    m3mbt$Cov <- gsub("UNOBTAINABLE",0,m3mbt$Cov)
    m1min <- mean(as.numeric(m1mbt$Min))*100
    m1cov <- mean(as.numeric(m1mbt$Cov))*100
    m3min <- mean(as.numeric(m3mbt$Min))*100
    m3cov <- mean(as.numeric(m3mbt$Cov))*100

  # Write to tables
    ipredbioq <- bioq.func(ipredftable, limitlo, limithi, "IPRED")
    ncabioq <- bioq.func(ncaftable, limitlo, limithi, "NCA")
    m1f1bioq <- bioq.func(m1ftable, limitlo, limithi, "M1F1")
    m1phbioq <- bioq.func(m1phftable, limitlo, limithi, "M1PH")
    m3f1bioq <- bioq.func(m3ftable, limitlo, limithi, "M3F1")
    m3phbioq <- bioq.func(m3phftable, limitlo, limithi, "M3PH")
    bioqtable <- rbind(ipredbioq, ncabioq, m1f1bioq, m1phbioq, m3f1bioq, m3phbioq)
    write.csv(bioqtable, file = paste(SIM.file, "BIOQtable.csv", sep = "_"), row.names = F)

    ipredcrat <- crat.func(ipredctable, limitlo, limithi, "IPRED")
    ncacrat <- crat.func(ncactable, limitlo, limithi, "NCA")
    m1crat <- crat.func(m1ctable, limitlo, limithi, "M1")
    m3crat <- crat.func(m3ctable, limitlo, limithi, "M3")
    crattable <- rbind(ipredcrat, ncacrat, m1crat, m3crat)
    write.csv(crattable, file=paste(SIM.file, "CRATtable.csv", sep = "_"), row.names = F)

    ncatbioq <- errortype.process(ncabioq, ipredbioq, 0, nsim, "NCA")
    m1f1tbioq <- errortype.process(m1f1bioq, ipredbioq, m1.fitfail, nsim, "M1F1")
    m3f1tbioq <- errortype.process(m3f1bioq, ipredbioq, m3.fitfail, nsim, "M3F1")
    m1phtbioq <- errortype.process(m1phbioq, ipredbioq, m1.fitfail, nsim, "M1PH")
    m3phtbioq <- errortype.process(m3phbioq, ipredbioq, m3.fitfail, nsim, "M3PH")
    write.csv(ncatbioq, file = paste(SIM.file, "NCA_BIOQtable.csv", sep = "_"), row.names = F)
    write.csv(m1f1tbioq, file = paste(SIM.file, "M1F1_BIOQtable.csv", sep = "_"), row.names = F)
    write.csv(m3f1tbioq, file = paste(SIM.file, "M3F1_BIOQtable.csv", sep = "_"), row.names = F)
    write.csv(m1phtbioq, file = paste(SIM.file, "M1PH_BIOQtable.csv", sep = "_"), row.names = F)
    write.csv(m3phtbioq, file = paste(SIM.file, "M3PH_BIOQtable.csv", sep = "_"), row.names = F)

    ncatcrat <- errortype.process(ncacrat, ipredcrat, 0, nsim, "NCA")
    m1tcrat <- errortype.process(m1crat, ipredcrat, m1.fitfail, nsim, "M1")
    m3tcrat <- errortype.process(m3crat, ipredcrat, m3.fitfail, nsim, "M3")
    write.csv(ncatcrat, file=paste(SIM.file, "NCA_CRATtable.csv" ,sep="_"), row.names = F)
    write.csv(m1tcrat, file=paste(SIM.file, "M1_CRATtable.csv" ,sep="_"), row.names = F)
    write.csv(m3tcrat, file=paste(SIM.file, "M3_CRATtable.csv" ,sep="_"), row.names = F)

    finaltable <- data.frame(
      mean(ipredbioq$p), mean(ncabioq$p), mean(m1f1bioq$p), mean(m1phbioq$p), mean(m3f1bioq$p),
      mean(m3phbioq$p), mean(ipredcrat$p), mean(ncacrat$p), mean(m1crat$p), mean(m3crat$p),
      mean(as.numeric(ncatbioq$pT1)-1), mean(as.numeric(ncatbioq$pT2)-1),
      mean(as.numeric(m1f1tbioq$pT1)-1), mean(as.numeric(m1f1tbioq$pT2)-1),
      mean(as.numeric(m3f1tbioq$pT1)-1), mean(as.numeric(m3f1tbioq$pT2)-1),
      mean(as.numeric(m1phtbioq$pT1)-1), mean(as.numeric(m1phtbioq$pT2)-1),
      mean(as.numeric(m3phtbioq$pT1)-1), mean(as.numeric(m3phtbioq$pT2)-1),
      mean(as.numeric(ncatcrat$pT1)-1), mean(as.numeric(ncatcrat$pT2)-1),
      mean(as.numeric(m1tcrat$pT1)-1), mean(as.numeric(m1tcrat$pT2)-1),
      mean(as.numeric(m3tcrat$pT1)-1), mean(as.numeric(m3tcrat$pT2)-1),
      per.bloq, trunc.blq, ruv.prop, ruv.add,
      m1min, m3min, m1cov, m3cov, m1.nsim, m3.nsim)
    colnames(finaltable) <- c(
      "IPRED_PBIOQ", "NCA_PBIOQ", "M1F1_PBIOQ", "M1PH_PBIOQ", "M3F1_PBIOQ",
      "M3PH_PBIOQ", "IPRED_PCRAT", "NCA_PCRAT", "M1_PCRAT", "M3_PCRAT",
      "NCA_FT1", "NCA_FT2", "M1F1_FT1", "M1F1_FT2",
      "M3F1_FT1", "M3F1_FT2", "M1PH_FT1", "M1PH_FT2",
      "M3PH_FT1", "M3PH_FT2", "NCA_CT1", "NCA_CT2",
      "M1_CT1", "M1_CT2", "M3_CT1", "M3_CT2",
      "%BLOQ", "LLOQ", "RUVprop", "RUVadd",
      "%M1msuc", "%M3msuc", "%M1csuc", "%M3csuc", "m1nsim", "m3nsim")
    rownames(finaltable) <- SIM.name.out
    write.csv(finaltable, file=paste(SIM.file,"FINALtable.csv",sep="_"),row.names=FALSE)

  }, vec = runvec)
