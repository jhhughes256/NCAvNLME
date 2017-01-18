### ------------------------------ NLMEvNCA ------------------------------- ###
#          (NLME vs. NCA Bioequivalence Analysis Stimulation Tool)            #

# As used in:
#  Hughes JH, Upton RU, Foster DJ (2016)
#  Comparison of Non-Compartmental and Mixed Effect Modelling Methods for
#  Establishing Bioequivalence for the Case of Two Compartment Kinetics and
#  Censored Concentrations

# The code below enables the user to compare the ability of NLME program NONMEM
#   and the automated NCA methods used by WinNonlin to determine the
#   bioequivalence of a drug. The tool is made up of four files which are
#   provided in the supplementary material of the paper. These files are:
#     - main_script.r
#     - function_utility.r
#     - NONMEM_ref_M1.ctl
#     - NONMEM_ref_m3.ctl
#
# DEPENDENCIES:
#   - NONMEM installation (original tool was built with Version 7.2)
#   - Wings for NONMEM (original tool was built with Version 720)
#
# To use the code, simply alter the data.frames and vectors below as
#   specified by the comments.
#
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

### INPUTS FOR ALTERATION -----------------------------------------------------

# Source functions file and NONMEM .ctl reference txt
#   Enter your desired directory as the master.dir
#   The tool will create further directories within this directory for each
#   scenario to be tested. It will also save the final table containing
#   collated results into master directory (master.dir).
  master.dir <- "Please enter directory"
  setwd(master.dir)

#   If you have changed the names of the tool files alter them here
  source("functions_utility.r")
  ctlm1 <- readLines("NONMEM_ref_M1.ctl")  #No BSV or BOV on Q & V3
  ctlm3 <- readLines("NONMEM_ref_M3.ctl")

# This refers to the Wings for NONMEM batch file, located in your local/remote
#   NONMEM installation. The directory below is an example of what the path may
#   look like.
  wfn.dir <- "c:/nm73ifort/wfn7/bin/wfn.bat"

# Load libraries
  library(doBy) # used for ordering of data.frames
  library(reshape2)
  library(stringr)
  library(plyr) # used for repeating functions
  library(MASS) # used for covariance matrix
  library(MBESS) # used for covariance matrix

# Specify run characteristics
  #set up values that are variable between runs
  rundf <- data.frame(
    RUN = rep(1:7, each = 15), #rep(number.of.runs, each=number.of.scenarios)
    SCEN = rep(1:15, times = 7), #opposite of above
    RUV.TYPE = rep(c(rep(1, 3),rep(2, 2)),21), #changing calculation of RUV
    RUV.BLQ = rep(c(0.2, 0.15, 0.1, 0.1, 0.5), 21),
    # changing Frel of generic and Frel BSV
    F1.POP = rep(c(rep(1.0, 5), rep(0.9, 5), rep(1.11, 5)), 7),
    F1.BSV = rep(c(rep(0.1225, 5), rep(0.0484, 5), rep(0.0529, 5)), 7),
    BLQ = c(rep(0.01, 15), rep(0.1, 15),rep(0.01, 75)),  #base LLOQ
    RUV.PROP = c(rep(0.05, 30), rep(0.09, 15), rep(0.05, 60)),#proportional RUV
    # Run 4 - reduced sampling schedule
    SS.TYPE = c(rep(1, 45), rep(2, 15), rep(1, 45)),
    # Run 5 - 20% lower generic KA
    KA.TYPE = c(rep(1, 60), rep(2, 15), rep(1, 30)),
    # Run 6 & 7 - BOV testing
    BOV.TYPE = c(rep(1, 75), rep(2, 15), rep(3, 15)))

#  Set up values that are constant between runs
  runvec <- c(
    NID = 24,  #number of patients per study
    NSIM = 20,  #number of simulations
    LIMITLO = 0.8,  #lower limit for bioequivalence
    LIMITHI = 1.25,  #upper limit for bioequivlence
    NCORE.M1 = 5,  #number of simultaneous instances of NONMEM running M1
    NCORE.M3 = 20,  #number of simultaneous instances of NONMEM running M3
    AMT = 125,  #dose of drug
    # Population parameters (THETA)
    CL.POP = 20,  #clearance
    V2.POP = 100,  #central volume
    Q.POP = 35,  #central/peripheral distribution
    V3.POP = 400,  #peripheral volume
    KA.POP = 0.5,  #(innovator) absorption rate constant
    #   Only used if KA.TYPE is set to 2 - gives generic different absorption
    KA.GEN = 0.3,  #generic absorption rate constant
    # Between subject variability parameters (ETA)
    CL.BSV = 0.045,  #clearance
    V2.BSV = 0.045,  #central volume
    V3.BSV = 0.045,   #peripheral volume
    Q.BSV = 0.045,  #central/peripheral distribution
    KA.BSV = 0.01,  #absorption rate
    # Between occasion variability parameters (ETA)
    CL.BOV = 0.045,  #clearance
    V2.BOV = 0.045,  #central volume
    V3.BOV = 0.045,  #peripheral volume
    Q.BOV = 0.045,  #central/peripheral distribution
    #   Only used if BOV.TYPE is set to 2 or 3 - effectively removes BOV
    ALT.BOV = 0.0001)

  # Set up time vector for simulation
  timevec <- c(
    seq(0, 3, 0.05),
    seq(3.25, 6, 0.25),
    seq(6.5, 12, 0.5),
    seq(13, 96, 1))
  TIME <- timevec  #TIME required in .GlobalEnv for simulate.2comp.abs()

  # Set up time vector for truncated samples
  #   These are the times of the samples that will be analysed by NCA & NLME
  #   T1 is used if SS.TYPE == 1, T2 is used if SS.TYPE == 2
  sstimelist <- list(
    T1 = c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96),
    T2 = c(0,0.25,0.5,1,2,4,8,16,36,96))

  # Set up correlation vector for simulation
  corvec <- c(
    1, 0.3, 0.3, 0.3,
    0.3, 1, 0.3, 0.3,
    0.3, 0.3, 1, 0.3,
    0.3, 0.3, 0.3 ,1)

### END OF INPUTS FOR ALTERATION ----------------------------------------------

### FIRST HALF ----------------------------------------------------------------
  # First loop includes
  # 1. Simulation of data
  # 2. Initial analysis of simulated data
  # 3. Non-Compartmental Analysis
  # 4. Initial analysis of NCA data
  # 5. Non-Linear Mixed Effects

  ddply(rundf, .(RUN, SCEN), function(df, vec, time, cor, limtime) {
### Object names within ddply function
#     rundf -> df     (index using dollar sign df$index)
#     runvec -> vec   (index using square brackets ["index"])
#     timevec -> time
#     corvec -> cor
#     sstimelist -> timelist (index using dollar sign)

### 1. Simulation of data
  # Set the seed (very important!)
  #   Ensures that the same set of patients are used between scenarios
  #   Simply change the seed to observe a different set of patients
  #   For any given scenario:
  #     Total unique patients = number of sims * number of subjects per study
    set.seed(1234)

  # Set working directory and file names
    SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
    SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
    SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
    EST.dir <- paste(SIM.dir,"ctl",sep="/")
    EST.file <- paste(EST.dir,SIM.name.out,sep="/")
    FIT.dir <- paste(SIM.dir,"fit",sep="/")

  # Define study design
  #   This takes many of the values specified above
    nid <- vec["NID"]
    nsim <- vec["NSIM"]
    nsub <- nid*nsim  #total number of subjects
    nobs <- length(time)  #number of observations per simulation
  #   SS.TYPE determines which sampling time is used as stated above
    if (df$SS.TYPE == 1) {
      sstime <- limtime$T1
    } else {
      sstime <- limtime$T2
    }

  # Define random unexplained variability (SIGMA) values
  #   RUV.TYPE determines which method of determining additive RUV is used
    ruv.prop <- df$RUV.PROP
    ruv.blq <- df$RUV.BLQ
    blq <- df$BLQ
  #   additive RUV and CV at the LOQ change between Scenarios
  #   lloq stays constant
    if (df$RUV.TYPE == 1) {  #Scenarios 1-9
      ruv.add <- (ruv.blq - ruv.prop) * blq
      trunc.blq <- blq
    }
  #   lloq and CV at the LOQ change between Scenarios
  #   additive RUV stays constant
    if (df$RUV.TYPE == 2) {  #Scenarios 10-15
      ruv.add <- (0.2 - ruv.prop) * blq
      trunc.blq <- ruv.add/(ruv.blq - ruv.prop)
    }
  #   object to be placed into nonmem control stream template
    ruv.add.nm <- ifelse(blq == 0 || ruv.add == 0,
      paste(ruv.add, "FIX"),
      ruv.add
    )

  # Set cores
    ncore.m1 <- vec["NCORE.M1"]
    ncore.m3 <- vec["NCORE.M3"]

  # Correlation between parameters
    cormat <- matrix(cor, nrow = 4, ncol = 4)
    std.bsv <- c(
      sqrt(vec["CL.BSV"]),sqrt(vec["V2.BSV"]),
      sqrt(vec["Q.BSV"]),sqrt(vec["V3.BSV"]))
    omega <- cor2cov(cormat, std.bsv)

  # Account for random BSV
  #   Produce samples from multi-variate normal distribution as specified by
  #   omega block created above.
    etamat <- mvrnorm(n = nsub, mu = c(0, 0, 0, 0), omega)
    ETA1 <- rep(etamat[, 1], each = 2)  # CLbsv
    ETA2 <- rep(etamat[, 2], each = 2)  # V2bsv
    ETA3 <- rep(etamat[, 3], each = 2)  # Qbsv
    ETA4 <- rep(etamat[, 4], each = 2)  # V3bsv
  #   Produce samples from normal distribution as specifed in the function
  #   If statements for runs with differing BSV & BOV.TYPE
    if (vec["KA.BSV"] == 0) {
      ETA5 <- rep(0, times = nsub)  # KA (no bsv)
    } else {
      ETA5 <- rnorm(nsub, mean = 0, sd = sqrt(vec["KA.BSV"]))  # KA (bsv)
    }
    if (df$F1.BSV == 0) {
      ETA6 <- rep(0, times = nsub)  # F1 (no bsv)
    } else {
      ETA6 <- rnorm(nsub, mean = 0, sd = sqrt(df$F1.BSV))  # F1 (bsv)
    }
    if (df$BOV.TYPE != 2) {
      ETA7 <- rnorm(nsub * 2, mean = 0, sd = sqrt(vec["CL.BOV"]))  # CLbov
    } else {
      ETA7 <- rnorm(nsub * 2, mean = 0, sd = sqrt(vec["ALT.BOV"]))  # CLbov
    }
    if (df$BOV.TYPE == 1) {
      ETA8 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["V2.BOV"]))  # V2bov
      ETA9 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["Q.BOV"]))  # Qbov
      ETA10 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["V3.BOV"]))  # V3bov
    } else {
      ETA8 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["ALT.BOV"]))  # V2bov
      ETA9 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["ALT.BOV"]))  # Qbov
      ETA10 <- rnorm(nsub*2,mean=0,sd=sqrt(vec["ALT.BOV"]))  # V3bov
    }
    if (vec["KA.BSV"]==0){
      ETA11 <- rep(0, times=nsub)  # KA (no bsv)
    } else {
      ETA11 <- rnorm(nsub,mean=0,sd=sqrt(vec["KA.BSV"]))  # KA (bsv)
    }

  # Create/clear directories before run
    dir.setup(SIM.dir)
    dir.setup(EST.dir)
    dir.setup(FIT.dir)
    nm.clear(EST.dir, nsim*2)

  # Determine KA individual values -> dependent on KA.TYPE
    if (df$KA.TYPE == 1) {
      ka.val <- rep(vec["KA.POP"] * exp(ETA5), each = 2)
    } else {
      ka.val <- as.vector(rbind(vec["KA.POP"] * exp(ETA5), vec["KA.GEN"] * exp(ETA11)))
    }

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
  #   Order of operations specified below
    simdata <- orderBy(~UID + FORM + TIME,  #SECOND: order by these columns
      merge(simdataIPRED, simdataPRED, all = T)  #FIRST: merge dataframes
    )[  #THIRD: reorder columns
      c("UID", "ID", "SIM", "TIME", "AMT", "FORM",
      "CL", "V2", "Q", "V3", "KA", "F1", "PRED", "IPRED")
    ]

  # Add residual error
    CP <- simdata$IPRED
    tnobs <- length(CP)
    EPS1 <- rnorm(tnobs, mean = 0, sd = ruv.prop)
    EPS2 <- rnorm(tnobs, mean = 0, sd = ruv.add)
    simdata$DV <- CP * (1 + EPS1) + EPS2
  # Remove residual error from TIME == 0
    simdata$DV <- ifelse(simdata$TIME == 0, 0, simdata$DV)

  # Save simulation file
    write.csv(simdata, file = paste0(SIM.file, "_RAW.csv"), row.names = F)

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
  # Replaces phrases in the .ctl reference file to values stated at the
  #   beginning of the script
    ctlm1 <- ctl.prep(ctlm1, df, vec, omega, ruv.prop, ruv.add.nm, trunc.blq)
    ctlm3 <- ctl.prep(ctlm3, df, vec, omega, ruv.prop, ruv.add.nm, trunc.blq)

  # Create .bat files for NONMEM shell and run in NONMEM (number of .bats
  #   determined by ncore value)
    nmbat1 <- rep("", times =nsim)
    nmbat1 <- nm.prep(
      nsim, ncore.m1, nid, EST.file, 1, SIM.name.out,
      limdata, limobs, trunc.blq, ctlm1, nmbat1, vec["AMT"])
    nmbat2 <- rep("", times = nsim)
    nmbat2 <- nm.prep(
      nsim, ncore.m3, nid, EST.file, 3, SIM.name.out,
      limdata, limobs, trunc.blq, ctlm3, nmbat2, vec["AMT"])
    setwd(EST.dir)
    cd.EST.dir <- paste("cd", EST.dir)

  # Run batch files in multiple instances of NONMEM
    for(i in 1:ncore.m1) {
      nmbat.split1 <- nmbat1[1:(nsim/ncore.m1)+(nsim/ncore.m1)*(i-1)]
  	  nmbat.split1 <- c(paste0("call ",wfn.dir), "E:", cd.EST.dir, nmbat.split1)
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
  	  if (i!=ncore.m3) {
        system(cmd, input=nmbat2, invisible=F, show.output.on.console=F, wait=F)
    # If final batch file wait to prevent other NONMEM instances finishing first
      } else {
        Sys.sleep(300)
        system(cmd, input=nmbat2, invisible=F, show.output.on.console=F, wait=F)
        wait.file <- paste0(SIM.name.out,"3_model",nsim)
        start.time <- Sys.time()
    # Wait for final batch file to produce .fit file
        while(!file.exists(paste0(
          EST.dir,"/",wait.file,".nm7/",tolower(wait.file),".smr"))) {
            Sys.sleep(60)
            print(Sys.time() - start.time)
        }
  	  }
    }
    setwd(master.dir)
    print(paste(SIM.name.out,"complete"))
  }, vec = runvec, time = timevec, cor = corvec, limtime = sstimelist)

### SECOND HALF ----------------------------------------------------------------

  bioqtable <- ddply(rundf, .(RUN, SCEN), function(df, vec, time, limtime) {
  # Set working directory
    SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
    SIM.dir <- paste(master.dir,SIM.name.out,sep="/")k

    SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
    EST.dir <- paste(SIM.dir,"ctl",sep="/")
    EST.file <- paste(EST.dir,SIM.name.out,sep="/")
    FIT.dir <- paste(SIM.dir,"fit",sep="/")

  # Define study design
    nid <- vec["NID"]
    nsim <- vec["NSIM"]
    nsub <- nid*nsim
    nobs <- length(time)
    if (df$SS.TYPE == 1) {
      sstime <- limtime$T1
    } else {
      sstime <- limtime$T2
    }

  # Define random unexplained variability (SIGMA) values
    ruv.prop <- df$RUV.PROP
    ruv.blq <- df$RUV.BLQ
    blq <- df$BLQ
    if (df$RUV.TYPE == 1) {
      ruv.add <- (ruv.blq - ruv.prop) * blq
      trunc.blq <- blq
    }
    if (df$RUV.TYPE == 2) {
      ruv.add <- (0.2 - ruv.prop) * blq
      trunc.blq <- ruv.add/(ruv.blq - ruv.prop)
    }
    ruv.add.nm <- ifelse(blq == 0 || ruv.add == 0,
      paste(ruv.add, "FIX"),
      ruv.add)

  # Load data files for processing
  #   Only load RAW if you need it for troubleshooting, RAW.csv can be 500MB+
    #simdata <- read.csv(paste(SIM.file,"_RAW.csv", sep=""))

    ipredresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_IPREDresult.csv"))
    ipredfresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_IPREDFresult.csv"))
    ipredcresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_IPREDCresult.csv"))

    ncaresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_NCAresult.csv"))
    ncafresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_NCAFresult.csv"))
    ncacresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_NCACresult.csv"))

    trunc.file <- paste(SIM.name.out,"TRUNCATED.csv",sep="_")
    limdata <- read.csv(paste(SIM.dir,trunc.file,sep="/"))
    limobs <- length(sstime)  #number of observations for truncated sampling times
    per.bloq <- percent.blq(limdata$DV,limdata$TIME,trunc.blq)  #percent BLOQ

  # Setup mbt .bat file
    setwd(EST.dir)
    cd.EST.dir <- paste("cd", EST.dir)
    mbtcall <- c(paste("call", wfn.dir), "E:", cd.EST.dir, "call nmmbt")
    mbtbat <- "nmmbtrun.bat"
    mbtbat.dir <- paste(EST.dir, mbtbat, sep="/")
    writeLines(mbtcall, mbtbat.dir)  #save .bat file
    system(mbtbat, show.output.on.console=F)  #run .bat file
    setwd(master.dir)

  # Process fit files into results table
    m1nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M1",nsim,nid,limobs)
    m1.nsim <- m1nlme.fitout[1]  #number of completed runs
    m1.nsub <- m1.nsim*nid
    m1.fitfail <- m1nlme.fitout[-1]  #list of failed runs (not terminated)
    m1sim.in <- read.csv(paste(SIM.file,"M1_NMTHETAS.csv", sep="_"))
    m1sim.time <- nlme.simtime(m1sim.in)
    #Simulation using estimated parameters from models
    m1sim.out <- nlme.sim(m1sim.in,m1sim.time,nid,m1.nsim)
    data.process(m1sim.out,m1sim.time,nid,m1.nsim,SIM.file,trunc.blq,mode=3)
    #Load processed simulation output
    m1result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M1result.csv"))
    m1fresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M1Fresult.csv"))
    m1cresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M1Cresult.csv"))
    m1phfresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M1PHFresult.csv"))

    m3nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M3",nsim,nid,limobs)
    m3.nsim <- m3nlme.fitout[1]  #number of completed runs
    m3.nsub <- m3.nsim*nid
    m3.fitfail <- m3nlme.fitout[-1]  #list of failed runs (not terminated)
    m3sim.in <- read.csv(paste(SIM.file, "M3_NMTHETAS.csv", sep="_"))
    m3sim.time <- nlme.simtime(m3sim.in)
    #Simulation using estimated parameters from models
    m3sim.out <- nlme.sim(m3sim.in,m3sim.time,nid,m3.nsim)
    data.process(m3sim.out,m3sim.time,nid,m3.nsim,SIM.file,trunc.blq,mode=4)
    #Load processed simulation output
    m3result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M3result.csv"))
    m3fresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M3Fresult.csv"))
    m3cresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M3Cresult.csv"))
    m3phfresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M3PHFresult.csv"))

  # Find percentage of successful runs
  #   Done by loading up nmmbt run made previously in the script.
  #   Scrape output to determine which models minimised successfully
    mbt <- read.table(file=paste(EST.dir,"nmmbt.nm7.txt",sep="/"),header=TRUE)
    mbt <- orderBy(~Run,mbt)
    mbt$Method <- rep(1:2,each=nsim)
    mbt$ModelNum <- rep(order(as.character(1:nsim)), times = 2)
    mbt$Success <- gsub("SUCCESSFUL",1,mbt$Min)
    mbt$Success <- gsub("TERMINATED",0,mbt$Success)
    mbt <- orderBy(~Method+ModelNum,mbt)
    #list of terminated runs
    m1.term <- which(mbt$Success == 0 & mbt$Method == 1)
    m3.term <- which(mbt$Success == 0 & mbt$Method == 2) - nsim

    #determine number of successfully minimised runs + covariance steps
    m1mbt <- mbt[mbt$Method == 1, c(1,5,6)] #m1
    m1mbt$Min <- gsub("SUCCESSFUL",1,m1mbt$Min)
    m1mbt$Min <- gsub("TERMINATED",0,m1mbt$Min)
    m1mbt$Cov <- gsub("NONE",0,m1mbt$Cov)
    m1mbt$Cov <- gsub("OK",1,m1mbt$Cov)
    m1mbt$Cov <- gsub("ABORTED",0,m1mbt$Cov)
    m1mbt$Cov <- gsub("UNOBTAINABLE",0,m1mbt$Cov)
    m3mbt <- mbt[mbt$Method == 2, c(1,5,6)] #m3
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

  # Determine for bioequvalence
  #   Create a data.frame ready for use with ddply
    aovprep <- data.frame(matrix(NA, nrow=nsub*4+m1.nsub*4+m3.nsub*4, ncol=8))
    colnames(aovprep) <- c("METH","ID","SIM","FORM","AUC","CMAX","IDf","FORMf")
    aovprep$METH <- as.factor(c(
      rep("IPRED", nsub*2),  #IPRED
      rep("NCA", nsub*2),  #NCA
      rep("M1PH", m1.nsub*2),  #M1PH
      rep("M1F1", m1.nsub*2),  #M1F1
      rep("M3PH", m3.nsub*2),  #M3PH
      rep("M3F1", m3.nsub*2)))  #M3F1
    aovprep$ID <- c(
      limdata[limdata$TIME == 0, 2],  #IPRED
      rep(ncaresult$STUD_ID, each = 2),  #NCA
      rep(m1sim.in$STUD_ID,2),  #M1PH & M1F1
      rep(m3sim.in$STUD_ID,2))  #M3PH & M3F1
    aovprep$SIM <- c(
      limdata[limdata$TIME == 0, 3],  #IPRED
      rep(ncaresult$SIM_ID, each = 2),  #NCA
      if (m1.fitfail == 0) {  #M1PH & M1F1
        rep(m1sim.in$SIM_ID,2)
      } else {
        rep((1:nsim)[-m1.fitfail], each = nid*2, times = 2)
      },
      if (m3.fitfail == 0) { #M3PH & M1F3
        rep(m3sim.in$SIM_ID,2)
      } else {
        rep((1:nsim)[-m3.fitfail], each = nid*2, times = 2)
      })
    aovprep$FORM <- c(
      limdata[limdata$TIME == 0, 6],  #IPRED
      rep(1:2, times = nsub),  #NCA
      rep(1:2, times = m1.nsub*2), #M1PH & M1F1
      rep(1:2, times = m3.nsub*2)) #M3PH &M3F1
    aovprep$AUC <- c(
      limdata[limdata$TIME == 0, 12],  #IPRED
      as.vector(rbind(ncaresult$INN_AUC, ncaresult$GEN_AUC)),  #NCA
      as.vector(rbind(m1result$INN_AUC, m1result$GEN_AUC)),  #M1PH
      m1sim.in$F1,  #M1F1
      as.vector(rbind(m3result$INN_AUC, m3result$GEN_AUC)),  #M3PH
      m3sim.in$F1)  #M3F1
    aovprep$CMAX <- c(
      as.vector(rbind(ipredresult$INN_CMAX,ipredresult$GEN_CMAX)),  #IPRED
      as.vector(rbind(ncaresult$INN_CMAX, ncaresult$GEN_CMAX)),  #NCA
      rep(as.vector(rbind(m1result$INN_CMAX, m1result$GEN_CMAX)),2),  #M1
      rep(as.vector(rbind(m3result$INN_CMAX, m3result$GEN_CMAX)),2))  #M3
    aovprep$IDf <- as.factor(aovprep$ID)
    aovprep$FORMf <- as.factor(aovprep$FORM)

  # Use ANOVA to determine bioequivalent studies
    faov <- ddply(aovprep, .(METH, SIM), function(df) runaov2(df, USE = "AUC"))
    caov <- ddply(aovprep, .(METH, SIM), function(df) runaov2(df, USE = "CMAX"))

  # Calculate proportion of bioequivalent studies
    fbioq <- ddply(faov, .(METH), function(df) mean(df$BE))
    cbioq <- ddply(caov, .(METH), function(df) mean(df$BE))

  # Remove simulations with termination status not successful
    faov.termstat <- faov[!(
      faov$METH == "M1F1" & faov$SIM %in% m1.term |
      faov$METH == "M1PH" & faov$SIM %in% m1.term |
      faov$METH == "M3F1" & faov$SIM %in% m3.term |
      faov$METH == "M3PH" & faov$SIM %in% m3.term), ]
    caov.termstat <- caov[!(
      caov$METH == "M1F1" & caov$SIM %in% m1.term |
      caov$METH == "M3F1" & caov$SIM %in% m3.term |
      caov$METH == "M1PH" | caov$METH == "M3PH"), ]

  # Separate methods from true bioequivalent studies (IPRED)
    meth.faov <- faov.termstat[faov.termstat$METH != "IPRED", ]
    meth.faov$IPRED.BE <- c(
      rep(
        faov.termstat$BE[faov.termstat$METH == "IPRED" &
        !faov.termstat$SIM %in% m1.term &
        !faov.termstat$SIM %in% m1.fitfail], 2),
      rep(
        faov.termstat$BE[faov.termstat$METH == "IPRED" &
        !faov.termstat$SIM %in% m3.term &
        !faov.termstat$SIM %in% m3.fitfail], 2),
      faov.termstat$BE[faov.termstat$METH == "IPRED"])

    meth.caov <- caov.termstat[caov.termstat$METH != "IPRED", ]
    meth.caov$IPRED.BE <- c(
      caov.termstat$BE[caov.termstat$METH == "IPRED" &
        !caov.termstat$SIM %in% m1.term &
        !caov.termstat$SIM %in% m1.fitfail],
      caov.termstat$BE[caov.termstat$METH == "IPRED" &
        !caov.termstat$SIM %in% m3.term &
        !caov.termstat$SIM %in% m3.fitfail],
      caov.termstat$BE[caov.termstat$METH == "IPRED"])

  # Determine percentage of bioequivalent studies
    fbioq <- ddply(faov.termstat, .(METH), function(df) mean(df$BE))
    cbioq <- ddply(caov.termstat, .(METH), function(df) mean(df$BE))
    ipred.fbioq <- ddply(meth.faov, .(METH), function(df) mean(df$IPRED.BE))
    ipred.cbioq <- ddply(meth.caov, .(METH), function(df) mean(df$IPRED.BE))

  # Determine which studies have made a Type I or Type II error
    ferror <- ddply(meth.faov, .(METH),
      function(df) errortype.func2(df$IPRED.BE, df$BE))
    cerror <- ddply(meth.caov, .(METH),
      function(df) errortype.func2(df$IPRED.BE, df$BE))

  # Determine percentage of studies with Type I or Type II error
    ferror.t1 <- ddply(ferror, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT1))))
    ferror.t2 <- ddply(ferror, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT2))))
    cerror.t1 <- ddply(cerror, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT1))))
    cerror.t2 <- ddply(cerror, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT2))))

  # Print to console that processing is complete
    print(paste(SIM.name.out,"processed"))

    if (m1.fitfail != 0) {
      m1.term <- unique(c(m1.fitfail, m1.term))
    }
    if (m3.fitfail != 0) {
      m3.term <- unique(c(m3.fitfail, m3.term))
    }

    m1.ifterm <- length(m3.term) >= nsim
    m3.ifterm <- length(m3.term) >= nsim

  # Collate data table for output
    aovbioqtable <- data.frame(
      IPRED.BE = fbioq$V1[fbioq$METH == "IPRED"]*100,
      NCA.BE = fbioq$V1[fbioq$METH == "NCA"]*100,
      M1F1.BE = ifelse(!m1.ifterm,
        fbioq$V1[fbioq$METH == "M1F1"]*100,0),
      M1PH.BE = ifelse(!m1.ifterm,
        fbioq$V1[fbioq$METH == "M1PH"]*100,0),
      M3F1.BE = ifelse(!m3.ifterm,
        fbioq$V1[fbioq$METH == "M3F1"]*100,0),
      M3PH.BE = ifelse(!m3.ifterm,
        fbioq$V1[fbioq$METH == "M3PH"]*100,0),
      IPRED.CM = cbioq$V1[cbioq$METH == "IPRED"]*100,
      NCA.CM = cbioq$V1[cbioq$METH == "NCA"]*100,
      M1.CM = ifelse(!m1.ifterm,
        cbioq$V1[cbioq$METH == "M1F1"]*100,0),
      M3.CM = ifelse(!m3.ifterm,
        cbioq$V1[cbioq$METH == "M3F1"]*100,0),
      NCAF.T1 = ferror.t1$V1[ferror.t1$METH == "NCA"]*100,
      M1F1F.T1 = ifelse(!m1.ifterm,
        ferror.t1$V1[ferror.t1$METH == "M1F1"]*100,0),
      M1PHF.T1 = ifelse(!m1.ifterm,
        ferror.t1$V1[ferror.t1$METH == "M1PH"]*100,0),
      M3F1F.T1 = ifelse(!m3.ifterm,
        ferror.t1$V1[ferror.t1$METH == "M3F1"]*100,0),
      M3PHF.T1 = ifelse(!m3.ifterm,
        ferror.t1$V1[ferror.t1$METH == "M3PH"]*100,0),
      NCAF.T2 = ferror.t2$V1[ferror.t2$METH == "NCA"]*100,
      M1F1F.T2 = ifelse(!m1.ifterm,
        ferror.t2$V1[ferror.t2$METH == "M1F1"]*100,0),
      M1PHF.T2 = ifelse(!m1.ifterm,
        ferror.t2$V1[ferror.t2$METH == "M1PH"]*100,0),
      M3F1F.T2 = ifelse(!m3.ifterm,
        ferror.t2$V1[ferror.t2$METH == "M3F1"]*100,0),
      M3PHF.T2 = ifelse(!m3.ifterm,
        ferror.t2$V1[ferror.t2$METH == "M3PH"]*100,0),
      NCAC.T1 = cerror.t1$V1[cerror.t1$METH == "NCA"]*100,
      M1C.T1 = ifelse(!m1.ifterm,
        cerror.t1$V1[cerror.t1$METH == "M1F1"]*100,0),
      M3C.T1 = ifelse(!m3.ifterm,
        cerror.t1$V1[cerror.t1$METH == "M3F1"]*100,0),
      NCAC.T2 = cerror.t2$V1[cerror.t2$METH == "NCA"]*100,
      M1C.T2 = ifelse(!m1.ifterm,
        cerror.t2$V1[cerror.t2$METH == "M1F1"]*100,0),
      M3C.T2 = ifelse(!m3.ifterm,
        cerror.t2$V1[cerror.t2$METH == "M3F1"]*100,0),
      PERBLOQ = per.bloq,
      TRUNCBLQ = trunc.blq,
      RUVPROP = ruv.prop,
      RUVADD = ruv.add,
      NSIM = nsim,
      M1MIN = m1min,
      M3MIN = m3min,
      M1COV = m1cov,
      M3COV = m3cov,
      M1NSIM = m1.nsim,
      M1SUCC = nsim - length(m1.term),
      M1IPREDBE = ifelse(!m1.ifterm,
        ipred.fbioq$V1[ipred.fbioq$METH == "M1F1"]*100,0),
      M1IPREDCM = ifelse(!m1.ifterm,
        ipred.cbioq$V1[ipred.cbioq$METH == "M1F1"]*100,0),
      M3NSIM = m3.nsim,
      M3SUCC = nsim - length(m3.term),
      M3IPREDBE = ifelse(!m3.ifterm,
        ipred.fbioq$V1[ipred.fbioq$METH == "M3F1"]*100,0),
      M3IPREDCM = ifelse(!m3.ifterm,
        ipred.cbioq$V1[ipred.cbioq$METH == "M3F1"]*100,0))
    aovbioqtable
  }, vec = runvec, time = timevec, limtime = sstimelist)
  write.csv(bioqtable,
    file = paste(master.dir, "collated_bioq_table.csv", sep = "/"),
    row.names = F)

  error.df <- bioqtable[c(1:2,13:28)]
  error.df.l <- melt(error.df, c("RUN", "SCEN"))
  error.df.w <- dcast(data.frame(
    RUN = error.df.l$RUN,
    SCEN = error.df.l$SCEN,
    colsplit(error.df.l$variable, "\\.", c("Method", "Error")),
    value = error.df.l$value), RUN+SCEN+Method ~ Error)

  ferror.df <- error.df.w[str_detect(error.df.w$Method, "F"), ]
  ferror.df$NSIM <- as.vector(rbind(
    bioqtable$M1SUCC, bioqtable$M1SUCC,
    bioqtable$M3SUCC, bioqtable$M3SUCC,
    bioqtable$NSIM
  ))
  f.ipredbe <- as.vector(rbind(
    bioqtable$M1IPREDBE, bioqtable$M1IPREDBE,
    bioqtable$M3IPREDBE, bioqtable$M3IPREDBE,
    bioqtable$IPRED.BE
  ))
  ferror.df$NBIOQ <- f.ipredbe/100*c(
    rep(bioqtable$M1SUCC, each = 2), rep(bioqtable$M3SUCC, each = 2), bioqtable$NSIM
  )

  cerror.df <- error.df.w[!str_detect(error.df.w$Method, "F"), ]
  cerror.df$NSIM <- as.vector(rbind(
    bioqtable$M1SUCC,
    bioqtable$M3SUCC,
    bioqtable$NSIM
  ))
  c.ipredbe <- as.vector(rbind(
    bioqtable$M1IPREDCM,
    bioqtable$M3IPREDCM,
    bioqtable$IPRED.CM
  ))
  cerror.df$NBIOQ <- c.ipredbe/100*c(
    bioqtable$M1SUCC, bioqtable$M3SUCC, bioqtable$NSIM
  )

  ferror.final.l <- melt(
    ddply(ferror.df, .(RUN, SCEN, Method), function (x) error.matrix.fun(x)),
    c("RUN", "SCEN", "Method"))
  ferror.final.l$variable <- paste0(ferror.final.l$Method, ferror.final.l$variable)
  ferror.final.l$Method <- NULL
  ferror.factors <- unique(ferror.final.l$variable)
  ferror.final.l$variable <- factor(ferror.final.l$variable, levels = ferror.factors)
  ferror.final <- dcast(ferror.final.l, RUN+SCEN ~ variable)

  cerror.final.l <- melt(
    ddply(cerror.df, .(RUN, SCEN, Method), function (x) error.matrix.fun(x)),
    c("RUN", "SCEN", "Method"))
  cerror.final.l$variable <- paste0(cerror.final.l$Method, cerror.final.l$variable)
  cerror.final.l$Method <- NULL
  cerror.factors <- unique(cerror.final.l$variable)
  cerror.final.l$variable <- factor(cerror.final.l$variable, levels = cerror.factors)
  cerror.final <- dcast(cerror.final.l, RUN+SCEN ~ variable)

  file.name <- "_collated_results.csv"
  write.csv(ferror.final, paste0("Frel", file.name))
  write.csv(cerror.final, paste0("Cmax", file.name))
