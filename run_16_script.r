### ---------------------------- NLMEvNCA BEAST ---------------------------- ###
#           (NLME vs. NCA Bioequivalence Analysis Stimulation Tool)            #

# As used in:
#  Hughes JH, Upton RU, Foster DJ (2016)
#  Comparison of Non-Compartmental and Mixed Effect Modelling Methods for
#  Establishing Bioequivalence for the Case of Two Compartment Kinetics and
#  Censored Concentrations


# The code below enables the user to compare the ability of NLME program NONMEM and
#   the automated NCA methods used by WinNonlin to determine the bioequivalence of a
#   drug. To run the code you require the five files provided in the supplementary
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Source functions file and NONMEM .ctl reference txt
  master.dir <- "E:/hscpw-df1/Data1/Jim Hughes/2016/01_MDV" ### Directory containing source files
  setwd(master.dir)
  source("functions_NCAvNLME_2016.r")
  ctlm1 <- readLines("NONMEM_ref_M1.ctl")  #No BSV or BOV on Q & V3
  ctlm3 <- readLines("NONMEM_ref_M3.ctl")
  rscript <- readLines("Rscript.r")
  wfn.dir <- "c:/nm72/wfn7/bin/wfn.bat"

# Load libraries
  library(doBy)
  library(plyr)
  library(MASS)
  library(MBESS)
  library(stringr)
  library(reshape2)
  set.seed(1234)

# Specify run characteristics
  #set up values that are variable between runs
  run.order <- c(rep(1, times = 9), rep(2, times = 6))
  f1.order <- c(1, 0.9, 1.11)
  bsv.order <- c(0.1225, 0.0484, 0.0529)
  rundf <- data.frame(
    RUN = unlist(llply(0:2, function(i) run.order + i*2)), #rep(number.of.runs, each=number.of.scenarios)
    SCEN = rep(c(1:9,1:6), times = 3), #opposite of above
    RUV.TYPE = rep(c(rep(1, 9),rep(2, 6)), 3),
    RUV.BLQ = rep(c(rep(c(0.2, 0.15, 0.1), 3),rep(c(0.1, 0.5), 3)), 3),
    F1.POP = rep(c(rep(f1.order, each = 3), rep(f1.order, each = 2)), 3), #changing Frel of generic
    F1.BSV = rep(c(rep(bsv.order, each = 3), rep(bsv.order, each = 2)), 3), #changing BSV on frel
    BLQ = c(rep(0.01, 15), rep(0.1, 15), rep(0.01, 15)),  #Run 2 - raised LLOQ
    RUV.PROP = c(rep(0.05, 30), rep(0.09, 15)),  #Run 3 - increased proportional RUV
    SS.TYPE = rep(1, 45),  #Run 4 - reduced sampling schedule
    KA.TYPE = rep(1, 45),  #Run 5 - 20% lower generic KA
    BOV.TYPE = rep(1, 45))  #Run 6 & 7 - BOV testing

  #set up values that are constant between runs
  runvec <- c(
    NID = 24,
    NSIM = 500,
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
    Q.BOV = 0.045,
    ALT.BOV = 0.0001)

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

### SECOND HALF ----------------------------------------------------------------

  bioqtable <- ddply(rundf, .(RUN, SCEN), function(df, vec, time) {
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
    if (df$SS.TYPE == 1) {
      sstime <- c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96)
    } else {
      sstime <- c(0,0.25,0.5,1,2,4,8,16,36,96)
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
    #simdata <- read.csv(paste(SIM.file,"_RAW.csv", sep="")) #Only do this if you need it for troubleshooting, RAW.csv can be 500MB+
    ipredresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_IPREDresult.csv"))
    ncaresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_NCAresult.csv"))

    trunc.file <- paste(SIM.name.out,"TRUNCATED.csv",sep="_")
    limdata <- read.csv(paste(SIM.dir,trunc.file,sep="/"))
    limobs <- length(sstime)
    per.bloq <- percent.blq(limdata$DV,limdata$TIME,trunc.blq)

   # Process fit files into results table (see functions utility)
    m1nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M1",nsim,nid,limobs)
    m1.nsim <- m1nlme.fitout[1]
    m1.nsub <- m1.nsim*nid
    m1.fitfail <- m1nlme.fitout[-1]
    m1sim.in <- read.csv(paste(SIM.file,"M1_NMTHETAS.csv", sep="_"))
    m1result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M1result.csv"))

    m3nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M3",nsim,nid,limobs)
    m3.nsim <- m3nlme.fitout[1]
    m3.nsub <- m3.nsim*nid
    m3.fitfail <- m3nlme.fitout[-1]
    m3sim.in <- read.csv(paste(SIM.file,"M3_NMTHETAS.csv", sep="_"))
    m3result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M3result.csv"))

  # Find percentage of successful runs
    mbt <- read.table(file=paste(EST.dir,"nmmbt.nm7.txt",sep="/"),header=TRUE)
    mbt <- orderBy(~Run,mbt)
    mbt$Method <- rep(1:2,each=nsim)
    mbt$ModelNum <- rep(order(as.character(1:nsim)), times = 2)
    mbt$Success <- gsub("SUCCESSFUL",1,mbt$Min)
    mbt$Success <- gsub("TERMINATED",0,mbt$Min)
    mbt <- orderBy(~Method+ModelNum,mbt)
    m1.term <- which(mbt$Success == 0 & mbt$Method == 1)
    m3.term <- which(mbt$Success == 0 & mbt$Method == 2) - nsim

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

  ### Determine for bioequvalence
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

    faov <- ddply(aovprep, .(METH, SIM), function(df) runaov2(df, USE = "AUC"))
    caov <- ddply(aovprep, .(METH, SIM), function(df) runaov2(df, USE = "CMAX"))
    meth.faov <- faov[faov$METH != "IPRED", ]
    meth.faov$IPRED.BE <- c(
      rep(faov$BE[faov$METH == "IPRED" & !faov$SIM %in% m1.fitfail], 2),
      rep(faov$BE[faov$METH == "IPRED" & !faov$SIM %in% m3.fitfail], 2),
      faov$BE[faov$METH == "IPRED"])
    meth.caov <- caov[caov$METH != "IPRED", ]
    meth.caov$IPRED.BE <- c(
      rep(caov$BE[caov$METH == "IPRED" & !caov$SIM %in% m1.fitfail], 2),
      rep(caov$BE[caov$METH == "IPRED" & !caov$SIM %in% m3.fitfail], 2),
      caov$BE[caov$METH == "IPRED"])

    fbioq <- ddply(faov, .(METH), function(df) mean(df$BE))
    cbioq <- ddply(caov, .(METH), function(df) mean(df$BE))
    ipred.fbioq <- ddply(meth.faov, .(METH), function(df) mean(df$IPRED.BE))
    ipred.cbioq <- ddply(meth.caov, .(METH), function(df) mean(df$IPRED.BE))

    ferror <- ddply(meth.faov, .(METH), function(df)
      errortype.func2(df$IPRED.BE, df$BE))

    ferror.t1 <- ddply(ferror, .(METH), function(df)
      mean(as.numeric(as.vector(df$pT1))))
    ferror.t2 <- ddply(ferror, .(METH), function(df)
      mean(as.numeric(as.vector(df$pT2))))

    cerror <- ddply(meth.caov, .(METH), function(df)
      errortype.func2(df$IPRED.BE, df$BE))

    cerror.t1 <- ddply(cerror, .(METH), function(df)
      mean(as.numeric(as.vector(df$pT1))))
    cerror.t2 <- ddply(cerror, .(METH), function(df)
      mean(as.numeric(as.vector(df$pT2))))

    faov.termstat <- faov[!(
      faov$METH == "M1F1" & faov$SIM %in% m1.term |
      faov$METH == "M1PH" & faov$SIM %in% m1.term |
      faov$METH == "M3F1" & faov$SIM %in% m3.term |
      faov$METH == "M3PH" & faov$SIM %in% m3.term), ]
    caov.termstat <- caov[!(
      caov$METH == "M1F1" & caov$SIM %in% m1.term |
      caov$METH == "M3F1" & caov$SIM %in% m3.term |
      caov$METH == "M1PH" | caov$METH == "M3PH"), ]

    meth.faov.termstat <- faov.termstat[faov.termstat$METH != "IPRED", ]
    meth.faov.termstat$IPRED.BE <- c(
      rep(
        faov.termstat$BE[faov.termstat$METH == "IPRED" &
        !faov.termstat$SIM %in% m1.term &
        !faov.termstat$SIM %in% m1.fitfail], 2),
      rep(
        faov.termstat$BE[faov.termstat$METH == "IPRED" &
        !faov.termstat$SIM %in% m3.term &
        !faov.termstat$SIM %in% m3.fitfail], 2),
      faov.termstat$BE[faov.termstat$METH == "IPRED"])

    meth.caov.termstat <- caov.termstat[caov.termstat$METH != "IPRED", ]
    meth.caov.termstat$IPRED.BE <- c(
      caov.termstat$BE[caov.termstat$METH == "IPRED" &
        !caov.termstat$SIM %in% m1.term &
        !caov.termstat$SIM %in% m1.fitfail],
      caov.termstat$BE[caov.termstat$METH == "IPRED" &
        !caov.termstat$SIM %in% m3.term &
        !caov.termstat$SIM %in% m3.fitfail],
      caov.termstat$BE[caov.termstat$METH == "IPRED"])

    fbioq.termstat <- ddply(faov.termstat, .(METH),
      function(df) mean(df$BE))
    cbioq.termstat <- ddply(caov.termstat, .(METH),
      function(df) mean(df$BE))
    ipred.fbioq.termstat <- ddply(meth.faov.termstat, .(METH),
      function(df) mean(df$IPRED.BE))
    ipred.cbioq.termstat <- ddply(meth.caov.termstat, .(METH),
      function(df) mean(df$IPRED.BE))

    ferror.termstat <- ddply(meth.faov.termstat, .(METH),
      function(df) errortype.func2(df$IPRED.BE, df$BE))

    ferror.termstat.t1 <- ddply(ferror.termstat, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT1))))
    ferror.termstat.t2 <- ddply(ferror.termstat, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT2))))

    cerror.termstat <- ddply(meth.caov.termstat, .(METH),
      function(df) errortype.func2(df$IPRED.BE, df$BE))

    cerror.termstat.t1 <- ddply(cerror.termstat, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT1))))
    cerror.termstat.t2 <- ddply(cerror.termstat, .(METH),
      function(df) mean(as.numeric(as.vector(df$pT2))))

    print(paste(SIM.name.out,"processed"))

    aovbioqtable <- data.frame(
      TERMSTAT = c("All","Only Success"),
      IPREDBE = c(fbioq$V1[1]*100,fbioq.termstat$V1[1]*100),
      NCABE = c(fbioq$V1[6]*100,fbioq.termstat$V1[6]*100),
      M1F1BE = c(fbioq$V1[2]*100,fbioq.termstat$V1[2]*100),
      M1PHBE = c(fbioq$V1[3]*100,fbioq.termstat$V1[3]*100),
      M3F1BE = c(fbioq$V1[4]*100,fbioq.termstat$V1[4]*100),
      M3PHBE = c(fbioq$V1[5]*100,fbioq.termstat$V1[5]*100),
      IPREDCM = c(cbioq$V1[1]*100,cbioq.termstat$V1[1]*100),
      NCACM = c(cbioq$V1[6]*100,cbioq.termstat$V1[4]*100),
      M1CM = c(cbioq$V1[2]*100,cbioq.termstat$V1[2]*100),
      M3CM = c(cbioq$V1[4]*100,cbioq.termstat$V1[3]*100),
      NCAFT1 = c(ferror.t1$V1[5]*100,ferror.termstat.t1$V1[5]*100),
      M1F1FT1 = c(ferror.t1$V1[1]*100,ferror.termstat.t1$V1[1]*100),
      M1PHFT1 = c(ferror.t1$V1[2]*100,ferror.termstat.t1$V1[2]*100),
      M3F1FT1 = c(ferror.t1$V1[3]*100,ferror.termstat.t1$V1[3]*100),
      M3PHFT1 = c(ferror.t1$V1[4]*100,ferror.termstat.t1$V1[4]*100),
      NCAFT2 = c(ferror.t2$V1[5]*100,ferror.termstat.t2$V1[5]*100),
      M1F1FT2 = c(ferror.t2$V1[1]*100,ferror.termstat.t2$V1[1]*100),
      M1PHFT2 = c(ferror.t2$V1[2]*100,ferror.termstat.t2$V1[2]*100),
      M3F1FT2 = c(ferror.t2$V1[3]*100,ferror.termstat.t2$V1[3]*100),
      M3PHFT2 = c(ferror.t2$V1[4]*100,ferror.termstat.t2$V1[4]*100),
      NCACT1 = c(cerror.t1$V1[5]*100,cerror.termstat.t1$V1[3]*100),
      M1CT1 = c(cerror.t1$V1[1]*100,cerror.termstat.t1$V1[1]*100),
      M3CT1 = c(cerror.t1$V1[3]*100,cerror.termstat.t1$V1[2]*100),
      NCACT2 = c(cerror.t2$V1[5]*100,cerror.termstat.t2$V1[3]*100),
      M1CT2 = c(cerror.t2$V1[1]*100,cerror.termstat.t2$V1[1]*100),
      M3CT2 = c(cerror.t2$V1[3]*100,cerror.termstat.t2$V1[2]*100),
      PERBLOQ = per.bloq,
      TRUNCBLQ = trunc.blq,
      RUVPROP = ruv.prop,
      RUVADD = ruv.add,
      M1MIN = m1min,
      M3MIN = m3min,
      M1COV = m1cov,
      M3COV = m3cov,
      M1NSIM = c(m1.nsim, m1.nsim - length(m1.term)),
      M1IPREDBE = c(ipred.fbioq$V1[1]*100,ipred.fbioq.termstat$V1[1]*100),
      M1IPREDCM = c(ipred.cbioq$V1[1]*100,ipred.cbioq.termstat$V1[1]*100),
      M3NSIM = c(m3.nsim, m3.nsim - length(m3.term)),
      M3IPREDBE = c(ipred.fbioq$V1[3]*100,ipred.fbioq.termstat$V1[3]*100),
      M3IPREDCM = c(ipred.cbioq$V1[3]*100,ipred.cbioq.termstat$V1[3]*100))
    aovbioqtable
  }, vec = runvec, time = timevec)
  write.csv(bioqtable,
    file = paste(master.dir, "collated_bioq_table.csv", sep = "/"),
    row.names = F)

    collate.SHK <- function(dir.name, work.dir) {
      nm.dir <- "nm7"
      shk.file.name <- gsub(nm.dir, "shk", dir.name)
      shk.file.path <- paste(work.dir, dir.name, shk.file.name, sep="/")
      #Scrape data from the *.shk file
      if (file.exists(shk.file.path) == T) {  #screen for missing file
        shk.data <- read.table(shk.file.path, skip = 1, header = T)  #read all the lines of the shk file
    	  shk.data <- subset(shk.data, SUBPOP == 1)  #Only take SUBPOP 1
        shrink.vals <- subset(shk.data, TYPE == 4)
      } else {
        shrink.vals <- rep(0, 10)
      }
      data.frame(
        Sim = if(!is.na(str_extract(dir.name, "l[0-9]{3}"))) {
          str_extract(str_extract(dir.name, "l[0-9]{3}"), "[0-9]{3}")
        } else if (!is.na(str_extract(dir.name, "l[0-9]{2}"))) {
          str_extract(str_extract(dir.name, "l[0-9]{2}"), "[0-9]{2}")
        } else if (!is.na(str_extract(dir.name, "l[0-9]"))) {
          str_extract(str_extract(dir.name, "l[0-9]"), "[0-9]")
        } else {
          NA
        },
        Method = str_extract(str_extract(dir.name, "[0-9]_m"), "[0-9]"),
        shrink.vals[-(1:2)])
    }

    nm.dir <- "nm7"
    search.term <- paste("*",nm.dir, sep="")
    shrink.data <- ddply(rundf, .(RUN, SCEN), function(df) {
    #Set working directory
      SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
      SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
      SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
      EST.dir <- paste(SIM.dir,"ctl",sep="/")
      EST.file <- paste(EST.dir,SIM.name.out,sep="/")
      FIT.dir <- paste(SIM.dir,"fit",sep="/")
      dir.names <- dir(path=EST.dir,pattern=glob2rx(search.term))
      shk.out <- mdply(dir.names, collate.SHK, work.dir = EST.dir)
      print(paste(SIM.name.out,"processed"))
      shk.out
    })
    shk.filename <- paste(master.dir, "collated_shrinkage_data.csv", sep = "/")
    write.csv(shrink.data,
      file = shk.filename,
      row.names = F)

    shrink.data <- arrange(read.csv(shk.filename), RUN, SCEN, Method, Sim)
    shrink.data$X1 <- NULL

    shrink.sum <- ddply(rundf, .(RUN, SCEN), function(df, vec, time, shk) {
    #Set working directory
      SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
      SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
      SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
      EST.dir <- paste(SIM.dir,"ctl",sep="/")
      EST.file <- paste(EST.dir,SIM.name.out,sep="/")
      FIT.dir <- paste(SIM.dir,"fit",sep="/")

      nid <- vec["NID"]
      nsim <- vec["NSIM"]
      nsub <- nid*nsim
      nobs <- length(time)
      if (df$SS.TYPE == 1) {
        sstime <- c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96)
      } else {
        sstime <- c(0,0.25,0.5,1,2,4,8,16,36,96)
      }

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

      sub.shk <- shk[shk$RUN == df$RUN & shk$SCEN == df$SCEN, ]
      l.shk <- dim(sub.shk)[2]
      names(sub.shk)[5:l.shk] <- c(
        "BSVCL", "BSVV2", "BSVKA", "BSVF1",
        "BOVCL1", "BOVCL2", "BOVV21", "BOVV22")

    #Load data for processing
      #simdata <- read.csv(paste(SIM.file,"_RAW.csv", sep="")) #Only do this if you need it for troubleshooting, RAW.csv can be 500MB+
      #ipredresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_IPREDresult.csv"))
      #ncaresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_NCAresult.csv"))

      trunc.file <- paste(SIM.name.out,"TRUNCATED.csv",sep="_")
      limdata <- read.csv(paste(SIM.dir,trunc.file,sep="/"))
      limobs <- length(sstime)
      per.bloq <- percent.blq(limdata$DV,limdata$TIME,trunc.blq)

     # Process fit files into results table (see functions utility)
      m1nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M1",nsim,nid,limobs)
      m1.nsim <- m1nlme.fitout[1]
      m1.nsub <- m1.nsim*nid
      m1.fitfail <- m1nlme.fitout[-1]
      #m1sim.in <- read.csv(paste(SIM.file,"M1_NMTHETAS.csv", sep="_"))
      #m1result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M1result.csv"))

      m3nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M3",nsim,nid,limobs)
      m3.nsim <- m3nlme.fitout[1]
      m3.nsub <- m3.nsim*nid
      m3.fitfail <- m3nlme.fitout[-1]
      #m3sim.in <- read.csv(paste(SIM.file,"M3_NMTHETAS.csv", sep="_"))
      #m3result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M3result.csv"))

    # Find percentage of successful runs
      mbt <- read.table(file=paste(EST.dir,"nmmbt.nm7.txt",sep="/"),header=TRUE)
      mbt <- orderBy(~Run,mbt)
      mbt$Method <- rep(1:2,each=nsim)
      mbt$ModelNum <- rep(order(as.character(1:nsim)), times = 2)
      mbt$Success <- gsub("SUCCESSFUL",1,mbt$Min)
      mbt$Success <- gsub("TERMINATED",0,mbt$Min)
      mbt <- orderBy(~Method+ModelNum,mbt)
      m1.term <- which(mbt$Success == 0 & mbt$Method == 1)
      m3.term <- which(mbt$Success == 0 & mbt$Method == 2) - nsim

      term.shk <- sub.shk[
        sub.shk$Method == 1 &
          !sub.shk$Sim %in% m1.fitfail &
          !sub.shk$Sim %in% m1.term |
        sub.shk$Method == 3 &
          !sub.shk$Sim %in% m3.fitfail &
          !sub.shk$Sim %in% m3.term, ]

      sub.shk$Term <- "All"
      term.shk$Term <- "Success"
      comb.shk <- rbind(sub.shk, term.shk)

      #insert alternate ending here
    }, vec = runvec, time = timevec, shk = shrink.data)
    shrink.sum <- arrange(shrink.sum, RUN, SCEN, Method, Sim)
    shk.filename <- paste(master.dir, "collatedterm_shrinkage_data.csv", sep = "/")
    write.csv(shrink.sum,
      file = shk.filename,
      row.names = F)
      #removing everything between

  collate.ERR <- function(dir.name, work.dir) {
    nm.dir <- "nm7"
    lst.file.name <- gsub(nm.dir, "lst", dir.name)
    lst.file.path <- paste(work.dir, dir.name, lst.file.name, sep = "/")

    #Scrape data from the *.lst file
    if (file.exists(lst.file.path) == T) {  #screen for missing file
      lst.lines <- readLines(lst.file.path)
      first.line <- grep("#TERM:", lst.lines)
      last.line <- grep(" NO. OF SIG. DIGITS", lst.lines)

      if (length(last.line) == 0 | length(first.line) == 0) {
        info <- c("Function has not detected", "a termination status in the",
        ".lst file. Please check", "models.lst file manually.", dir.name)
        browser("Unknown error, type info for more information, type Q to quit.")
      } else {
        error.lines <- lst.lines[first.line:last.line[1]]
      }
      err.code <- "None"
      desc.code <- "None"

      #Determine termination status and error status
      if (length(grep("TERMINATED", error.lines[2])) != 0) {
        term.code <- "Unsuccessful"
        if (length(grep("ROUNDING", error.lines)) != 0) {
          err.code <- "Rounding errors"
          if (length(grep("0PARAMETER ESTIMATE IS NEAR", lst.lines[last.line[1] + 1])) != 0) {
            desc.code <- "Parameter estimate near boundary"
          } else {
            desc.code <- "Rounding errors"
          }
        } else if (length(grep("ZERO GRADIENT", error.lines)) != 0) {
          err.code <- "Zero gradient"
          desc.code <- "Zero gradient"
        } else if (length(grep("OBJ. FUNC. IS INFINITE", error.lines)) != 0) {
          err.code <- "Obj close to infinity"
          if (length(grep("LAST ITERATION", error.lines)) != 0) {
            desc.code <- "Due to last iteration"
          } else if (length(grep("NEXT ITERATION", error.lines)) != 0) {
            desc.code <- "Due to next iteration"
          } else {
            desc.code <- error.lines[3]
          }
        } else if (length(grep("MAX. NO. OF FUNCTION", error.lines)) != 0) {
          err.code <- "Reached max evaluations"
          if (length(grep("0PARAMETER ESTIMATE IS NEAR", lst.lines[last.line[1] + 1])) != 0) {
            desc.code <- "Parameter estimate near boundary"
          }
        } else {
          err.code <- error.lines[3]
        }
      } else if (length(grep("SUCCESSFUL", error.lines[2])) != 0) {
        term.code <- "Successful"
        if (length(grep("HOWEVER, PROBLEMS", error.lines[3])) != 0) {
          desc.code <- "Problems occurred during minimisation"
        }
        if (length(grep("0PARAMETER ESTIMATE IS NEAR", lst.lines[last.line[1] + 1])) != 0) {
          desc.code <- "COV Parameter estimate near boundary"
        }
        if (length(grep("0R MATRIX", lst.lines)) != 0) {
          desc.code <- "R matrix algorithmically singular"
        }
        if (length(grep("0S MATRIX", lst.lines)) != 0) {
          desc.code <- "S matrix unobtainable"
        }
        if (length(grep("0PRED EXIT CODE = 1", lst.lines)) != 0) {
          desc.code <- "Problems with individual"
        }
      } else {
        term.code <- error.lines[2]
      }

      #Determine covariance step status
      cov.line <- grep("STANDARD ERROR OF ESTIMATE", lst.lines)
      if (length(cov.line) > 0 & term.code == "Successful") {
        cov.code <- "Passed"
        if(desc.code != "None") desc.code <- paste("COVPASS",desc.code)
      } else if (length(cov.line) == 0 | term.code == "Unsuccessful") {
        cov.code <- "Failed"
        if(desc.code != "None" & term.code == "Successful") desc.code <- paste("COVFAIL",desc.code)
      }
    }
    data.frame(
      Model = dir.name,
      TermStat = term.code,
      CovCode = cov.code,
      ErrCode = err.code,
      TermDesc = desc.code)
  }

  nm.dir <- "nm7"
  search.term <- paste("*",nm.dir, sep="")
  error.data <- ddply(rundf, .(RUN, SCEN), function(df) {
  #Set working directory
    SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
    SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
    SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
    EST.dir <- paste(SIM.dir,"ctl",sep="/")
    EST.file <- paste(EST.dir,SIM.name.out,sep="/")
    FIT.dir <- paste(SIM.dir,"fit",sep="/")
    dir.names <- dir(path=EST.dir,pattern=glob2rx(search.term))
    err.out <- mdply(dir.names, collate.ERR, work.dir = EST.dir)
    print(paste(SIM.name.out,"processed"))
    err.out
  })
  err.filename <- paste(master.dir, "collated_error_data.csv", sep = "/")
  write.csv(error.data,
    file = err.filename,
    row.names = F)

  shrink.data <- arrange(read.csv(shk.filename), RUN, SCEN, Method, Sim)
  shrink.data$X1 <- NULL

  shrink.sum <- ddply(rundf, .(RUN, SCEN), function(df, vec, time, shk) {
  #Set working directory
    SIM.name.out <- paste0("Run", df$RUN, "_Scen", df$SCEN)
    SIM.dir <- paste(master.dir,SIM.name.out,sep="/")
    SIM.file <- paste(SIM.dir,SIM.name.out,sep="/")
    EST.dir <- paste(SIM.dir,"ctl",sep="/")
    EST.file <- paste(EST.dir,SIM.name.out,sep="/")
    FIT.dir <- paste(SIM.dir,"fit",sep="/")

    nid <- vec["NID"]
    nsim <- vec["NSIM"]
    nsub <- nid*nsim
    nobs <- length(time)
    if (df$SS.TYPE == 1) {
      sstime <- c(0,0.25,0.5,1,2,4,6,8,12,16,24,36,48,72,96)
    } else {
      sstime <- c(0,0.25,0.5,1,2,4,8,16,36,96)
    }

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

    sub.shk <- shk[shk$RUN == df$RUN & shk$SCEN == df$SCEN, ]
    l.shk <- dim(sub.shk)[2]
    names(sub.shk)[5:l.shk] <- c(
      "BSVCL", "BSVV2", "BSVKA", "BSVF1",
      "BOVCL1", "BOVCL2", "BOVV21", "BOVV22")

  #Load data for processing
    #simdata <- read.csv(paste(SIM.file,"_RAW.csv", sep="")) #Only do this if you need it for troubleshooting, RAW.csv can be 500MB+
    #ipredresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_IPREDresult.csv"))
    #ncaresult <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_NCAresult.csv"))

    trunc.file <- paste(SIM.name.out,"TRUNCATED.csv",sep="_")
    limdata <- read.csv(paste(SIM.dir,trunc.file,sep="/"))
    limobs <- length(sstime)
    per.bloq <- percent.blq(limdata$DV,limdata$TIME,trunc.blq)

   # Process fit files into results table (see functions utility)
    m1nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M1",nsim,nid,limobs)
    m1.nsim <- m1nlme.fitout[1]
    m1.nsub <- m1.nsim*nid
    m1.fitfail <- m1nlme.fitout[-1]
    #m1sim.in <- read.csv(paste(SIM.file,"M1_NMTHETAS.csv", sep="_"))
    #m1result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M1result.csv"))

    m3nlme.fitout <- nlme.fit(SIM.name.out,FIT.dir,EST.dir,"M3",nsim,nid,limobs)
    m3.nsim <- m3nlme.fitout[1]
    m3.nsub <- m3.nsim*nid
    m3.fitfail <- m3nlme.fitout[-1]
    #m3sim.in <- read.csv(paste(SIM.file,"M3_NMTHETAS.csv", sep="_"))
    #m3result <- read.csv(paste0(SIM.dir,"/",SIM.name.out,"_M3result.csv"))

  # Find percentage of successful runs
    mbt <- read.table(file=paste(EST.dir,"nmmbt.nm7.txt",sep="/"),header=TRUE)
    mbt <- orderBy(~Run,mbt)
    mbt$Method <- rep(1:2,each=nsim)
    mbt$ModelNum <- rep(order(as.character(1:nsim)), times = 2)
    mbt$Success <- gsub("SUCCESSFUL",1,mbt$Min)
    mbt$Success <- gsub("TERMINATED",0,mbt$Min)
    mbt <- orderBy(~Method+ModelNum,mbt)
    m1.term <- which(mbt$Success == 0 & mbt$Method == 1)
    m3.term <- which(mbt$Success == 0 & mbt$Method == 2) - nsim

    term.shk <- sub.shk[
      sub.shk$Method == 1 &
        !sub.shk$Sim %in% m1.fitfail &
        !sub.shk$Sim %in% m1.term |
      sub.shk$Method == 3 &
        !sub.shk$Sim %in% m3.fitfail &
        !sub.shk$Sim %in% m3.term, ]

    sub.shk$Term <- "All"
    term.shk$Term <- "Success"
    comb.shk <- rbind(sub.shk, term.shk)

    #insert alternate ending here
  }, vec = runvec, time = timevec, shk = shrink.data)
  shrink.sum <- arrange(shrink.sum, RUN, SCEN, Method, Sim)
  shk.filename <- paste(master.dir, "collatedterm_shrinkage_data.csv", sep = "/")
  write.csv(shrink.sum,
    file = shk.filename,
    row.names = F)
    #removing everything between
# Alternate ending
#      ddply(comb.shk, .(Term, Method), function(x) {
#        y <- data.frame(
#          Stat = c("Min", "Q1", "Med", "Mean", "Q3", "Max"),
#          colwise(summary)(x[5:l.shk]))
#        z <- dcast(melt(y, "Stat"), variable~Stat)[c(1,5,6,4,3,7,2)]
#      })
#    }, vec = runvec, time = timevec, shk = shrink.data)
#    shrink.sum <- arrange(shrink.sum, RUN, SCEN, Method, variable)
#    shk.filename <- paste(master.dir, "summaryall_shrinkage_data.csv", sep = "/")
#    write.csv(shrink.sum,
#      file = shk.filename,
#      row.names = F)
#
#    shrink.sum <- arrange(shrink.sum, RUN, Method, variable)
#    shk.filename <- paste(master.dir, "summarysub_shrinkage_data.csv", sep = "/")
#    write.csv(shrink.sum[shrink.sum$Term == "Success", ],
#      file = shk.filename,
#      row.names = F)
