# Functions required by the code found in NCA_v_NLME.r
# Create/Clear directories before a run
  # Warnings are suppressed so warnings about file.paths already existing + permission for deletion of folders denied are not printed to console

  dir.setup <- function(dir) {
    suppressWarnings(dir.create(file.path(dir)))
    suppressWarnings(do.call(file.remove,list(list.files(dir,full.names=TRUE))))
  }

# Clear NONMEM folders
  nm.clear <- function(dir, nsim)
  {
    NM.dir.all <- list.files(dir,full.names=TRUE)
    for (i in 1:nsim)
    {
      NM.dir <- NM.dir.all[i]
      suppressWarnings(do.call(file.remove,list(list.files(NM.dir,full.names=TRUE))))
    }
  }
#--------------------------------------------------------------------------------------------
# SIMULATION FUNCTIONS
# Simulation function for 2 compartmental model with first order absorption kinetics
  simulate.2comp.abs <- function(ID,AMT,CL,Q,V2,V3,KA,F1) {
    k10 <- CL/V2
	  k12 <- Q/V2
	  k21 <- Q/V3
	  apb <- k10+k12+k21            # alpha + beta
	  amb <- k10*k21                # alpha * beta
    alpha <- ((apb)+sqrt((apb)^2-4*amb))/2
	  beta <- ((apb)-sqrt((apb)^2-4*amb))/2
    A <- KA*(k21-alpha)/(V2*(KA-alpha)*(beta-alpha))
	  B <- KA*(k21-beta)/(V2*(KA-beta)*(alpha-beta))
	  Cplasma <- AMT*F1*(A*exp(-alpha*TIME)+B*exp(-beta*TIME)-(A+B)*exp(-KA*TIME))
	  result <- data.frame("TIME"=TIME, "CP"=Cplasma)
	  result
  }
#--------------------------------------------------------------------------------------------
# PROCESSING FUNCTIONS (for IPRED, NCA, and for fit files from NLME)
data.process <- function(simdata,sstime,nid,nsim,SIM.file,blq,mode)  { # 1 <- IPRED   2 <- NCA    3 <- M1    4 <- M3
  #Setup both mode-dependent and independent objects
  nobs <- length(sstime)
  nsub <- nsim*nid
  if(mode==1) {
    file.tag <- "_IPRED"
	  col.names <- c("UNQ_ID","STUD_ID","SIM_ID","CL","Q","V2","V3","KA","INN_CMAX","GEN_CMAX","INN_TMAX","GEN_TMAX","FREL","CRATIO")
	  rcol <- length(col.names)
    df.in <- simdata
  }
  if(mode==2) {
    file.tag <- "_NCA"
	  col.names <- c("UNQ_ID","STUD_ID","SIM_ID","CL","Q","V2","V3","KA","INN_CMAX","GEN_CMAX","INN_TMAX","GEN_TMAX","INN_AUC","GEN_AUC","FREL","CRATIO")
	  rcol <- length(col.names)
    totalobs <- length(unique(simdata$TIME))
	  ssrows  <- match(sstime,simdata$TIME)
    limtimes <- as.vector(outer(ssrows, 0:(nsub * 2 - 1) * length(TIME), "+"))
    df.in <- simdata[limtimes,]
    rownames(df.in) <- NULL
    write.csv(df.in, file=paste(SIM.file,"TRUNCATED.csv", sep="_"), row.names=F)
  }
  if(mode>=3) {
	  if(mode==3) {file.tag <- "_M1"}
	  if(mode==4) {file.tag <- "_M3"}
	  col.names <- c("UNQ_ID","STUD_ID","SIM_ID","CL1","CL2","Q1","Q2","V2_1","V2_2","V3_1","V3_2","KA","INN_CMAX","GEN_CMAX","INN_TMAX","GEN_TMAX","INN_AUC","GEN_AUC","FREL","PHFREL","CRATIO")
	  rcol <- length(col.names)
    df.in <- simdata
  }

  #Setting up data.frames as empty as opposed to rbinding to save computation
  result <- data.frame(matrix(NA, nrow = nsub, ncol = rcol))
  colnames(result) <- col.names
  fresult <- data.frame(matrix(NA, nrow = nsim, ncol = 4))
  colnames(fresult) <- c("SIM_ID","FREL_GMEAN","FREL_GLO95","FREL_GHI95")
  cresult <- data.frame(matrix(NA, nrow = nsim, ncol = 4))
  colnames(cresult) <- c("SIM_ID","CMAX_GMEAN","CMAX_GLO95","CMAX_GHI95")
  if(mode>=3) {
    phfresult <- data.frame(matrix(NA, nrow = nsim, ncol = 4))
    colnames(phfresult) <- c("SIM_ID","PHFREL_GMEAN","PHFREL_GLO95","PHFREL_GHI95")
  }

  #Begin loop
  for(i in 1:nsim) {
  # Subset the data
    subdf   <- subset(df.in, SIM==i)
	  formdf1 <- subset(subdf, subdf$FORM==1)
    formdf2 <- subset(subdf, subdf$FORM==2)

  # Reshape table -> column names are time, dv."1:maxid"
  # Create data frame only containing dependent variable
    if(mode==2) { #if NCA
      time1 <- c(formdf1$TIME)
      id1 <- c(formdf1$ID)
      dv1 <- c(formdf1$DV)
      inndf <- data.frame(id1, time1, dv1)
      inndf <- reshape(inndf,idvar="time1",timevar="id1",direction="wide")

      time2 <- c(formdf2$TIME)
      id2 <- c(formdf2$ID)
      dv2 <- c(formdf2$DV)
      gendf <- data.frame(id2, time2, dv2)
      gendf <- reshape(gendf,idvar="time2",timevar="id2",direction="wide")

	# Censor data using blq, reinstate c0 as 0 if sstime includes a sample at t=0
      inndf$time1 <- NULL
      aucinndf <- inndf
      inndf[inndf<=blq] <- blq
      aucinndf[aucinndf<=blq] <- NA
      gendf$time2 <- NULL
      aucgendf <- gendf
      gendf[gendf<=blq] <- blq
      aucgendf[aucgendf<=blq] <- NA
      if(sstime[1]==0) {
        inndf[1,] <- 0
        aucinndf[1,] <- 0
        gendf[1,] <- 0
        aucgendf[1,] <- 0
      }
	  }else{ #if not NCA
	    time1 <- c(formdf1$TIME)
      id1 <- c(formdf1$ID)
      if(mode==1) {ipred1 <- c(formdf1$IPRED)}
	    if(mode>=3) {ipred1 <- c(formdf1$CP)}
      inndf <- data.frame(id1, time1, ipred1)
      inndf <- reshape(inndf,idvar="time1",timevar="id1",direction="wide")

	    time2 <- c(formdf2$TIME)
      id2 <- c(formdf2$ID)
      if(mode==1) {ipred2 <- c(formdf2$IPRED)}
	    if(mode>=3) {ipred2 <- c(formdf2$CP)}
      gendf <- data.frame(id2, time2, ipred2)
      gendf <- reshape(gendf,idvar="time2",timevar="id2",direction="wide")

	    inndf$time1 <- NULL
      gendf$time2 <- NULL
	  } #end if(mode==NCA)
  # Calculate cmax for each ID set, columns named dv.i where i=ID number
    inncmax <- as.numeric(apply(inndf, 2, max))
    gencmax <- as.numeric(apply(gendf, 2, max))
	  cratio <- gencmax/inncmax

  # Calculate tmax for each ID set
    inntmax <- as.numeric(apply(inndf, 2, tmax, time=sstime))
	  gentmax <- as.numeric(apply(gendf, 2, tmax, time=sstime))

  # Create result vector
    uid_unq <- unique(subdf$UID)
    uid <- unique(subdf$ID)
    cl <- subdf$CL[(1:nid)*nobs*2]
	  q <- subdf$Q[(1:nid)*nobs*2]
    v2 <- subdf$V2[(1:nid)*nobs*2]
	  v3 <- subdf$V3[(1:nid)*nobs*2]
    ka <- subdf$KA[(1:nid)*nobs*2]

	  if(mode==1) {
	    frel <- subdf$F1[(1:nid)*nobs*2]
      tempresult <- data.frame(uid_unq,uid,i,cl,q,v2,v3,ka,inncmax,gencmax,inntmax,gentmax,frel,cratio)
	  }
    if(mode==2) {
	# Calculate AUC(0-inf)
      innauc <- as.numeric(apply(aucinndf, 2, linlogAUCfunc, time.na=sstime, loq=blq))
      genauc <- as.numeric(apply(aucgendf, 2, linlogAUCfunc, time.na=sstime, loq=blq))
	    frel <- ifelse(innauc!=0, genauc/innauc, 0)

	    tempresult <- data.frame(uid_unq,uid,i,cl,q,v2,v3,ka,inncmax,gencmax,inntmax,gentmax,innauc,genauc,frel,cratio)
	  }
	  if(mode>=3) {
	    cl1 <- subdf$CL[(1:nid)*nobs*2-nobs]
	    q1 <- subdf$Q[(1:nid)*nobs*2-nobs]
	    v2.1 <- subdf$V2[(1:nid)*nobs*2-nobs]
	    v3.1 <- subdf$V3[(1:nid)*nobs*2-nobs]
	    frel <- subdf$F1[(1:nid)*nobs*2]
	    innauc <- subdf$AUC[(1:nid)*nobs*2-nobs]
	    genauc <- subdf$AUC[(1:nid)*nobs*2]
	    phfrel <- genauc/innauc
	    tempresult <- data.frame(uid_unq,uid,i,cl1,cl,q1,q,v2.1,v2,v3.1,v3,ka,inncmax,gencmax,inntmax,gentmax,innauc,genauc,frel,phfrel,cratio)
	  }
	  result[(1:nid)+nid*(i-1),] <- tempresult

  # Create average frel dataframe
    sumstatsfrel <- as.data.frame(geomeansemCI(frel))
    tempresult2<-data.frame(i,sumstatsfrel[1,],sumstatsfrel[2,],sumstatsfrel[3,])
	  fresult[i,] <- tempresult2
    sumstatscrat <- as.data.frame(geomeansemCI(cratio))
    tempresult3<-data.frame(i,sumstatscrat[1,],sumstatscrat[2,],sumstatscrat[3,])
	  cresult[i,] <- tempresult3
    if(mode>=3) {
      sumstatsphfrel <- as.data.frame(geomeansemCI(phfrel))
      tempresult4<-data.frame(i,sumstatsphfrel[1,],sumstatsphfrel[2,],sumstatsphfrel[3,])
      phfresult[i,] <- tempresult4
    }
  }
  write.csv(result, file=paste(SIM.file,file.tag,"result.csv", sep=""), row.names=FALSE)
  write.csv(fresult, file=paste(SIM.file,file.tag,"Fresult.csv", sep=""), row.names=FALSE)
  write.csv(cresult, file=paste(SIM.file,file.tag,"Cresult.csv", sep=""), row.names=FALSE)
  if(mode>=3) {write.csv(phfresult, file=paste(SIM.file,file.tag,"PHFresult.csv", sep=""), row.names=FALSE)}
}

#--------------------------------------------------------------------------------------------
# NON-COMPARTMENTAL FUNCTIONS
# AUC Function with Linear Up/Logarithmic Down Trapezoidal Method
  linlogAUCfunc <- function(dv.na,time.na,loq) {
  # Find AUC using trapezoidal method
  # Define values to be chosen for AUC calculation and give base value for AUC0t
	  n1 <- 1
	  n2 <- 2
	  n3 <- 3
	  n4 <- 4
	  AUC0t <- 0
	  dvtimedf <- na.omit(data.frame(dv.na,time.na))
	  dv <- c(unlist(dvtimedf[1]))
	  time <- c(unlist(dvtimedf[2]))
	  numobs <- length(time)

    #start loop to find AUC0t
	  while(n1<numobs) {
	  # Define variables to be used in trapezoidal method
	    c1 <- dv[n1]
		  c2 <- dv[n2]
		  t1 <- time[n1]
		  t2 <- time[n2]
		# Find sum of AUC0-t1 and the new AUCt1-t2
		  if(c2>c1)     # if second data point is larger -> Linear Trapezoidal = (t2-t1)*(C1+C2)/2
		  {
		    AUCtemp <- (t2-t1)*(c1+c2)/2
		    AUC0t <- sum(AUC0t,AUCtemp)
		  }
		  if(c2<c1)     # if second data point is smaller -> Logarithmic Trapezoidal = (t2-t1)*(C2-C1)/ln(C2/C1)
		  {
		    AUCtemp <- (t2-t1)*(c2-c1)/log(c2/c1)
		    AUC0t <- sum(AUC0t,AUCtemp)
		  }
      # Define next values to be chosen for AUC calculation and LOQ searching
		  n1 <- n1+1
		  n2 <- n2+1
		  n3 <- n3+1
		  n4 <- n4+1
	  } #end AUC0t

	# Define R2 to allow comparison, truncate time and dv to remove BLQ values
    ntail <- 3
    bestR2 <- 0
    bestk <- 0
	  whichtime <- which(dv==max(dv)) #designate point of cmax
	  adjdv <- dv[(whichtime-1):numobs] #delete all values before cmax as these are not needed for terminal phase, -1 to allow inclusion of Cmax in regression if required (as this is what WinNonLin does oddly)
	  flag <- 0

  	# While number of values being used to determine the slope is less than the number of total values* do the following
    #  - add one to ntail
    #  - define tail values for dv and time
    #  - fit line to log(dv) and time
    #  - define k as the slope and the R2
    #  - if R2 is better than the best R2 so far replace best R2 with new R2 and also replace best k with the new k
    #subtraction of one as C0 may cause an error when logged
	  if(length(adjdv)<ntail)
	  {
	    AUCtinf <- 0
	  }else{
      while(ntail<(length(adjdv)) && flag!=1) {
		    dvtail   <- unlist(tail(adjdv,ntail))
        timetail <- tail(time,ntail)
        fittail  <- lm(log(dvtail) ~ timetail)

        k  <- -1*fittail$coefficients["timetail"]
        R2 <- as.numeric(summary(fittail)["r.squared"])
		    adjR2 <- 1-((1-R2)*(ntail-1)/(ntail-2)) #adjusted R2 as per WinNonLin user guide
        if(adjR2>(bestR2-0.0001) && k>0) { #if statement as per WinNonLin user guide with added precautions against -ve k vals
          bestR2 <- adjR2
          bestk  <- k
		      ntail <- ntail+1
        }else{
		      flag <- 1
		    }
	    }
	    if(bestk == 0) {
	      AUCtinf <- 0
	    }else{
	# Calculate AUC from final time to infinite and then add it to AUC from time zero to final time
        Clast   <- tail(dv,1)
        AUCtinf <- Clast/bestk
      }
	  }
    AUC0inf <- AUC0t+AUCtinf
  }
#--------------------------------------------------------------------------------------------
# DATA ANALYSIS FUNCTIONS

#Define a function for geometric mean and 90% CI of the sem
  geomeansemCI <- function(x, na.rm=F) {
  #Note x cannot be negative, zero
    logx <- log(x)
    logmean <- mean(logx)
    n <- length(x)
    logsem <- sd(logx)/sqrt(n)
  #Critical value of the t-distribution for two one-sided p=0.05
    critt <- qt(.95, df=(n-1))

    loglo95 <- logmean - critt*logsem
    loghi95 <- logmean + critt*logsem
    gmean <- exp(logmean)
    glo95 <- exp(loglo95)
    ghi95 <- exp(loghi95)
    result <- c("gmean"=gmean, "glo95"=glo95, "ghi95"=ghi95, "crit.t"=critt)
    result
  }

  # Finds mean and 90% CI
  sumfunc95 <- function(x)   # Called sumfunc95 as this is a quick fix
  {
    stat1 <-  mean(unlist(x), na.rm=T)
    stat2 <-  quantile(x, probs=0.05, na.rm=T, names=F)  #90%CI
    stat3 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    result <- c("mean"=stat1, "low95"=stat2, "hi95"=stat3)
    result
  }

# Percent below the limit of quantification (not including t=0)
  percent.blq <- function(dv,time,blq)
  {
    timedv <- data.frame("TIME"=time,"DV"=dv)
    totaldv <- ifelse(timedv$TIME==0,NA,timedv$DV)
    totaldv <- totaldv[!is.na(totaldv)]
    loqdv <- totaldv[totaldv>=blq]
    percent.bloq <- (1-length(loqdv)/length(totaldv))*100
	  percent.bloq
  }

# Tmax
  tmax <- function(dv, time) { # computes the time of Cmax
    cmax <- max(dv)
    tindex <- which(dv==cmax)
    tmax <- time[tindex]
    head(tmax, n=1)   #as there can be 2 or more equal Cmax's, choose the first
  }

# Finds 95% CI corrected for mean (this is a ratio, not usable with bioequivalence)
  glohipercent.func <- function(indata,method,variable) {
    indata[5] <- indata[3]/indata[2]*100
    indata[6] <- indata[4]/indata[2]*100
    colnames(indata)[5:6] <- c(paste(variable,"_GLO95_PERCENT",sep=""),paste(variable,"_GHI95_PERCENT",sep=""))
    indata <- data.frame("Method"=method,indata)
    indata
   }

# Finds mean and 95% CI for specific dataset
  confint.func <- function(indata,method,variable) {
    mean1 <- as.data.frame(sumfunc95(indata[[2]]))
	  lo95 <- as.data.frame(sumfunc95(indata[[3]]))
    hi95 <- as.data.frame(sumfunc95(indata[[4]]))
    VARIABLE <- c(paste(variable,"_LO95",sep=""),paste(variable,"_MEAN",sep=""),paste(variable,"_HI95",sep=""))
	  ANALYSIS <- method
	  CILO95 <- c(lo95[2,], mean1[2,], hi95[2,])
	  MEAN <- c(lo95[1,], mean1[1,], hi95[1,])
	  CIHI95 <- c(lo95[3,], mean1[3,], hi95[3,])
	  all <- data.frame(VARIABLE,ANALYSIS,CILO95,MEAN,CIHI95)
	  all
  }

  # Determine bioequivalence
  bioq.func <- function(outputdf,limitlo,limithi,ctl.name)
  {
    outputdf$PF_FREL <- ifelse(outputdf$FREL_GLO95 < limitlo | outputdf$FREL_GHI95 > limithi,1,0)
    probtable <- ddply(outputdf, .(SIM_ID), function(df) CalcProb(df$PF_FREL))
    probtable <- data.frame("Metric"="FREL","Data"=ctl.name, probtable)
  }

  crat.func <- function(outputdf,limitlo,limithi,ctl.name) {
    outputdf$PF_CMAX <- ifelse(outputdf$CMAX_GLO95 < limitlo | outputdf$CMAX_GHI95 > limithi,1,0)
    probtable <- ddply(outputdf, .(SIM_ID), function(df) CalcProb(df$PF_CMAX))
    probtable <- data.frame("Metric"="CRAT","Data"=ctl.name, probtable)
  }

### Assess Bioequivalence
 # Assign Pass/Fail flag to Confidence Intervals
   # 0 is CI within limits, 1 is CI outside limits
  CalcProb <- function(x) {
   # Probability of 0 for binary events coded as 0 and 1
    prob <- sum(x)/length(x)
    prob <- 1-prob
    c("p"=prob)  #p of zero
  }

# Checking for type1 and type2 error against reference data
  # 1 is positive, 0 is negative
	errortype.func <- function(ref,test) {  # type=1 -> Type1 Error  type=2 -> Type2 Error
	  T1error <- gsub(TRUE,1,ref>test)	 # REF>TEST <- Type1 Error
	  T1error <- gsub(FALSE,0,T1error)
	  T1error <- gsub("NA",0,T1error)
	  T2error <- gsub(TRUE,1,ref<test)    # REF<TEST <- Type2 Error
	  T2error <- gsub(FALSE,0,T2error)
	  T2error <- gsub("NA",0,T2error)
	  data.frame("pT1"=T1error,"pT2"=T2error)
  }
      # Finds true SIMID for bioqtables, of particular use in M3 where NONMEM can fail to create fit file
	  # Then uses errortype.func to find type1 and type2 error
  errortype.process <- function(test,ref,fitfail,nsim,tag) {
    nfail <- length(fitfail)
	  trueSIMID <- (1:nsim)[-c(fitfail)]
    if(fitfail[1]!=0) {
      test$SIM_ID <- trueSIMID
	    missingdf <- data.frame(Metric=rep("FREL",times=nfail),"Data"=rep(tag,times=nfail),"SIM_ID"=fitfail,"p"=rep("NA",times=nfail))
	    truetable <- rbind(test,missingdf)
	    truetable <- orderBy(~SIM_ID,truetable)
	    result <- truetable
    }else{
	    result <- test
	  }
	  final <- data.frame(result,errortype.func(ref$p,result$p))
	  final
  }

### NONMEM Preparation and Processing Functions
# Control Stream creating function
# Can be used for mass replacement of specific parts of a .txt file
# pat and repl are vectors that represent each pattern to be replaced and each replacement respectively
  ctl.rep <- function(pat,repl,ctl) {
    tempctl <- ctl
    if(length(pat)==length(repl)) {
      for(i in 1:length(pat)) {
	      tempctl <- gsub(pat[i],repl[i],tempctl)
	    }
	    tempctl
	  }else{
	    warning("length(pattern)!=length(replacement): Amount of values to be replaced is not equal to amount of values given")
	  }
  }

# NONMEM Batch File creating function (also directs NONMEM to the .csv file by changing "dataname" in the controls stream)
  nm.prep <- function(nsim,ncore,nid,EST.file,ctl.name,SIM.name.out,limdata,limobs,blq,ctlref,nmbat,amt) {
    for (i in 1:nsim) {
	    file.name <- paste(EST.file,ctl.name,"_model",i,".ctl",sep="")
	    data.name <- paste(SIM.name.out,"_model",i,".csv",sep="")
      tempctl <- sub("dataname",data.name,ctlref)
	    writeLines(tempctl,file.name)
  # Create .csv for NONMEM input for each SIM
      subdata <- subset(limdata,SIM==i)
	    SIM <- rep(i,times=2*limobs*nid)
      EVID <- ifelse(subdata$TIME==0,4,0)
      BLQ <- ifelse(subdata$DV<blq,1,0)
	    BLQ <- ifelse(subdata$TIME==0,0,BLQ)
      AMT <- ifelse(EVID==4,amt,".")
	    MDV <- ifelse(AMT==amt,1,0)
	    MDV <- ifelse(subdata$DV<blq,1,MDV)
      nlmedata <- data.frame(subdata$ID,subdata$TIME,AMT,EVID,subdata$FORM,subdata$DV,MDV,SIM,BLQ)
      colnames(nlmedata) <- c("#ID","TIME","AMT","EVID","FORM","DV","MDVX","SIM","BLQ")
      file.name <- paste(EST.file,"_model",i,".csv",sep="")
      write.csv(nlmedata,file.name,quote=FALSE,row.names=FALSE)
  # Create command lines to be split into seperate .bats
      tempbat <- paste("call nmgo ",SIM.name.out,ctl.name,"_model",i,".ctl",sep="")
	    nmbat[i] <- tempbat
    }
    nmbat
  }

# Process output from NONMEM
  nlme.fit <- function(SIM.name.out,FIT.dir,EST.dir,ctl.name,nsim,nid,limobs) {
    fit2sim <- data.frame(matrix(NA, nrow = nsub*2, ncol = 11))
    colnames(fit2sim) <- c("UNQ_ID","STUD_ID","SIM_ID","AMT","F1","CL","Q","V2","V3","KA","AUC")
	  ctl.num <- ifelse(ctl.name=="M1",1,3)
	  rnum <- 1
	  fail.fit <- 0
    for (i in 1:nsim) {
      fit.name.in <- paste(SIM.name.out,ctl.num,"_model",i,sep="")
	    fit.dir.out <- paste(FIT.dir,"/",ctl.name,"_fit",rnum,".csv",sep="")
	    fit.file <- paste(EST.dir,"/",fit.name.in,".nm7/",fit.name.in,".fit",sep="")
	    if(file.exists(fit.file)) {  # Due to M3 terminating ~0.5% of the time
	      fitdata <- read.table(file=fit.file, sep="", skip=1, header=T, na.strings=c("NA","***********","1.#INFE+00"))
	      fitdata <- cbind(rep(1:nid+nid*(i-1),each=limobs*2),fitdata)
	      colnames(fitdata)[1] <- "UNQ_ID"
	      write.csv(fitdata, file=fit.dir.out, row.names=FALSE)
	      fit.temp <- read.csv(paste(FIT.dir,"/",ctl.name,"_fit",rnum,".csv",sep=""))
	      fit.temp <- subset(fit.temp,fit.temp$TIME==0)
	      fit.temp <- fit.temp[-c(4,6,7,8,9,17,18,19,20,21)]
		    fit.temp[3] <- rep(rnum,times=nid*2)
	      fit2sim[(1:(nid*2))+nid*2*(rnum-1),] <- fit.temp
		    rnum <- rnum+1                                                              # Record how many sims were successful to ensure safe use of data with universal functions
	    }else{
	      fail.fit <- c(fail.fit,i)                                                   # Record rows that did not give a .fit file
	    }
    }
	  if(length(fail.fit)>1) {
	    fail.fit <- fail.fit[-1]                                                      # If no sims failed to give .fit file, output==0, if sims did fail to give .fit file, remove 0 from output
	  }
    write.csv(na.omit(fit2sim), file=paste(SIM.file,ctl.name,"NMTHETAS.csv", sep="_"), row.names=FALSE)
	  c(rnum-1,fail.fit)
  }

# Simulate for NONMEM parameters
# Should be used together, see NCA_v_NLME for details on its usage
  nlme.simtime <- function(tdf) {
   	minKA <- (log(2)/min(tdf$KA))*6
	  TIME <- seq(from=0,to=minKA,by=minKA/200)
  }

  nlme.sim <- function(thetadf,TIME,nid,nsim) {
	  nobs <- length(TIME)
	  c1 <- match("UNQ_ID",colnames(thetadf))
	  c2 <- match("SIM_ID",colnames(thetadf))
	  c3 <- match("AUC",colnames(thetadf))
	  col1 <- unlist(thetadf[c1])
	  col2 <- unlist(thetadf[c2])
	  col3 <- unlist(thetadf[c3])
	  form <- rep(1:2,each=nobs,times=nid*nsim)
	  uid <- rep(col1,each=nobs)
	  sim <- rep(col2,each=nobs)
	  auc <- rep(col3,each=nobs)
	  tdf <- thetadf[-c(c1,c2,c3)]
	  colnames(tdf)[1] <- "ID"
	  sim.tdf <- mdply(tdf,simulate.2comp.abs)
	  result <- data.frame("UID"=uid,"SIM"=sim,"FORM"=form,sim.tdf,"AUC"=auc)
  }
 #####################################################
 #####################################################
                  # 2016 Additions #
 #####################################################
 #####################################################

  runaov2 <- function(df, USE) {
  #debug
  #df <- EXP.data[EXP.data$REP==1,]
  #USE <- "AUC"
  #uselog <- T

  #BE aov for study replicate i
    df$metric <- df[,USE]
    result <- aov(log(metric)~FORM, data=df)
    result2 <- summary(result)

  #This function avoids changing the contrasts for the aov
    aovtable <- model.tables(result,"means", se=T)
    int <- confint(result, level=0.9) # FUNCTION TO CALCULATE CONFIDENCE INTERVAL OF LN DATA

  #Extracting aov results
    Xt <- aovtable$tables$FORM[2]
    Xr <- aovtable$tables$FORM[1]

    pointestimate <- exp(Xt-Xr)
    pointestimate

    CI90lovalue <- exp(int[2,1])
    CI90lovalue

    CI90hivalue <- exp(int[2,2])
    CI90hivalue

  #Bioequivalence test 1=bioequivalent
    bioflag <- 0
    if (CI90lovalue > 0.8 & CI90hivalue < 1.25) bioflag <- 1

    BEresultsi <- data.frame("Metric"=USE,"pointestimate"=pointestimate,"lowerCI"=CI90lovalue,"upperCI"=CI90hivalue,"BE"=bioflag)
  }

  errortype.process2 <- function(test,ref,fitfail,nsim,tag) {
    test2 <- data.frame(test[,1],test[,2],test[,6])
	  colnames(test2) <- c("SIM","Metric","BE")
    if(fitfail[1]!=0) {
	    nfail <- length(fitfail)
	    trueSIMID <- (1:nsim)[-c(fitfail)]
      test2[1] <- trueSIMID
	    missingdf <- data.frame("SIM"=fitfail,"Metric"=rep("AUC",times=nfail),"BE"=rep("NA",times=nfail))
	    truetable <- rbind(test2,missingdf)
	    truetable <- orderBy(~SIM,truetable)
	    result <- truetable
    }else{
	    result <- test2
	  }
	  final <- errortype.func2(ref$BE,result$BE)
	  final
  }

  errortype.func2 <- function(ref,test) {
	  T1error <- gsub(TRUE,1,ref<test)
	  T1error <- gsub(FALSE,0,T1error)
	  T1error <- gsub("NA",0,T1error)
	  T2error <- gsub(TRUE,1,ref>test)
	  T2error <- gsub(FALSE,0,T2error)
	  T2error <- gsub("NA",0,T2error)
	  data.frame("pT1"=T1error,"pT2"=T2error)
  }
