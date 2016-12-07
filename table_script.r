# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

  source("functions_NCAvNLME_2016.r")

  library(plyr)
  library(ggplot2)
  library(reshape2)

# Choose the directory you wish to compile results from
  master.dir <- "E:/hscpw-df1/Data1/Jim Hughes/2016"
  setwd(master.dir)
  file.name <- "collated_bioq_table.csv"
  dir.r16 <- "01_MDV"
  df.r16 <- read.csv(paste(dir.r16, file.name, sep="/"))
  dir.r78 <- "04_SS"
  df.r78 <- read.csv(paste(dir.r78, file.name, sep="/"))
  dir.r910 <- "06_genKA"
  df.r910 <- read.csv(paste(dir.r910, file.name, sep="/"))
  dir.r1114 <- "07_noBOV"
  df.r1114 <- read.csv(paste(dir.r1114, file.name, sep="/"))
  df.all <- rbind(df.r16, df.r78, df.r910, df.r1114)
  df.sub1 <- df.all[1:29]
  df.sub2 <- df.all[30:43]

  shk.file.name <- "collatedterm_shrinkage_data.csv"
  shk.r16 <- read.csv(paste(dir.r16, shk.file.name, sep="/"))
  shk.r78 <- read.csv(paste(dir.r78, shk.file.name, sep="/"))
  shk.r910 <- read.csv(paste(dir.r910, shk.file.name, sep="/"))
  shk.r1114 <- read.csv(paste(dir.r1114, shk.file.name, sep="/"))
  names(shk.r1114)[9:10] <- c("BOVV21", "BOVV22")
  shk.all <- rbind(shk.r16, shk.r78, shk.r910, shk.r1114)

  nsim <- 500
  nbioq <- df.all$IPREDBE*nsim/100
  nbioq.m1 <- df.all$M1IPREDBE*df.all$M1NSIM/100
  nbioq.m3 <- df.all$M3IPREDBE*df.all$M3NSIM/100
  nnonb <- nsim - nbioq
  nnonb.m1 <- df.all$M1NSIM - nbioq.m1
  nnonb.m3 <- df.all$M3NSIM - nbioq.m3

  cnbioq <- df.all$IPREDCM*nsim/100
  cnbioq.m1 <- df.all$M1IPREDCM*df.all$M1NSIM/100
  cnbioq.m3 <- df.all$M3IPREDCM*df.all$M3NSIM/100
  cnnonb <- nsim - cnbioq
  cnnonb.m1 <- df.all$M1NSIM - cnbioq.m1
  cnnonb.m3 <- df.all$M3NSIM - cnbioq.m3

  ncaTP <- nbioq - df.all$NCAFT2*nsim/100
  ncaTN <- nnonb - df.all$NCAFT1*nsim/100
  cncaTP <- cnbioq - df.all$NCACT2*nsim/100
  cncaTN <- cnnonb - df.all$NCACT1*nsim/100

  m1f1TP <- nbioq.m1 - df.all$M1F1FT2*df.all$M1NSIM/100
  m1f1TN <- nnonb.m1 - df.all$M1F1FT1*df.all$M1NSIM/100
  cm1TP <- cnbioq.m1 - df.all$M1CT2*df.all$M1NSIM/100
  cm1TN <- cnnonb.m1 - df.all$M1CT1*df.all$M1NSIM/100

  m1phTP <- nbioq.m1 - df.all$M1PHFT2*df.all$M1NSIM/100
  m1phTN <- nnonb.m1 - df.all$M1PHFT1*df.all$M1NSIM/100

  m3f1TP <- nbioq.m3 - df.all$M3F1FT2*df.all$M3NSIM/100
  m3f1TN <- nnonb.m3 - df.all$M3F1FT1*df.all$M3NSIM/100
  cm3TP <- cnbioq.m3 - df.all$M3CT2*df.all$M3NSIM/100
  cm3TN <- cnnonb.m3 - df.all$M3CT1*df.all$M3NSIM/100

  m3phTP <- nbioq.m3 - df.all$M3PHFT2*df.all$M3NSIM/100
  m3phTN <- nnonb.m3 - df.all$M3PHFT1*df.all$M3NSIM/100

  df.sub1$NCAFSENS <- ncaTP/nbioq*100
  df.sub1$M1F1FSENS <- m1f1TP/nbioq.m1*100
  df.sub1$M1PHFSENS <- m1phTP/nbioq.m1*100
  df.sub1$M3F1FSENS <- m3f1TP/nbioq.m3*100
  df.sub1$M3PHFSENS <- m3phTP/nbioq.m3*100

  df.sub1$NCAFSPEC <- ncaTN/nnonb*100
  df.sub1$M1F1FSPEC <- m1f1TN/nnonb.m1*100
  df.sub1$M1PHFSPEC <- m1phTN/nnonb.m1*100
  df.sub1$M3F1FSPEC <- m3f1TN/nnonb.m3*100
  df.sub1$M3PHFSPEC <- m3phTN/nnonb.m3*100

  df.sub1$NCAFACC <- (ncaTP+ncaTN)/nsim*100
  df.sub1$M1F1FACC <- (m1f1TP+m1f1TN)/df.all$M1NSIM*100
  df.sub1$M1PHFACC <- (m1phTP+m1phTN)/df.all$M1NSIM*100
  df.sub1$M3F1FACC <- (m3f1TP+m3f1TN)/df.all$M3NSIM*100
  df.sub1$M3PHFACC <- (m3phTP+m3phTN)/df.all$M3NSIM*100

  df.sub1$NCACSENS <- cncaTP/nbioq*100
  df.sub1$M1CSENS <- cm1TP/nbioq.m1*100
  df.sub1$M3CSENS <- cm3TP/nbioq.m3*100

  df.sub1$NCACSPEC <- cncaTN/nnonb*100
  df.sub1$M1CSPEC <- cm1TN/nnonb.m1*100
  df.sub1$M3CSPEC <- cm3TN/nnonb.m3*100

  df.sub1$NCACACC <- (cncaTP+cncaTN)/nsim*100
  df.sub1$M1CACC <- (cm1TP+cm1TN)/df.all$M1NSIM*100
  df.sub1$M3CACC <- (cm3TP+cm3TN)/df.all$M3NSIM*100

  df.final <- data.frame(df.sub1, df.sub2)

  df.final.all <- df.final[df.final$TERMSTAT == "All", ]
  df.final.suc <- df.final[df.final$TERMSTAT == "Only Success", ]

  write.csv(df.final.suc,"results_successonly.csv",row.names=F)
  write.csv(df.final.suc,"results_all.csv",row.names=F)

# PLOT SET 1
  #Determine median, upper lower bounds of the difference between using
  #all runs and only using successfully minimised runs

  df.diff <- data.frame(
    M1F1.F.SPEC = df.final.suc$M1F1FSPEC - df.final.all$M1F1FSPEC,
    M1PH.F.SPEC = df.final.suc$M1PHFSPEC - df.final.all$M1PHFSPEC,
    M3F1.F.SPEC = df.final.suc$M3F1FSPEC - df.final.all$M3F1FSPEC,
    M3PH.F.SPEC = df.final.suc$M3PHFSPEC - df.final.all$M3PHFSPEC,
    M1F1.F.SENS = df.final.suc$M1F1FSENS - df.final.all$M1F1FSENS,
    M1PH.F.SENS = df.final.suc$M1PHFSENS - df.final.all$M1PHFSENS,
    M3F1.F.SENS = df.final.suc$M3F1FSENS - df.final.all$M3F1FSENS,
    M3PH.F.SENS = df.final.suc$M3PHFSENS - df.final.all$M3PHFSENS,
    M1F1.F.ACC = df.final.suc$M1F1FACC - df.final.all$M1F1FACC,
    M1PH.F.ACC = df.final.suc$M1PHFACC - df.final.all$M1PHFACC,
    M3F1.F.ACC = df.final.suc$M3F1FACC - df.final.all$M3F1FACC,
    M3PH.F.ACC = df.final.suc$M3PHFACC - df.final.all$M3PHFACC,
    M1.C.SPEC = df.final.suc$M1CSPEC - df.final.all$M1CSPEC,
    M3.C.SPEC = df.final.suc$M3CSPEC - df.final.all$M3CSPEC,
    M1.C.SENS = df.final.suc$M1CSENS - df.final.all$M1CSENS,
    M3.C.SENS = df.final.suc$M3CSENS - df.final.all$M3CSENS,
    M1.C.ACC = df.final.suc$M1CACC - df.final.all$M1CACC,
    M3.C.ACC = df.final.suc$M3CACC - df.final.all$M3CACC)
  #colwise(summary)(df.diff)

  df.final.all.M1F1FSPEC <- df.final.all$M1F1FSPEC
  df.final.all.M1F1FSPEC[df.final.all.M1F1FSPEC == 0] <- 0.0000001
  df.ratio <- data.frame(
    M1F1.F.SPEC = df.final.suc$M1F1FSPEC/df.final.all.M1F1FSPEC,
    M1PH.F.SPEC = df.final.suc$M1PHFSPEC/df.final.all$M1PHFSPEC,
    M3F1.F.SPEC = df.final.suc$M3F1FSPEC/df.final.all$M3F1FSPEC,
    M3PH.F.SPEC = df.final.suc$M3PHFSPEC/df.final.all$M3PHFSPEC,
    M1F1.F.SENS = df.final.suc$M1F1FSENS/df.final.all$M1F1FSENS,
    M1PH.F.SENS = df.final.suc$M1PHFSENS/df.final.all$M1PHFSENS,
    M3F1.F.SENS = df.final.suc$M3F1FSENS/df.final.all$M3F1FSENS,
    M3PH.F.SENS = df.final.suc$M3PHFSENS/df.final.all$M3PHFSENS,
    M1F1.F.ACC = df.final.suc$M1F1FACC/df.final.all$M1F1FACC,
    M1PH.F.ACC = df.final.suc$M1PHFACC/df.final.all$M1PHFACC,
    M3F1.F.ACC = df.final.suc$M3F1FACC/df.final.all$M3F1FACC,
    M3PH.F.ACC = df.final.suc$M3PHFACC/df.final.all$M3PHFACC,
    M1.C.SPEC = df.final.suc$M1CSPEC/df.final.all$M1CSPEC,
    M3.C.SPEC = df.final.suc$M3CSPEC/df.final.all$M3CSPEC,
    M1.C.SENS = df.final.suc$M1CSENS/df.final.all$M1CSENS,
    M3.C.SENS = df.final.suc$M3CSENS/df.final.all$M3CSENS,
    M1.C.ACC = df.final.suc$M1CACC/df.final.all$M1CACC,
    M3.C.ACC = df.final.suc$M3CACC/df.final.all$M3CACC)
  #colwise(summary)(df.ratio)
  df.diff.l <- melt(df.diff)
  df.diff.l <- data.frame(colsplit(df.diff.l$variable, pattern = "\\.",
    names = c("method", "bioq", "stat")), value = df.diff.l$value)
  df.ratio.l <- melt(df.ratio)
  df.ratio.l <- data.frame(colsplit(df.ratio.l$variable, pattern = "\\.",
    names = c("method", "bioq", "stat")), value = df.ratio.l$value)

  df.diff.l$statf <- factor(df.diff.l$stat)
  levels(df.diff.l$statf) <- c("Accuracy", "Sensitivity", "Specificity")
  df.ratio.l$statf <- factor(df.ratio.l$stat)
  levels(df.ratio.l$statf) <- c("Accuracy", "Sensitivity", "Specificity")

  theme_set(theme_bw())

  titletext1 <- expression(atop("Change in Results (Accuracy, Sensitivity, Specificity)",
    atop("Difference between using only successfully minimised runs and using all runs to determine bioequivalence")))
  plotobj1 <- NULL
  plotobj1 <- ggplot(data = df.diff.l[df.diff.l$bioq == "F", ])
  plotobj1 <- plotobj1 + ggtitle(titletext1)
  plotobj1 <- plotobj1 + geom_boxplot(aes(factor(method), value))
  plotobj1 <- plotobj1 + scale_x_discrete("\nMethod")
  plotobj1 <- plotobj1 + scale_y_continuous("Change in Percentage\n")
  plotobj1 <- plotobj1 + facet_wrap(~statf)
  plotobj1
  ggsave("DiffPlot_F1.png", width=20, height=16, units=c("cm"))

  CI90lo <- function(x) quantile(x,probs = 0.05)
	CI90hi <- function(x) quantile(x,probs = 0.95)
  df.diff.stat <- ddply(df.diff.l, .(method, bioq, statf), function(x) {
    c(CI90lo(x$value), mean(x$value), CI90hi(x$value))
  })
  names(df.diff.stat)[4:6] <- c("ci90lo", "mean", "ci90hi")

  titletext1 <- expression(atop("Change in Results (Accuracy, Sensitivity, Specificity)",
    atop("Confidence intervals of difference between use of successfully minimised runs and all runs to determine bioequivalence")))
  plotobj1a <- NULL
  plotobj1a <- ggplot(data = df.diff.stat[df.diff.stat$bioq == "F", ], aes(factor(method), mean))
  plotobj1a <- plotobj1a + ggtitle(titletext1)
  plotobj1a <- plotobj1a + geom_point()
  plotobj1a <- plotobj1a + geom_errorbar(aes(ymin = ci90lo, ymax = ci90hi), width = 0.5)
  plotobj1a <- plotobj1a + scale_x_discrete("\nMethod")
  plotobj1a <- plotobj1a + scale_y_continuous("Change in Percentage\n", lim = c(-8, 8))
  plotobj1a <- plotobj1a + facet_wrap(~statf)
  plotobj1a
  ggsave("CIDiffPlot_F1_lim.png", width=20, height=8, units=c("cm"))

#  plotobj2 <- NULL
#  plotobj2 <- ggplot(data = df.diff.l[df.diff.l$bioq == "C", ])
#  plotobj2 <- plotobj2 + ggtitle(titletext1)
#  plotobj2 <- plotobj2 + geom_boxplot(aes(factor(method), value))
#  plotobj2 <- plotobj2 + scale_x_discrete("\nMethod")
#  plotobj2 <- plotobj2 + scale_y_continuous("Change in Percentage\n")
#  plotobj2 <- plotobj2 + facet_wrap(~statf)
#  plotobj2
#  ggsave("DiffPlot_CM.png", width=20, height=16, units=c("cm"))

# PLOT SET 2
  #Identify patterns that show change in shrinkages relating to Type 1 error
  #Do this by finding median, upper and lower bounds of shrinkages

  shk.all$Scen <- shk.all$SCEN
  shk.all$Scen[shk.all$RUN %% 2 == 0] <- shk.all$SCEN[shk.all$RUN %% 2 == 0] + 9
  shk.all$Run <- ceiling(shk.all$RUN/2)
  shk.all.l <- melt(shk.all[shk.all$Term == "Success", ], c("RUN","SCEN","Run","Scen","Sim","Method","Term"))
  shk.all.l$RunVar <- paste("Run",shk.all.l$Run,shk.all.l$variable)

  titletext3 <- expression(atop("Box Plot for Shrinkages on each ETA",
    atop("Shrinkages grouped by Simulation Run")))
  plotobj3 <- NULL
  plotobj3 <- ggplot(data = shk.all.l)
  plotobj3 <- plotobj3 + ggtitle(titletext3)
  plotobj3 <- plotobj3 + geom_boxplot(aes(factor(Run), value))
  plotobj3 <- plotobj3 + scale_x_discrete("\nRun")
  plotobj3 <- plotobj3 + scale_y_continuous("Shrinkage (%)\n", lim = c(0, 100))
  plotobj3 <- plotobj3 + facet_wrap(~variable)
  plotobj3
  ggsave("SHKplot_byRun.png", width=20, height=16, units=c("cm"))

  titletext4 <- expression(atop("Box Plot for Shrinkages on each ETA",
    atop("Shrinkages grouped by Scenario")))
  plotobj4 <- NULL
  plotobj4 <- ggplot(data = shk.all.l)
  plotobj4 <- plotobj4 + ggtitle(titletext4)
  plotobj4 <- plotobj4 + geom_boxplot(aes(factor(Scen), value))
  plotobj4 <- plotobj4 + scale_x_discrete("\nScenario")
  plotobj4 <- plotobj4 + scale_y_continuous("Shrinkage (%)\n", lim = c(0, 100))
  plotobj4 <- plotobj4 + facet_wrap(~variable)
  plotobj4
  ggsave("SHKplot_byScen.png", width=20, height=16, units=c("cm"))

  titletext5 <- expression(atop("Box Plot for Shrinkages on each ETA",
    atop("Shrinkages grouped by Run and Scenario")))
  plotobj5 <- NULL
  plotobj5 <- ggplot(data = shk.all.l)
  plotobj5 <- plotobj5 + ggtitle(titletext5)
  plotobj5 <- plotobj5 + geom_boxplot(aes(factor(Scen), value), outlier.size = 0.5)
  plotobj5 <- plotobj5 + scale_x_discrete("\nScenario")
  plotobj5 <- plotobj5 + scale_y_continuous("Shrinkage (%)\n", lim = c(0, 100))
  plotobj5 <- plotobj5 + facet_wrap(~RunVar)
  #plotobj5
  ggsave("SHKplot_byScen_facetRun.png", width=30, height=24, units=c("cm"))

  titletext6 <- expression(atop("Box Plot for Shrinkages on each ETA",
    atop("Shrinkages grouped by Run and Scenario")))
  plotobj6 <- NULL
  plotobj6 <- ggplot(data = shk.all.l[shk.all.l$RunVar == "", ])
  plotobj6 <- plotobj6 + ggtitle(titletext6)
  plotobj6 <- plotobj6 + geom_boxplot(aes(factor(Scen), value), outlier.size = 0.5)
  plotobj6 <- plotobj6 + scale_x_discrete("\nScenario")
  plotobj6 <- plotobj6 + scale_y_continuous("Shrinkage (%)\n", lim = c(0, 100))
  #plotobj6 <- plotobj6 + facet_wrap(~RunVar)
  #plotobj6
  ggsave("SHKplot_byScen_facetRun.png", width=30, height=24, units=c("cm"))

# PLOT SET 3
  #Identify patterns that show change in shrinkages relating to Type 1 error
  #Do this by plotting median shrinkages against type 1 error percentage

  shk.median <- dcast(
    ddply(na.omit(shk.all.l), .(Run, Scen, Method, variable),
      function(x) summary(x$value)["Median"]),
    Run+Scen+Method ~ variable
  )
  shk.median.d <- arrange(rbind(shk.median, shk.median), Run, Scen, Method)
  shk.median.d$Meth <- c("M1F1","M1PH","M3F1","M3PH")

  df.suc.sub <- df.final.suc[c(1:2, 15:18)]
  names(df.suc.sub)[3:6] <- c("M1F1","M1PH","M3F1","M3PH")
  df.suc.l <- arrange(melt(df.suc.sub, c("RUN", "SCEN")), RUN, SCEN, variable)

  shk.t1.df <- data.frame(
    shk.median.d[-3],
    T1 = df.suc.l$value
  )
  shk.t1.df.l <- melt(shk.t1.df, c("Run", "Scen", "Meth", "T1"))

  shk.t1.r2 <- ddply(shk.t1.df.l, .(Meth, variable), function(x) summary(lm(x$T1 ~ x$value))$r.squared)

  titletext7 <- expression(atop("Percent ETA Shrinkage versus Type 1 Error",
    atop("Using M1 and the \"F estimate\" method ")))
  plotobj7 <- NULL
  plotobj7 <- ggplot(data = shk.t1.df.l[shk.t1.df.l$Meth == "M1F1", ])
  plotobj7 <- plotobj7 + ggtitle(titletext7)
  plotobj7 <- plotobj7 + geom_point(aes(x = value, y = T1))
  plotobj7 <- plotobj7 + scale_x_continuous("Shrinkage (%)")
  plotobj7 <- plotobj7 + scale_y_continuous("Type 1 Error (%)")
  plotobj7 <- plotobj7 + facet_wrap(~ variable, nrow = 2, ncol = 4)
  plotobj7

  titletext8 <- expression(atop("Percent ETA Shrinkage versus Type 1 Error",
    atop("Using M1 and the \"Post-Hoc\" method ")))
  plotobj8 <- NULL
  plotobj8 <- ggplot(data = shk.t1.df.l[shk.t1.df.l$Meth == "M1PH", ])
  plotobj8 <- plotobj8 + ggtitle(titletext8)
  plotobj8 <- plotobj8 + geom_point(aes(x = value, y = T1))
  plotobj8 <- plotobj8 + scale_x_continuous("Shrinkage (%)")
  plotobj8 <- plotobj8 + scale_y_continuous("Type 1 Error (%)")
  plotobj8 <- plotobj8 + facet_wrap(~ variable, nrow = 2, ncol = 4)
  plotobj8

  titletext9 <- expression(atop("Percent ETA Shrinkage versus Type 1 Error",
    atop("Using M3 and the \"F estimate\" method ")))
  plotobj9 <- NULL
  plotobj9 <- ggplot(data = shk.t1.df.l[shk.t1.df.l$Meth == "M3F1", ])
  plotobj9 <- plotobj9 + ggtitle(titletext9)
  plotobj9 <- plotobj9 + geom_point(aes(x = value, y = T1))
  plotobj9 <- plotobj9 + scale_x_continuous("Shrinkage (%)")
  plotobj9 <- plotobj9 + scale_y_continuous("Type 1 Error (%)")
  plotobj9 <- plotobj9 + facet_wrap(~ variable, nrow = 2, ncol = 4)
  plotobj9

  titletext0 <- expression(atop("Percent ETA Shrinkage versus Type 1 Error",
    atop("Using M3 and the \"Post-Hoc\" method ")))
  plotobj0 <- NULL
  plotobj0 <- ggplot(data = shk.t1.df.l[shk.t1.df.l$Meth == "M3PH", ])
  plotobj0 <- plotobj0 + ggtitle(titletext0)
  plotobj0 <- plotobj0 + geom_point(aes(x = value, y = T1))
  plotobj0 <- plotobj0 + scale_x_continuous("Shrinkage (%)")
  plotobj0 <- plotobj0 + scale_y_continuous("Type 1 Error (%)")
  plotobj0 <- plotobj0 + facet_wrap(~ variable, nrow = 2, ncol = 4)
  plotobj0

  shk.grid.df <- shk.t1.df.l[shk.t1.df.l$variable %in% c("BSVCL", "BSVV2", "BSVKA", "BSVF1"), ]
  shk.grid.df$MethVar <- factor(paste(shk.grid.df$Meth, shk.grid.df$variable),
    levels = c(
      "M1F1 BSVCL", "M1F1 BSVV2", "M1PH BSVCL", "M1PH BSVV2",
      "M1F1 BSVKA", "M1F1 BSVF1", "M1PH BSVKA", "M1PH BSVF1",
      "M3F1 BSVCL", "M3F1 BSVV2", "M3PH BSVCL", "M3PH BSVV2",
      "M3F1 BSVKA", "M3F1 BSVF1", "M3PH BSVKA", "M3PH BSVF1")
    )

  titletext10 <- expression(atop("Percent BSV ETA Shrinkage versus Type 1 Error",
    atop("Split into different methods of determining bioequivalence")))
  plotobj10 <- NULL
  plotobj10 <- ggplot(data = shk.grid.df)
  plotobj10 <- plotobj10 + ggtitle(titletext10)
  plotobj10 <- plotobj10 + geom_point(aes(x = value, y = T1), size = 0.5)
  plotobj10 <- plotobj10 + scale_x_continuous("Shrinkage (%)")
  plotobj10 <- plotobj10 + scale_y_continuous("Type 1 Error (%)")
  plotobj10 <- plotobj10 + facet_wrap(~ MethVar, nrow = 4, ncol = 4)
  plotobj10
  ggsave("SHKT1plot_grid.png", width=20, height=16, units=c("cm"))
