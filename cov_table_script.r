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
  file.name <- "cov_collated_bioq_table.csv"
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

  write.csv(df.final.suc,"cov_results_successonly.csv",row.names=F)
  write.csv(df.final.suc,"cov_results_all.csv",row.names=F)

# PLOT SET 1
  #Determine median, upper lower bounds of the difference between using
  #all runs and only using successfully minimised runs

  df.diff <- data.frame(
    M1F1.F.SPEC = df.final.suc$M1F1FSPEC - df.final.all$M1F1FSPEC,
    M1PH.F.SPEC = df.final.suc$M1PHFSPEC - df.final.all$M1PHFSPEC,
    M1F1.F.SENS = df.final.suc$M1F1FSENS - df.final.all$M1F1FSENS,
    M1PH.F.SENS = df.final.suc$M1PHFSENS - df.final.all$M1PHFSENS,
    M1F1.F.ACC = df.final.suc$M1F1FACC - df.final.all$M1F1FACC,
    M1PH.F.ACC = df.final.suc$M1PHFACC - df.final.all$M1PHFACC,
    M1.C.SPEC = df.final.suc$M1CSPEC - df.final.all$M1CSPEC,
    M1.C.SENS = df.final.suc$M1CSENS - df.final.all$M1CSENS,
    M1.C.ACC = df.final.suc$M1CACC - df.final.all$M1CACC)
  #colwise(summary)(df.diff)

  df.final.all.M1F1FSPEC <- df.final.all$M1F1FSPEC
  df.final.all.M1F1FSPEC[df.final.all.M1F1FSPEC == 0] <- 0.0000001
  df.ratio <- data.frame(
    M1F1.F.SPEC = df.final.suc$M1F1FSPEC/df.final.all.M1F1FSPEC,
    M1PH.F.SPEC = df.final.suc$M1PHFSPEC/df.final.all$M1PHFSPEC,
    M1F1.F.SENS = df.final.suc$M1F1FSENS/df.final.all$M1F1FSENS,
    M1PH.F.SENS = df.final.suc$M1PHFSENS/df.final.all$M1PHFSENS,
    M1F1.F.ACC = df.final.suc$M1F1FACC/df.final.all$M1F1FACC,
    M1PH.F.ACC = df.final.suc$M1PHFACC/df.final.all$M1PHFACC,
    M1.C.SPEC = df.final.suc$M1CSPEC/df.final.all$M1CSPEC,
    M1.C.SENS = df.final.suc$M1CSENS/df.final.all$M1CSENS,
    M1.C.ACC = df.final.suc$M1CACC/df.final.all$M1CACC)
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
  ggsave("cov_DiffPlot_F1.png", width=20, height=16, units=c("cm"))
  plotobj1b <- plotobj1 + scale_y_continuous("Change in Percentage\n", lim = c(-20, 20))
  plotobj1b
  ggsave("cov_DiffPlot_F1_lim.png", width=20, height=16, units=c("cm"))

  CI90lo <- function(x) quantile(x,probs = 0.05, na.rm = T)
	CI90hi <- function(x) quantile(x,probs = 0.95, na.rm = T)
  df.diff.stat <- ddply(df.diff.l, .(method, bioq, statf), function(x) {
    c(CI90lo(x$value), mean(x$value, na.rm = T), CI90hi(x$value))
  })
  names(df.diff.stat)[4:6] <- c("ci90lo", "mean", "ci90hi")

  titletext1 <- expression(atop("Change in Results (Accuracy, Sensitivity, Specificity)",
    atop("CI of the difference between use of covariance step passing runs and all runs for bioequivalence")))
  plotobj1a <- NULL
  plotobj1a <- ggplot(data = df.diff.stat[df.diff.stat$bioq == "F", ], aes(factor(method), mean))
  plotobj1a <- plotobj1a + ggtitle(titletext1)
  plotobj1a <- plotobj1a + geom_point()
  plotobj1a <- plotobj1a + geom_errorbar(aes(ymin = ci90lo, ymax = ci90hi), width = 0.5)
  plotobj1a <- plotobj1a + scale_x_discrete("\nMethod")
  plotobj1a <- plotobj1a + scale_y_continuous("Change in Percentage\n",
    labels = dollar_format(suffix = "%", prefix = ""), lim = c(-8, 8))
  plotobj1a <- plotobj1a + facet_wrap(~statf)
  plotobj1a
  ggsave("cov_CIDiffPlot_F1_lim.png", width=20, height=8, units=c("cm"))

  stat2 <- read.csv("dfdiffstat.csv")
  df.diff.stat2 <- stat2[stat2$bioq == "F" & (stat2$method == "M1F1" | stat2$method == "M1PH"),]

  plotobj1c <- plotobj1a + geom_point(data = df.diff.stat2, colour = "grey50")
  plotobj1c <- plotobj1c + geom_errorbar(aes(ymin = ci90lo, ymax = ci90hi), data = df.diff.stat2, width = 0.5, colour = "grey50")
  ggsave("compar_CIDiffPlot_F1.png", width=20, height=8, units=c("cm"))

#  plotobj2 <- NULL
#  plotobj2 <- ggplot(data = df.diff.l[df.diff.l$bioq == "C", ])
#  plotobj2 <- plotobj2 + ggtitle(titletext1)
#  plotobj2 <- plotobj2 + geom_boxplot(aes(factor(method), value))
#  plotobj2 <- plotobj2 + scale_x_discrete("\nMethod")
#  plotobj2 <- plotobj2 + scale_y_continuous("Change in Percentage\n")
#  plotobj2 <- plotobj2 + facet_wrap(~statf)
#  plotobj2
#  ggsave("DiffPlot_CM.png", width=20, height=16, units=c("cm"))
