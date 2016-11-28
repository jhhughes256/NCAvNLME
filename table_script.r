# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

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
  df.sub2 <- df.all[30:39]

  nsim <- 500
  nbioq <- df.all$IPREDBE*nsim/100
  nbioq.m1 <- df.all$IPREDBE*nsim/100
  nbioq.m3 <- df.all$IPREDBE*nsim/100
  nnonb <- nsim - nbioq
  nnonb.m1 <- nsim - nbioq.m1
  nnonb.m3 <- nsim - nbioq.m3

  cnbioq <-
  cnnonb <- nsim - nbioq

  ncaTP <- nbioq - df.all$NCAFT2*nsim/100
  ncaTN <- nnonb - df.all$NCAFT1*nsim/100
  cncaTP <- nbioq - df.all$NCACT2*nsim/100
  cncaTN <- nnonb - df.all$NCACT1*nsim/100

  m1f1TP <- nbioq - df.all$M1F1FT2*nsim/100
  m1f1TN <- nnonb - df.all$M1F1FT1*nsim/100
  cm1TP <- nbioq - df.all$M1CT2*nsim/100
  cm1TN <- nnonb - df.all$M1CT1*nsim/100

  m1phTP <- nbioq - df.all$M1PHFT2*nsim/100
  m1phTN <- nnonb - df.all$M1PHFT1*nsim/100

  m3f1TP <- nbioq - df.all$M3F1FT2*nsim/100
  m3f1TN <- nnonb - df.all$M3F1FT1*nsim/100
  cm3TP <- nbioq - df.all$M3CT2*nsim/100
  cm3TN <- nnonb - df.all$M3CT1*nsim/100

  m3phTP <- nbioq - df.all$M3PHFT2*nsim/100
  m3phTN <- nnonb - df.all$M3PHFT1*nsim/100

  df.sub1$NCAFSENS <- ncaTP/nbioq
  df.sub1$M1F1FSENS <- m1f1TP/nbioq
  df.sub1$M1PHFSENS <- m1phTP/nbioq
  df.sub1$M3F1FSENS <- m3f1TP/nbioq
  df.sub1$M3PHFSENS <- m3phTP/nbioq

  df.sub1$NCAFSPEC <- ncaTN/nnonb
  df.sub1$M1F1FSPEC <- m1f1TN/nnonb
  df.sub1$M1PHFSPEC <- m1phTN/nnonb
  df.sub1$M3F1FSPEC <- m3f1TN/nnonb
  df.sub1$M3PHFSPEC <- m3phTN/nnonb

  df.sub1$NCAFACC <- (ncaTP+ncaTN)/nsim
  df.sub1$M1F1FACC <- (m1f1TP+m1f1TN)/nsim
  df.sub1$M1PHFACC <- (m1phTP+m1phTN)/nsim
  df.sub1$M3F1FACC <- (m3f1TP+m3f1TN)/nsim
  df.sub1$M3PHFACC <- (m3phTP+m3phTN)/nsim

  df.sub1$NCACSENS <- cncaTP/nbioq
  df.sub1$M1CSENS <- cm1TP/nbioq
  df.sub1$M3CSENS <- cm3TP/nbioq

  df.sub1$NCACSPEC <- cncaTN/nnonb
  df.sub1$M1CSPEC <- cm1TN/nnonb
  df.sub1$M3CSPEC <- cm3TN/nnonb

  df.sub1$NCACACC <- (cncaTP+cncaTN)/nsim
  df.sub1$M1CACC <- (cm1TP+cm1TN)/nsim
  df.sub1$M3CACC <- (cm3TP+cm3TN)/nsim
