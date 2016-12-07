# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

  source("functions_NCAvNLME_2016.r")

  library(plyr)
  library(ggplot2)
  library(reshape2)
  library(stringr)

# Choose the directory you wish to compile results from
  master.dir <- "E:/hscpw-df1/Data1/Jim Hughes/2016"
  setwd(master.dir)
  file.name <- "collated_error_data.csv"
  dir.r16 <- "01_MDV"
  df.r16 <- read.csv(paste(dir.r16, file.name, sep="/"))
  dir.r78 <- "04_SS"
  df.r78 <- read.csv(paste(dir.r78, file.name, sep="/"))
  dir.r910 <- "06_genKA"
  df.r910 <- read.csv(paste(dir.r910, file.name, sep="/"))
  dir.r1114 <- "07_noBOV"
  df.r1114 <- read.csv(paste(dir.r1114, file.name, sep="/"))
  df.all <- rbind(df.r16, df.r78, df.r910, df.r1114)

  df.all$Scen <- df.all$SCEN
  df.all$Scen[df.all$RUN %% 2 == 0] <- df.all$SCEN[df.all$RUN %% 2 == 0] + 9
  df.all$Run <- ceiling(df.all$RUN/2)
  meth.str1 <- str_extract(df.all$Model, "Run[0-9]_Scen[0-9]{2}")
  meth.str2 <- str_extract(df.all$Model, "Run[0-9]{2}_Scen[0-9]{2}")
  meth.str <- c(meth.str1[!is.na(meth.str1)], meth.str2[!is.na(meth.str2)])
  df.all$Method <- as.numeric(substr(meth.str, nchar(meth.str), nchar(meth.str)))

  df.m1 <- df.all[df.all$Method == 1,]
  df.m3 <- df.all[df.all$Method == 3,]

  search.termstat <- as.character(unique(df.all$TermStat))
  search.covcode <- as.character(unique(df.all$CovCode))
  search.errcode <- as.character(unique(df.all$ErrCode))
  search.termdesc <- as.character(unique(df.all$TermDesc))

  per.onelevel <- function(x, col, search) {
    unlist(llply(1:length(search), function(i, x, col, search) {
      length(x[col][x[col] == search[i]])/dim(x)[1]*100
    }, x = x, col = col, search = search))
  }

  per.twolevel <- function(x, col, s.grid) {
    unlist(llply(1:dim(s.grid)[1], function(i, x, col, search) {
      length(x[col[1]][x[col[1]] == search[i,1] &
        x[col[2]] == search[i,2]])/dim(x)[1]*100
    }, x = x, col = col, search = s.grid))
  }

  #Percentages of successful vs. unsuccessful minimisation
  term.all.all <- per.onelevel(df.all, "TermStat", search.termstat)
  term.all.m1 <- per.onelevel(df.m1, "TermStat", search.termstat)
  term.all.m3 <- per.onelevel(df.m3, "TermStat", search.termstat)
  #Run
  term.run.all <- ddply(df.all, .(Run), function(x) per.onelevel(x, "TermStat", search.termstat))
  term.run.m1 <- ddply(df.m1, .(Run), function(x) per.onelevel(x, "TermStat", search.termstat))
  term.run.m3 <- ddply(df.m3, .(Run), function(x) per.onelevel(x, "TermStat", search.termstat))
  #Scen
  term.scen.all <- ddply(df.all, .(Scen), function(x) per.onelevel(x, "TermStat", search.termstat))
  term.scen.m1 <- ddply(df.m1, .(Scen), function(x) per.onelevel(x, "TermStat", search.termstat))
  term.scen.m3 <- ddply(df.m3, .(Scen), function(x) per.onelevel(x, "TermStat", search.termstat))
  #Run+Scen
  term.runscen.all <- ddply(df.all, .(Run, Scen), function(x) per.onelevel(x, "TermStat", search.termstat))
  term.runscen.m1 <- ddply(df.m1, .(Run, Scen), function(x) per.onelevel(x, "TermStat", search.termstat))
  term.runscen.m3 <- ddply(df.m3, .(Run, Scen), function(x) per.onelevel(x, "TermStat", search.termstat))

  #Percentages of successful vs. unsuccessful covariance step
  cov.all.all <- per.onelevel(df.all, "CovCode", search.covcode)
  cov.all.m1 <- per.onelevel(df.m1, "CovCode", search.covcode)
  cov.all.m3 <- per.onelevel(df.m3, "CovCode", search.covcode)
  #Run
  cov.run.all <- ddply(df.all, .(Run), function(x) per.onelevel(x, "CovCode", search.covcode))
  cov.run.m1 <- ddply(df.m1, .(Run), function(x) per.onelevel(x, "CovCode", search.covcode))
  cov.run.m3 <- ddply(df.m3, .(Run), function(x) per.onelevel(x, "CovCode", search.covcode))
  #Scen
  cov.scen.all <- ddply(df.all, .(Scen), function(x) per.onelevel(x, "CovCode", search.covcode))
  cov.scen.m1 <- ddply(df.m1, .(Scen), function(x) per.onelevel(x, "CovCode", search.covcode))
  cov.scen.m3 <- ddply(df.m3, .(Scen), function(x) per.onelevel(x, "CovCode", search.covcode))
  #Run+Scen
  cov.runscen.all <- ddply(df.all, .(Run, Scen), function(x) per.onelevel(x, "CovCode", search.covcode))
  cov.runscen.m1 <- ddply(df.m1, .(Run, Scen), function(x) per.onelevel(x, "CovCode", search.covcode))
  cov.runscen.m3 <- ddply(df.m3, .(Run, Scen), function(x) per.onelevel(x, "CovCode", search.covcode))

  #Percentages of different errors - "None", "Rounding errors", "Obj close to infinity", "Zero gradient", "Reached max evaluations"
  err.all.all <- per.onelevel(df.all, "ErrCode", search.errcode)
  err.all.m1 <- per.onelevel(df.m1, "ErrCode", search.errcode)
  err.all.m3 <- per.onelevel(df.m3, "ErrCode", search.errcode)
  #Run
  err.run.all <- ddply(df.all, .(Run), function(x) per.onelevel(x, "ErrCode", search.errcode))
  err.run.m1 <- ddply(df.m1, .(Run), function(x) per.onelevel(x, "ErrCode", search.errcode))
  err.run.m3 <- ddply(df.m3, .(Run), function(x) per.onelevel(x, "ErrCode", search.errcode))
  #Scen
  err.scen.all <- ddply(df.all, .(Scen), function(x) per.onelevel(x, "ErrCode", search.errcode))
  err.scen.m1 <- ddply(df.m1, .(Scen), function(x) per.onelevel(x, "ErrCode", search.errcode))
  err.scen.m3 <- ddply(df.m3, .(Scen), function(x) per.onelevel(x, "ErrCode", search.errcode))
  #Run+Scen
  err.runscen.all <- ddply(df.all, .(Run, Scen), function(x) per.onelevel(x, "ErrCode", search.errcode))
  err.runscen.m1 <- ddply(df.m1, .(Run, Scen), function(x) per.onelevel(x, "ErrCode", search.errcode))
  err.runscen.m3 <- ddply(df.m3, .(Run, Scen), function(x) per.onelevel(x, "ErrCode", search.errcode))

  #Percentages of different descriptions - "None", "Parameter estimate near boundary",
  #"R matrix algorithmically singular", "Problems occurred during minimisation",
  #"S matrix unobtainable", "Due to last iteration",
  #"Problems with individual", "Due to next iteration"
  desc.all.all <- per.onelevel(df.all, "TermDesc", search.termdesc)
  desc.all.m1 <- per.onelevel(df.m1, "TermDesc", search.termdesc)
  desc.all.m3 <- per.onelevel(df.m3, "TermDesc", search.termdesc)
  #Run
  desc.run.all <- ddply(df.all, .(Run), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  desc.run.m1 <- ddply(df.m1, .(Run), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  desc.run.m3 <- ddply(df.m3, .(Run), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  #Scen
  desc.scen.all <- ddply(df.all, .(Scen), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  desc.scen.m1 <- ddply(df.m1, .(Scen), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  desc.scen.m3 <- ddply(df.m3, .(Scen), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  #Run+Scen
  desc.runscen.all <- ddply(df.all, .(Run, Scen), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  desc.runscen.m1 <- ddply(df.m1, .(Run, Scen), function(x) per.onelevel(x, "TermDesc", search.termdesc))
  desc.runscen.m3 <- ddply(df.m3, .(Run, Scen), function(x) per.onelevel(x, "TermDesc", search.termdesc))

  qplot(factor(TermStat), data = df.all, geom = "bar", facets = ~Run)
  qplot(factor(CovCode), data = df.all, geom = "bar", facets = ~Run)
  qplot(factor(ErrCode), data = df.all, geom = "bar", facets = ~Run)
  qplot(factor(TermDesc), data = df.all, geom = "bar", facets = ~Run)
  qplot(factor(TermStat), data = df.all, geom = "bar", facets = ~Method)
  qplot(factor(CovCode), data = df.all, geom = "bar", facets = ~Method)
  qplot(factor(ErrCode), data = df.all, geom = "bar", facets = ~Method)
  qplot(factor(TermDesc), data = df.all, geom = "bar", facets = ~Method)

  err.all.m1[2]/term.all.m1[2]
  err.all.m3[2]/term.all.m3[2]

  err.run.m1[3]/term.run.m1[3]
  err.run.m3[3]/term.run.m3[3]

  err.scen.m1[3]/term.scen.m1[3]
  err.scen.m3[3]/term.scen.m3[3]

  err.runscen.m1[4]/term.runscen.m1[4]
  err.runscen.m3[4]/term.runscen.m3[4]


  desc.all.m1[2]/(cov.all.m1[2] - term.all.m1[2])*100
  desc.all.m3[2]/(cov.all.m3[2] - term.all.m3[2])*100

  desc.all.m1[3]/(cov.all.m1[2] - term.all.m1[2])*100
  desc.all.m3[3]/(cov.all.m3[2] - term.all.m3[2])*100




  #Determine search terms for
  s.grid.errdesc.all <- expand.grid(search.covcode, search.termdesc, stringsAsFactors = F)
  err.desc.all.all <- per.twolevel(df.all, c("CovCode","TermDesc"), s.grid.errdesc.all)
  s.grid.errdesc <- s.grid.errdesc.all[err.desc.all.all != 0, ]
  cols.errdesc <- paste(s.grid.errdesc$Var1, s.grid.errdesc$Var2)

  #Run
  errdesc.run.all <- ddply(df.all, .(Run), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  errdesc.run.m1 <- ddply(df.m1, .(Run), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  errdesc.run.m3 <- ddply(df.m3, .(Run), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  #Scen
  errdesc.scen.all <- ddply(df.all, .(Scen), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  errdesc.scen.m1 <- ddply(df.m1, .(Scen), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  errdesc.scen.m3 <- ddply(df.m3, .(Scen), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  #Run+Scen
  errdesc.runscen.all <- ddply(df.all, .(Run, Scen), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  errdesc.runscen.m1 <- ddply(df.m1, .(Run, Scen), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
  errdesc.runscen.m3 <- ddply(df.m3, .(Run, Scen), function(x) per.twolevel(x, c("ErrCode","TermDesc"), s.grid.errdesc))
