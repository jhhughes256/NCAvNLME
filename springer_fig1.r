### Remaking Figure 1 in R

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set directories
  source.dir <- "E:/hscpw-df1/Data1/Jim Hughes/2016/01_MDV"
  setwd(source.dir)
  source("functions_RU_ANOVA_2016.r")
  master.dir <- "E:/hscpw-df1/Data1/Jim Hughes/SpringerFigs"      ### Directory containing source files
  setwd(master.dir)

  library(ggplot2)
  theme_set(theme_bw()) 
  
  loq <- 0.14
  
  d1 <- c("ID"=0,"AMT"=270,"CL"=20,"Q"=30,"V2"=100,"V3"=400,"KA"=0.5,"F1"=1)
  d2 <- c("ID"=0,"AMT"=270,"CL"=25,"Q"=40,"V2"=100,"V3"=600,"KA"=0.4,"F1"=1)
  
  TIME <- c(0,0.25,0.5,1,2,4,6,8,12,16,20,24)
  
  df1 <- data.frame("ID" = 1, simulate.2comp.abs("ID"=d1[1],"AMT"=d1[2],"CL"=d1[3],"Q"=d1[4],"V2"=d1[5],"V3"=d1[6],"KA"=d1[7],"F1"=d1[8]))
  df2 <- data.frame("ID" = 2, simulate.2comp.abs("ID"=d2[1],"AMT"=d2[2],"CL"=d2[3],"Q"=d2[4],"V2"=d2[5],"V3"=d2[6],"KA"=d2[7],"F1"=d2[8]))
  dfn <- rbind(df1,df2)
  dfn[dfn$CP<=0.1, ] <- NA
  
  df1[df1$CP<=loq, ] <- NA
  df1 <- na.omit(df1)
  fit.df1 <- coef(lm(log(tail(df1$CP,3)) ~ tail(df1$TIME,3)))
  plot.df1 <- rbind(c(1,10,exp(fit.df1[2]*10+fit.df1[1])),tail(df1,3),c(1,24,exp(fit.df1[2]*24+fit.df1[1])))
  
  df2[df2$CP<=loq, ] <- NA
  df2 <- na.omit(df2)
  fit.df2 <- coef(lm(log(tail(df2$CP,3)) ~ tail(df2$TIME,3)))
  plot.df2 <- rbind(c(2,4,exp(fit.df2[2]*4+fit.df2[1])),tail(df2,3),c(2,13.3,exp(fit.df2[2]*13.3+fit.df2[1])))
  
# Set up rich sampling times
  T.1 <- seq(from=0, to=3, by=0.05)
  T.25 <- seq(from=3.25, to=6, by=0.25)
  T.5 <- seq(from=6.5, to=12, by=0.5)
  T1 <- seq(from=13, to=24, by=1)
  TIME <- c(T.1,T.25,T.5,T1)
  
  df3 <- data.frame("ID" = 1, simulate.2comp.abs("ID"=d1[1],"AMT"=d1[2],"CL"=d1[3],"Q"=d1[4],"V2"=d1[5],"V3"=d1[6],"KA"=d1[7],"F1"=d1[8]))
  df4 <- data.frame("ID" = 2, simulate.2comp.abs("ID"=d2[1],"AMT"=d2[2],"CL"=d2[3],"Q"=d2[4],"V2"=d2[5],"V3"=d2[6],"KA"=d2[7],"F1"=d2[8]))
  dfx <- rbind(df3, df4)
  dfx[dfx$CP<=0.1, ] <- NA
  
  plotobj <- NULL
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x=TIME, y=CP), dfn[dfn$ID==1,], size=2, colour="blue")
  plotobj <- plotobj + geom_line(aes(x=TIME, y=CP), dfx[dfx$ID==1,], size=1, colour="blue")
  plotobj <- plotobj + geom_smooth(aes(x=TIME, y=CP), plot.df1, method="lm", se=FALSE, linetype="longdash", colour="black")
  plotobj <- plotobj + geom_line(aes(x=TIME, y=CP), dfx[dfx$ID==2,], size=1, colour="red")
  plotobj <- plotobj + geom_point(aes(x=TIME, y=CP), dfn[dfn$ID==2,], size=2, colour="red")
  plotobj <- plotobj + geom_smooth(aes(x=TIME, y=CP), plot.df2, method="lm", se=FALSE, linetype="longdash", colour="black")
  plotobj <- plotobj + geom_hline(yintercept=loq, linetype="dashed")
  plotobj <- plotobj + scale_y_log10("Concentration (ug/L)\n", limits=c(0.1,1.1))
  plotobj <- plotobj + scale_x_continuous("\nTime (h)", breaks=c(0,4,8,12,16,20,24))
  plotobj
  ggsave("NCA_Fig1.pdf")
  ggsave("NCA_Fig1.eps")
  
  