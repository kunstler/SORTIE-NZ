#######################################
#######################################
#### REANALYSE ALL DATA UNDER 600 M OF ELEVATION

Dave.m.ss.nor <- read.csv("data.adult.csv")

###### calcule données BA et basup bainf
Dave.dec <- Dave.m.ss.nor[!is.na(Dave.m.ss.nor$PLHDAltitude) & Dave.m.ss.nor$PLHDAltitude<601,]
Dave.dec <- Dave.dec[Dave.dec$dbh1>100,]

Dave.mdd2<- Dave.dec[Dave.dec$DMTRSpeciesCode %in% c("DACCUP","METUMB","NOTCLI","NOTMEN","PODHAL","PRUFER","WEIRAC"),]      ### c("SCHDIG","FUCEXC")
Dave.mdd2$DMTRSpeciesCode  <- factor(Dave.mdd2$DMTRSpeciesCode)
Dave.md2 <- Dave.mdd2

vecbasup <- rep(0,length(Dave.md2$dbh1))
vecbainf <- rep(0,length(Dave.md2$dbh1))


basal.plot <- tapply(pi*(Dave.dec$dbh1/10)^2, INDEX= Dave.dec$pol.idi.y,FUN=sum)
basal.plot.10 <- tapply(pi*(Dave.dec$dbh1[Dave.dec$dbh1>100]/10)^2, INDEX= Dave.dec$pol.idi.y[Dave.dec$dbh1>100],FUN=sum)
data.basal.plot <- data.frame(pol.idi.y = names(basal.plot), BA.p= basal.plot/1000,BA.p.10 = basal.plot.10 )

Dave.md2T <- data.frame(Dave.md2,bas=vecbasup/1000,bai=vecbainf/1000)

Dave.md2T2 <-  merge(Dave.md2T, data.basal.plot,by="pol.idi.y")



Dave.dec.m <- Dave.mdd2[,]


####################
###### Analysis
####################

library( neighparam )




##data
surv <- rep(1,length(Dave.dec.m$DMTRSpeciesCode))
surv[Dave.dec.m$dbh2==0] <- 0
Dave.dec.s <- data.frame(Dave.dec.m,surv=surv,year=Dave.dec.m$year2- Dave.dec.m$year1)
split.ns <- split(Dave.dec.s,  Dave.dec.s$DMTRSpeciesCode)

length(split.ns[[1]][,1])
#################
##### ESTIMATION

## function
model.s1 <- function( max.s,x.0,x.b,diam2,year)
{
 pre <- (max.s*exp(-1/2*(log(diam2/x.0)/x.b)^2 ))^year
 return( pre)
}

###
list.Adult.S.1b <- vector("list",7)
for (i in 1:7)

{
print(i)
#windows()
data.temp <- split.ns[[i]][(split.ns[[i]]$dbh1/10)>10,]
surv<-( data.temp$surv)
diam2 <- data.temp$dbh1/10
year <- data.temp$year
datatemp <- data.frame(surv=surv,diam2=diam2,year=year)

par <- list(max.s = .9, x.0 = mean(datatemp$diam2),x.b=10, diam2 = "diam2",year="year")
par_lo <- list(max.s =0.000 , x.0 = 2 ,x.b =0.0001)
par_hi <- list(max.s = 1 , x.0 = 350 ,x.b= 80)
par_step <- list(max.s=1, x.0 = 3 ,x.b= 2)
par$x <- "surv"
## Mean in normal PDF
par$prob <- "predicted"
par$size <- 1
## Have it calculate log likelihood
par$log <- TRUE
list.Adult.S.1b[[i]] <- neighanneal (model = model.s1, par,
 source_data= datatemp, par_lo, par_hi, par_step, pdf= dbinom, dep_var="surv", max_iter= 50000)
title(main= names(split.ns)[i])
}

##############
## a longer estimation lead to different estimate .....
################"


  ## param
  Matrix.Param.S.mALL <- matrix(0,nrow=7,ncol=5)
  rownames(Matrix.Param.S.mALL) <- names(split.ns)
  colnames(Matrix.Param.S.mALL) <-  c("max.s","x.0","x.b","R²","num")

  for (i in 1:7)
  {
   Matrix.Param.S.mALL[i,1] <-  list.Adult.S.1[[i]]$best_pars[[1]]
   Matrix.Param.S.mALL[i,2] <-  list.Adult.S.1[[i]]$best_pars[[2]]
   Matrix.Param.S.mALL[i,3] <-  list.Adult.S.1[[i]]$best_pars[[3]]
   Matrix.Param.S.mALL[i,4] <-  list.Adult.S.1[[i]]$R2
   Matrix.Param.S.mALL[i,5] <-  length(split.ns[[i]][(split.ns[[i]]$dbh1/10)>10,1])
  }

  Matrix.Param.S.mALLb <- matrix(0,nrow=7,ncol=5)
  rownames(Matrix.Param.S.mALLb) <- names(split.ns)
  colnames(Matrix.Param.S.mALLb) <-  c("max.s","x.0","x.b","R²","num")

  for (i in 1:7)
  {
   Matrix.Param.S.mALLb[i,1] <-  list.Adult.S.1b[[i]]$best_pars[[1]]
   Matrix.Param.S.mALLb[i,2] <-  list.Adult.S.1b[[i]]$best_pars[[2]]
   Matrix.Param.S.mALLb[i,3] <-  list.Adult.S.1b[[i]]$best_pars[[3]]
   Matrix.Param.S.mALLb[i,4] <-  list.Adult.S.1b[[i]]$R2
   Matrix.Param.S.mALLb[i,5] <-  length(split.ns[[i]][(split.ns[[i]]$dbh1/10)>10,1])
  }


write.table( Matrix.Param.S.mALLb,file="param.s.all.NEW.txt")

Matrix.Param.S.mALLb <- read.table("param.s.all.NEW.txt")
                                       
#######################################
###### PLOTS

 windows()
   
   Dd <- 10:150
   i <- 1
    plot(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),xlab="DBH (cm)",ylab="Annual probability of mortality",ylim=c(0,0.03),type="l")

   i <- 3
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=2,lwd=2)
   i <- 4
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=3)
   i <- 5
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=4,lwd=2)
   i <- 6

  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=5,lwd=1)
   i <- 7
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=1,lwd=2)
      i <- 2
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=6,lwd=1)
  
  vec.lwd <- c(1,2,1,2,1,2,1)
  vec.lty <- c(1:5,1,6)
  
  legend(locator(1),legend=names(split.ns)[c( 1,3:7,2)],lty=vec.lty,lwd=vec.lwd,bty="n")
  
                                       
                                       
                                       
#######################################
#######################################
#######################################
#### for sapling survival
### data no perturbé pb

Dave.dec <- Dave.m.ss.nor[!is.na(Dave.m.ss.nor$PLHDAltitude) & Dave.m.ss.nor$PLHDAltitude<601,]


Dave.mdd2<- Dave.dec[Dave.dec$DMTRSpeciesCode %in% c("DACCUP","METUMB","NOTCLI","NOTMEN","PODHAL","PRUFER","WEIRAC"),]      ### c("SCHDIG","FUCEXC")
Dave.mdd2$DMTRSpeciesCode  <- factor(Dave.mdd2$DMTRSpeciesCode)
Dave.dec.m <- Dave.mdd2[,]


####################
###### Analysis
####################

library( neighparam )




##data
surv <- rep(1,length(Dave.dec.m$DMTRSpeciesCode))
surv[Dave.dec.m$dbh2==0] <- 0
Dave.dec.s <- data.frame(Dave.dec.m,surv=surv,year=Dave.dec.m$year2- Dave.dec.m$year1)
split.ns <- split(Dave.dec.s,  Dave.dec.s$DMTRSpeciesCode)


vec.plot.pbb <- names.pbb [!(names.pbb %in% names.pbb2)]

Dave.dec.s2 <- Dave.dec.s[!(Dave.dec.s$pol.idi.y %in% plot.pb),]
Dave.dec.s2 <- Dave.dec.s2[!(Dave.dec.s2$pol.idi.y %in% vec.plot.pbb),]
## Exclude plot with to high mortalty rate (100%) proably due to e disturbance plots exluded (6 plots excluded).
plot.pb ## list of plot to exclude

split.ns2 <- split(Dave.dec.s2,  Dave.dec.s2$DMTRSpeciesCode)



for (i in 1:7)

{
#windows()
data.temp <- split.ns2[[i]][(split.ns2[[i]]$dbh1/10)<10 & (split.ns2[[i]]$dbh1/10)>5,]
print(  names(split.ns2)[i])
surv<-( data.temp$surv)
print(length(surv))
#print(length(table(factor(split.ns2[[i]]$PLHDPlotName.x))))
#print(range(data.temp$PLHDAltitude))
}


##### annual mortality rate for the small tree

list.Adult.S.0.inf10 <- vector("list",7)
for (i in 1:7)

{
#windows()
data.temp <- split.ns2[[i]][(split.ns2[[i]]$dbh1/10)<10 & (split.ns2[[i]]$dbh1/10)>5,]

surv<-( data.temp$surv)
diam2 <- data.temp$dbh1/10
year <- data.temp$year
datatemp <- data.frame(surv=surv,diam2=diam2,year=year)

par <- list(max.s = .9,year="year")
par_lo <- list(max.s =0.000)
par_hi <- list(max.s = 1)
par_step <- list(max.s=1)
par$x <- "surv"
## Mean in normal PDF
par$prob <- "predicted"
par$size <- 1
## Have it calculate log likelihood
par$log <- TRUE
list.Adult.S.0.inf10[[i]] <- neighanneal (model = model.s0, par,
 source_data= datatemp, par_lo, par_hi, par_step, pdf= dbinom, dep_var="surv", max_iter= 20000)
print(i)
}

list.Adult.S.0.inf10b <- vector("list",7)
for (i in 1:7)

{
#windows()
data.temp <- split.ns[[i]][(split.ns[[i]]$dbh1/10)<10 & (split.ns[[i]]$dbh1/10)>5,]

surv<-( data.temp$surv)
diam2 <- data.temp$dbh1/10
year <- data.temp$year
datatemp <- data.frame(surv=surv,diam2=diam2,year=year)

par <- list(max.s = .9,year="year")
par_lo <- list(max.s =0.000)
par_hi <- list(max.s = 1)
par_step <- list(max.s=1)
par$x <- "surv"
## Mean in normal PDF
par$prob <- "predicted"
par$size <- 1
## Have it calculate log likelihood
par$log <- TRUE
list.Adult.S.0.inf10b[[i]] <- neighanneal (model = model.s0, par,
 source_data= datatemp, par_lo, par_hi, par_step, pdf= dbinom, dep_var="surv", max_iter= 50000)
print(i)
}


Matrix.Param.S.SAP <- matrix(0,nrow=7,ncol=3)
rownames(Matrix.Param.S.SAP) <- names(split.ns)
colnames(Matrix.Param.S.SAP) <-  c("max.s","R²","Repl.")

for (i in 1:7)

{
 Matrix.Param.S.SAP[i,1] <-  1- list.Adult.S.0.inf10[[i]]$best_pars[[1]]
 Matrix.Param.S.SAP[i,2] <- list.Adult.S.0.inf10[[i]]$R2
 Matrix.Param.S.SAP[i,3] <- length(list.Adult.S.0.inf10[[i]]$source_data[,1])
}


Matrix.Param.S.SAPb <- matrix(0,nrow=7,ncol=3)
rownames(Matrix.Param.S.SAPb) <- names(split.ns)
colnames(Matrix.Param.S.SAPb) <-  c("max.s","R²","Repl.")

for (i in 1:7)

{
 Matrix.Param.S.SAPb[i,1] <-  1- list.Adult.S.0.inf10b[[i]]$best_pars[[1]]
 Matrix.Param.S.SAPb[i,2] <- list.Adult.S.0.inf10b[[i]]$R2
 Matrix.Param.S.SAPb[i,3] <- length(list.Adult.S.0.inf10b[[i]]$source_data[,1])
}


